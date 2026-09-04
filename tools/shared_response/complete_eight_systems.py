"""Complete missing matched Cartesian/Unified benchmarks in isolated directories.

This is a research driver, not a new public command. No source calculation
is modified; the two routes use identical seeds and execution profiles.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import time

import numpy as np
import spglib

from zstar import workflow
from zstar.shared_abacus import collect_shared_abacus, prepare_shared_abacus
from zstar.shared_response import DEFAULT_DISTANCE, read_structure, write_structure


CASES = {'hBN': 2, 'MoS2': 2, 'H2O': 0, 'CH4': 0}
ABACUS = '/home/zhuxd/Software/abacus/INSTALL/3.10.0-LTS/bin/abacus'


def save(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')


def prepare(root, case, source):
    from zstar.abacus_assets import prepare_stru_assets

    output = root / case
    if output.exists():
        raise FileExistsError(f'{output} exists; run it to resume, or use a new root')
    seed = output / 'seed'
    shutil.copytree(source, seed, ignore=shutil.ignore_patterns(
        'OUT.*', 'pyatb*', '*.cube', '*.csr', '*.log', 'results', 'work*'))
    for asset in list((seed / 'assets').rglob('*')) if (seed / 'assets').is_dir() else []:
        if asset.is_file():
            destination = seed / asset.name
            if destination.exists() and destination.read_bytes() != asset.read_bytes():
                raise ValueError(f'Ambiguous asset basename: {asset.name}')
            shutil.copy2(asset, destination)
    assets = prepare_stru_assets(seed / 'STRU', pp_dir=seed, orb_dir=seed,
                                output_dir=seed / '.zstar-assets')
    if assets.changed:
        shutil.copy2(assets.path, seed / 'STRU')
    for asset in assets.assets:
        if asset.resolve() != (seed / asset.name).resolve():
            shutil.copy2(asset, seed / asset.name)
    original = (seed / 'INPUT').read_text()
    text = original
    settings = {'calculation': 'scf', 'cal_force': '1', 'cal_stress': '0',
                'gamma_only': '0', 'scf_thr': '1e-8', 'symmetry': '1',
                'pseudo_dir': '.', 'orbital_dir': '.', 'init_chg': 'auto'}
    if CASES[case] == 0:
        settings['kspacing'] = '0'
        (seed / 'KPT').write_text('K_POINTS\n0\nGamma\n1 1 1 0 0 0\n')
    for key, value in settings.items():
        text = workflow._set_abacus_parameter(text, key, value)
    (seed / 'INPUT.source').write_text(original)
    (seed / 'INPUT').write_text(text)
    atoms = read_structure(seed / 'STRU')
    write_structure(seed / 'STRU', seed / 'STRU.direct', atoms)
    shutil.copy2(seed / 'STRU.direct', seed / 'STRU')
    ds = spglib.get_symmetry_dataset(atoms.totuple(), symprec=1e-5)
    plans = {}
    for scheme in ('unified', 'cartesian'):
        plans[scheme] = prepare_shared_abacus(
            seed / 'STRU', root=output / scheme, scf_input=seed / 'INPUT',
            dimension=CASES[case], symprec=1e-5,
            method='auto' if scheme == 'unified' else 'central',
            displacement_scheme='phonopy' if scheme == 'unified' else 'cartesian-control',
            displacement_angstrom=DEFAULT_DISTANCE,
        )
    plan = {'case': case, 'dimension': CASES[case], 'source': str(source),
            'protocol': 'Cartesian central joint-response control versus Phonopy-auto unified response',
            'scope': 'BEC/APT and Gamma Hessian; not Raman or full dispersion',
            'space_group_of_periodic_cell': ds.international,
            'symprec_angstrom': 1e-5, 'changes_to_source_input': settings,
            'structure': {'cell_angstrom': atoms.cell.tolist(),
                          'positions_angstrom': atoms.positions.tolist(), 'symbols': atoms.symbols},
            'stages': {name: len(p['stages']) for name, p in plans.items()},
            'seed_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in seed.iterdir() if p.is_file()},
            'polarization_mp_density': .08, 'allocated_cores': 40,
            'mpi_omp': {'ABACUS': [1, 40], 'PYATB': [40, 1]},
            'timing': 'monotonic seconds times allocated cores; successful solver calls only',
            'completed': False}
    save(output / 'plan.json', plan)
    print(json.dumps({'case': case, 'space_group': ds.international, 'stages': plan['stages']}))


def run(root, case, schemes):
    from zstar.pyatb_precision import precision_command

    case_root = root / case
    plan = json.loads((case_root / 'plan.json').read_text())
    binary = Path(sys.executable).parent
    abacus = f'mpirun -np 1 {ABACUS}'
    original = workflow._run_shell
    for scheme in schemes:
        output = case_root / scheme
        lock = output / '.worker.lock'
        with lock.open('x') as handle:
            handle.write(f'{socket.gethostname()} {os.getpid()}\n')
        worker = {'host': socket.gethostname(), 'pid': os.getpid(), 'status': 'running',
                  'python': sys.executable, 'abacus': ABACUS,
                  'cpu': subprocess.run(['lscpu'], capture_output=True, text=True).stdout}
        save(output / 'worker.json', worker)

        def timed(command, **kwargs):
            is_pyatb = '-m zstar.pyatb_precision' in command
            solver = command == abacus or is_pyatb
            kwargs['env'] = {**kwargs['env'], 'OMP_NUM_THREADS': '1' if is_pyatb else '40',
                             'MKL_NUM_THREADS': '1' if is_pyatb else '40',
                             'OPENBLAS_NUM_THREADS': '1', 'I_MPI_PIN_DOMAIN': 'omp'}
            record = {'command': command, 'cwd': str(kwargs['cwd']),
                      'host': socket.gethostname(), 'mpi': 40 if is_pyatb else 1,
                      'omp': 1 if is_pyatb else 40, 'allocated_cores': 40,
                      'kind': 'ABACUS' if command == abacus else 'PYATB' if is_pyatb else 'preparation',
                      'include_in_solver_core_hours': solver}
            start = time.monotonic()
            try:
                result = original(command, **kwargs)
                record['success'] = True
                return result
            except Exception as exc:
                record.update(success=False, error=repr(exc))
                raise
            finally:
                record['wall_seconds'] = time.monotonic() - start
                record['allocated_core_hours'] = record['wall_seconds'] * 40 / 3600
                with (output / 'component_times.jsonl').open('a') as handle:
                    handle.write(json.dumps(record) + '\n')
                print(json.dumps(record), flush=True)

        workflow._run_shell = timed
        try:
            workflow.run_serial_workflow(
                output, abacus_command=abacus,
                pyatb_input=str(binary / 'pyatb_input'), pyatb_executable=str(binary / 'pyatb'),
                pyatb_command=precision_command(f'mpirun -np 40 {binary / "pyatb"}'),
                omp_threads=40, dimensionality=plan['dimension'],
                mp_density=plan['polarization_mp_density'],
            )
            collect_shared_abacus(output)
            worker['status'] = 'completed'
            save(output / 'completed.json', {'host': socket.gethostname(), 'status': 'completed'})
        except Exception as exc:
            worker.update(status='failed', error=repr(exc))
            raise
        finally:
            save(output / 'worker.json', worker)
            workflow._run_shell = original
            lock.unlink()
    plan['completed'] = all((case_root / s / 'completed.json').is_file() for s in ('unified', 'cartesian'))
    save(case_root / 'plan.json', plan)


def launch(root, case):
    case_root = root / case
    marker = case_root / 'launch.json'
    if marker.exists():
        raise FileExistsError(f'Launch already recorded at {marker}; inspect before restarting')
    script = Path(__file__).with_name('eight_system_worker.sh')
    with (case_root / 'driver.log').open('ab') as log:
        process = subprocess.Popen(['bash', str(script), str(root), case],
                                   stdin=subprocess.DEVNULL, stdout=log,
                                   stderr=subprocess.STDOUT, start_new_session=True)
    data = {'host': socket.gethostname(), 'pid': process.pid, 'case': case,
            'worker_script': str(script), 'status': 'launched_not_yet_validated'}
    save(marker, data)
    print(json.dumps(data))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['prepare', 'run', 'launch'])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--case', choices=list(CASES), required=True)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--scheme', choices=['unified', 'cartesian', 'both'], default='both')
    args = parser.parse_args()
    if args.action == 'prepare':
        if args.source is None:
            parser.error('--source is required for prepare')
        prepare(args.root.resolve(), args.case, args.source.resolve())
    elif args.action == 'run':
        run(args.root.resolve(), args.case, ['unified', 'cartesian'] if args.scheme == 'both' else [args.scheme])
    else:
        launch(args.root.resolve(), args.case)
