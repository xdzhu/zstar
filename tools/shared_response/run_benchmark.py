"""Run matched shared/Cartesian controls and record measured component costs."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import time

import numpy as np
import spglib

from zstar import workflow
from zstar.shared_abacus import prepare_shared_abacus, collect_shared_abacus
from zstar.shared_response import DEFAULT_DISTANCE, read_structure

ROOT = Path('/home/zhuxd/abacus/agent-runs/20260904-shared-response-benchmark')
BIN = Path('/home/zhuxd/Software/anaconda3/envs/icu_copy/bin')


def run(case, scheme):
    if case == 'sic':
        source = Path('/home/zhuxd/abacus/agent-runs/20260903-spectra-backend-benchmark/sic_abacus/workflow/0.no-move')
        structure, dim = source / 'STRU', 3
    elif case == 'hfo2':
        source = Path('/home/zhuxd/abacus/zstar_validation/hfo2_pbesol_tzdp9_20260901/bec/0.no-move')
        structure, dim = source / 'STRU', 3
    else:
        source = ROOT / 'in2se3_nc2017/relax_symmetry_verified'
        log = source / 'OUT.SHARED_IN2SE3/running_relax.log'
        if 'Relaxation is converged!' not in log.read_text():
            raise RuntimeError('In2Se3 relaxation is not converged')
        structure, dim = source / 'OUT.SHARED_IN2SE3/STRU_ION_D', 2
    output = ROOT / case / scheme
    output.mkdir(parents=True, exist_ok=True)
    if not (output / 'shared_response.json').exists():
        atoms = read_structure(structure)
        dataset = spglib.get_symmetry_dataset(atoms.totuple(), symprec=1e-5)
        provenance = {'source': str(source), 'structure': str(structure),
                      'structure_sha256': hashlib.sha256(structure.read_bytes()).hexdigest(),
                      'space_group': dataset.international, 'symprec_A': 1e-5,
                      'cell_A': atoms.cell.tolist(), 'symbols': atoms.symbols,
                      'host': socket.gethostname(), 'control': scheme,
                      'comparison': 'same-inputs force-and-polarization SCFs; Cartesian central vs Phonopy auto'}
        (output / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
        shutil.copy2(source / 'KPT', output / 'KPT')
        extra = []
        # Archive STRU may contain relative filenames resolved at its old stage.
        from zstar.gen_polar import _abacus_assets_from_stru
        extra.extend(str(p) for p in _abacus_assets_from_stru(structure, source_dir=source))
        prepare_shared_abacus(structure, root=output, scf_input=source / 'INPUT',
                              dimension=dim, symprec=1e-5, input_sets=extra,
                              method='auto' if scheme.startswith('shared') else 'central',
                              displacement_angstrom=DEFAULT_DISTANCE / 2 if scheme.endswith('-half') else DEFAULT_DISTANCE,
                              displacement_scheme='phonopy' if scheme.startswith('shared') else 'cartesian-control')
    original = workflow._run_shell
    timing = output / 'component_times.jsonl'

    def timed(command, **kwargs):
        ranks = int(os.environ.get('ZSTAR_PYATB_MPI', '1'))
        is_pyatb = command.startswith('mpirun -np 1 ') and ('pyatb' in command or 'zstar.pyatb_precision' in command)
        if is_pyatb and ranks != 1:
            if ranks < 1 or 40 % ranks:
                raise ValueError('PYATB MPI count must divide the 40-core allocation')
            command = command.replace('mpirun -np 1 ', f'mpirun -np {ranks} ', 1)
            kwargs['env'] = {**kwargs['env'], 'OMP_NUM_THREADS': str(40//ranks),
                             'MKL_NUM_THREADS': str(40//ranks), 'OPENBLAS_NUM_THREADS': '1'}
        start = time.monotonic()
        record = {'command': command, 'cwd': str(kwargs['cwd']),
                  'host': socket.gethostname(), 'omp': 40//ranks if is_pyatb else 40,
                  'mpi': ranks if is_pyatb else 1, 'allocated_cores': 40 if command.startswith('mpirun') else 1}
        try:
            result = original(command, **kwargs)
            record['success'] = True
            return result
        except Exception as exc:
            record['success'] = False
            record['error'] = repr(exc)
            raise
        finally:
            record['wall_seconds'] = time.monotonic() - start
            record['allocated_core_hours'] = record['wall_seconds'] * record['allocated_cores'] / 3600
            with timing.open('a') as handle:
                handle.write(json.dumps(record) + '\n')

    workflow._run_shell = timed
    os.environ['MKL_NUM_THREADS'] = '40'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    try:
        workflow.run_serial_workflow(output,
            abacus_command='mpirun -np 1 /home/zhuxd/Software/abacus/INSTALL/3.10.0-LTS/bin/abacus',
            pyatb_input=str(BIN / 'pyatb_input'), pyatb_executable=str(BIN / 'pyatb'),
            pyatb_command=f'mpirun -np 1 {BIN / "pyatb"}', omp_threads=40,
            dimensionality=dim, mp_density=.08)
        collect_shared_abacus(output)
    finally:
        workflow._run_shell = original


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('case', choices=['sic', 'hfo2', 'in2se3'])
    parser.add_argument('scheme', choices=['shared', 'shared-half', 'cartesian', 'cartesian-half', 'both'])
    args = parser.parse_args()
    for scheme in (['shared', 'cartesian'] if args.scheme == 'both' else [args.scheme]):
        run(args.case, scheme)
