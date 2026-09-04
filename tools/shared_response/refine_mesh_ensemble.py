"""Refine Berry responses in separate, auditable result trees without new SCFs."""

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import socket
import subprocess
import time

from zstar.pyatb_compat import _OPTICAL_BLOCK_RE
from zstar.shared_abacus import collect_shared_abacus, load_manifest
from zstar.workflow import prepare_pyatb_assets

from polarization_mesh_diagnostic import PYTHON, ROOT
CASE = 'in2se3'
MESH = 88


def paths(scheme):
    return ROOT/CASE/scheme, ROOT/CASE/f'{scheme}-mesh{MESH}'


def prepare(scheme):
    source, dest = paths(scheme)
    if (dest/'refinement.json').is_file():
        return
    dest.mkdir(exist_ok=True)
    meta = load_manifest(source)
    for filename in ['shared_response.json', 'STRU', 'provenance.json', *meta['input_hashes']]:
        if filename == 'provenance.json' and not (source/filename).is_file():
            continue
        target = dest/filename
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source/filename, target)
    for name in ['0.no-move']+[s['name'] for s in meta['stages']]:
        stage = dest/name
        for pattern in ('OUT.*/running_scf.log', 'OUT.*/*CHG.cube'):
            for path in (source/name).glob(pattern):
                target = dest/path.relative_to(source)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, target)
        optical = source/name/'pyatb/Out/Optical_Conductivity'
        if optical.is_dir():
            shutil.copytree(optical, stage/'pyatb/Out/Optical_Conductivity', dirs_exist_ok=True)
    shutil.copytree(source/'.zstar', dest/'.zstar', dirs_exist_ok=True)
    (dest/'refinement.json').write_text(json.dumps({
        'source': str(source), 'polarization_mesh': [MESH,MESH,2], 'SCFs_repeated': False,
        'electronic_dielectric': 'unchanged original reference static calculation',
        'timing': 'original workflow and additional polarization-only refinement reported separately',
    }, indent=2)+'\n')


def run(scheme, start, stop):
    source, dest = paths(scheme)
    if not (dest/'refinement.json').is_file():
        raise ValueError('Prepare the refinement before starting workers')
    meta = load_manifest(dest)
    names = ['0.no-move']+[s['name'] for s in meta['stages']]
    for name in names[start:stop]:
        target = dest/name/'pyatb'
        marker = target/'refinement_timing.json'
        if marker.is_file():
            continue
        cached = ROOT/CASE/f'berry-mesh-{MESH}/atom-3'/scheme/name
        if (cached/'timing.json').is_file():
            for filename in ('Input', 'STRU', 'Out/input.json'):
                output = target/filename
                output.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(cached/filename, output)
            shutil.copytree(cached/'Out/Polarization', target/'Out/Polarization', dirs_exist_ok=True)
            record = json.loads((cached/'timing.json').read_text())
            record['reused_diagnostic'] = str(cached)
        else:
            target.mkdir(parents=True, exist_ok=True)
            text = _OPTICAL_BLOCK_RE.sub('', (source/name/'pyatb/Input').read_text())
            for axis in (1,2):
                text, count = re.subn(rf'(?m)^(\s*nk{axis}\s+)\d+', rf'\g<1>{MESH}', text)
                if count != 1:
                    raise ValueError('Missing or ambiguous Berry mesh parameter')
            (target/'Input').write_text(text)
            shutil.copy2(source/name/'STRU', target/'STRU')
            prepare_pyatb_assets(source/name, target)
            env = {**os.environ, 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'}
            begin = time.monotonic()
            with (target/'run.log').open('w') as log:
                subprocess.run(['mpirun', '-np', '40', PYTHON, '-m', 'zstar.pyatb_precision'],
                               cwd=target, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            seconds = time.monotonic()-begin
            record = {'wall_seconds': seconds, 'reserved_core_hours': seconds*40/3600, 'mpi':40, 'omp':1}
        record['worker_host'] = socket.gethostname()
        marker.write_text(json.dumps(record, indent=2)+'\n')
        print(scheme, name, f'mesh{MESH} complete', flush=True)


def collect(scheme):
    source, dest = paths(scheme)
    meta = load_manifest(dest)
    names = ['0.no-move']+[s['name'] for s in meta['stages']]
    records = [json.loads((dest/n/'pyatb/refinement_timing.json').read_text()) for n in names]
    records_path = dest/'refinement_cost.json'
    records_path.write_text(json.dumps({'stages':len(records),
        'additional_reserved_core_hours':sum(r['reserved_core_hours'] for r in records),
        'per_stage':dict(zip(names,records))}, indent=2)+'\n')
    # Keep the original measured run ledger verbatim. Refinement costs have
    # their own ledger instead of pretending that another SCF was performed.
    shutil.copy2(source/'component_times.jsonl', dest/'component_times.jsonl')
    collect_shared_abacus(dest)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['prepare','run','collect','shared-then-control','one-complete','launch'])
    parser.add_argument('--scheme', choices=['shared','unified','cartesian'], default='shared')
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--case', default=CASE)
    parser.add_argument('--mesh', type=int, default=MESH)
    parser.add_argument('--attempt', type=int, default=1)
    parser.add_argument('--start', type=int, default=0)
    parser.add_argument('--stop', type=int)
    args = parser.parse_args()
    ROOT, CASE, MESH = args.root, args.case, args.mesh
    if MESH < 2:
        parser.error('--mesh must be at least 2')
    if args.action == 'prepare':
        for scheme in ('shared','cartesian'):
            prepare(scheme)
    elif args.action == 'run':
        run(args.scheme, args.start, args.stop)
    elif args.action == 'collect':
        collect(args.scheme)
    elif args.action == 'one-complete':
        prepare(args.scheme)
        run(args.scheme, 0, None)
        collect(args.scheme)
    elif args.action == 'launch':
        _, dest = paths(args.scheme)
        log = dest.parent / f'{dest.name}.log'
        marker = log.with_suffix(f'.launch{args.attempt}.json')
        if marker.exists():
            raise FileExistsError(marker)
        command = [PYTHON, '-u', str(Path(__file__).resolve()), 'one-complete', '--root', str(ROOT),
                   '--case', CASE, '--scheme', args.scheme, '--mesh', str(MESH)]
        with log.open('ab') as handle:
            worker = subprocess.Popen(command, stdout=handle, stderr=subprocess.STDOUT,
                                      stdin=subprocess.DEVNULL, start_new_session=True)
        marker.write_text(json.dumps({'host': socket.gethostname(), 'pid': worker.pid, 'command': command})+'\n')
        print(marker.read_text())
    else:
        run('shared', 0, None)
        collect('shared')
        run('cartesian', 16, None)
