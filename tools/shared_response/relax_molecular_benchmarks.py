"""Relax the benchmark seed with its actual basis before matched response runs."""

import argparse
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import time

import numpy as np
from phonopy.interface.abacus import read_abacus_output

from complete_eight_systems import ABACUS, prepare, run
from zstar.workflow import _set_abacus_parameter


def calculate(root, source, case):
    relax = root / 'relaxation' / case
    if relax.exists():
        raise FileExistsError(f'Inspect existing relaxation before restarting: {relax}')
    shutil.copytree(source / case / 'seed', relax,
                    ignore=shutil.ignore_patterns('.zstar-assets', 'INPUT.source', 'STRU.direct'))
    text = (relax / 'INPUT').read_text()
    for key, value in {'calculation': 'relax', 'relax_nmax': '80', 'force_thr_ev': '0.003',
                       'out_mat_hs2': '0', 'out_mat_r': '0', 'out_chg': '0', 'init_chg': 'auto'}.items():
        text = _set_abacus_parameter(text, key, value)
    (relax / 'INPUT').write_text(text)
    env = {**os.environ, 'OMP_NUM_THREADS': '40', 'MKL_NUM_THREADS': '40',
           'OPENBLAS_NUM_THREADS': '1', 'I_MPI_PIN_DOMAIN': 'omp'}
    start = time.monotonic()
    with (relax / 'abacus.log').open('w') as log:
        subprocess.run(['mpirun', '-np', '1', ABACUS], cwd=relax, env=env,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    logs = list(relax.glob('OUT.*/running_relax.log'))
    structures = list(relax.glob('OUT.*/STRU_ION_D'))
    if len(logs) != 1 or len(structures) != 1:
        raise ValueError('Missing unique relaxation log or final STRU_ION_D')
    forces = np.asarray(read_abacus_output(str(logs[0])))
    if not forces.size or np.max(np.linalg.norm(forces, axis=1)) > .005:
        raise ValueError(f'Relaxation not sufficiently converged: {forces}')
    summary = {'case': case, 'host': socket.gethostname(),
               'maximum_force_eV_A': float(np.max(np.linalg.norm(forces, axis=1))),
               'wall_seconds': time.monotonic() - start,
               'included_in_response_benchmark_cost': False,
               'source_seed': str(source / case / 'seed')}
    (relax / 'relaxation.json').write_text(json.dumps(summary, indent=2)+'\n')
    seed = root / 'relaxed-seeds' / case
    shutil.copytree(source / case / 'seed', seed,
                    ignore=shutil.ignore_patterns('.zstar-assets', 'INPUT.source', 'STRU.direct'))
    shutil.copy2(structures[0], seed / 'STRU')
    prepare(root, case, seed)
    run(root, case, ['unified', 'cartesian'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['launch', 'run'])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--case', choices=['H2O', 'CH4'], required=True)
    args = parser.parse_args()
    if args.action == 'run':
        calculate(args.root, args.source, args.case)
    else:
        args.root.mkdir(parents=True, exist_ok=True)
        marker = args.root / f'{args.case}.launch.json'
        if marker.exists():
            raise FileExistsError(marker)
        command = [sys.executable, '-u', str(Path(__file__).resolve()), 'run', '--root', str(args.root),
                   '--source', str(args.source), '--case', args.case]
        with (args.root / f'{args.case}.log').open('ab') as log:
            process = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=log,
                                       stderr=subprocess.STDOUT, start_new_session=True)
        marker.write_text(json.dumps({'host': socket.gethostname(), 'pid': process.pid, 'command': command})+'\n')
        print(marker.read_text())
