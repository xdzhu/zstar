"""Test Berry-mesh symmetry leakage without repeating or replacing any SCF."""

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import socket
import subprocess
import time

import numpy as np

from zstar.deal_polar import _parse_pyatb_polar_file, _read_pyatb_geom
from zstar.pyatb_compat import _OPTICAL_BLOCK_RE
from zstar.shared_response import make_phonopy, read_structure, symmetry_operations
from zstar.workflow import prepare_pyatb_assets

ROOT = Path('/home/zhuxd/abacus/agent-runs/20260904-shared-response-benchmark')
PYTHON = '/home/zhuxd/Software/anaconda3/envs/icu_copy/bin/python'


def run(nk, atom, mpi=1):
    if mpi < 1 or 40 % mpi:
        raise ValueError('MPI ranks must divide the 40-core allocation')
    output = ROOT/'in2se3'/f'berry-mesh-{nk}'/f'atom-{atom}'
    output.mkdir(parents=True, exist_ok=True)
    data = {'mesh': [nk, nk, 2], 'atom_zero_based': atom,
            'purpose': 'polarization-only convergence; original SCF and responses unchanged',
            'host': socket.gethostname(), 'schemes': {}}
    for scheme in ('shared', 'cartesian'):
        root = ROOT/'in2se3'/scheme
        manifest = json.loads((root/'shared_response.json').read_text())
        observations = [s for s in manifest['stages'] if s['atom'] == atom]
        responses = {}
        for name in ['0.no-move']+[s['name'] for s in observations]:
            source = root/name
            dest = output/scheme/name
            dest.mkdir(parents=True, exist_ok=True)
            marker = dest/'timing.json'
            if not marker.is_file():
                text = _OPTICAL_BLOCK_RE.sub('', (source/'pyatb/Input').read_text())
                for axis in (1, 2):
                    text, count = re.subn(rf'(?m)^(\s*nk{axis}\s+)\d+', rf'\g<1>{nk}', text)
                    if count != 1:
                        raise ValueError('Expected exactly one polarization mesh entry')
                (dest/'Input').write_text(text)
                shutil.copy2(source/'STRU', dest/'STRU')
                prepare_pyatb_assets(source, dest)
                env = {**os.environ, 'OMP_NUM_THREADS': str(40//mpi), 'MKL_NUM_THREADS': str(40//mpi),
                       'OPENBLAS_NUM_THREADS': '1'}
                start = time.monotonic()
                with (dest/'run.log').open('w') as log:
                    subprocess.run(['mpirun', '-np', str(mpi), PYTHON, '-m', 'zstar.pyatb_precision'],
                                   cwd=dest, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
                elapsed = time.monotonic()-start
                marker.write_text(json.dumps({'wall_seconds': elapsed, 'reserved_core_hours': elapsed*40/3600,
                                             'mpi': mpi, 'omp': 40//mpi})+'\n')
            values = _parse_pyatb_polar_file(dest/'Out/Polarization/polarization.dat')
            responses[name] = np.array(values[:3]), np.array(values[3:])
            print(nk, scheme, name, 'complete', flush=True)
        transform, volume = _read_pyatb_geom(output/scheme/'0.no-move/Out/input.json')
        p0, q = responses['0.no-move']
        x, y = [], []
        for s in observations:
            delta = responses[s['name']][0]-p0
            delta -= np.rint(delta/q)*q
            x.append(s['displacement_A'])
            y.append(delta@transform*volume/1.602176634e-19/1e-10)
        atoms = read_structure(root/'0.no-move/STRU')
        ops = symmetry_operations(make_phonopy(atoms, symprec=manifest['symprec_A']), dimension=2)
        expanded_x, expanded_y = [], []
        for r, perm in ops:
            if perm[atom] == atom:
                expanded_x.extend(np.asarray(x)@r.T)
                expanded_y.extend(np.asarray(y)@r.T)
        constrained = np.linalg.lstsq(expanded_x, expanded_y, rcond=None)[0].T
        record = {'observations': [{'displacement_A': u, 'berry_dipole_change_e_A': v.tolist()} for u, v in zip(x, y)],
                  'site_symmetric_Berry_BEC': constrained.tolist()}
        if scheme == 'cartesian':
            record['unconstrained_Berry_BEC'] = np.linalg.lstsq(x, y, rcond=None)[0].T.tolist()
        data['schemes'][scheme] = record
    a, b = [np.array(data['schemes'][s]['site_symmetric_Berry_BEC']) for s in ('shared', 'cartesian')]
    data['max_in_plane_BEC_difference_e'] = float(np.max(np.abs(a[:2]-b[:2])))
    data['reserved_core_hours'] = sum(json.loads(p.read_text())['reserved_core_hours'] for p in output.glob('*/*/timing.json'))
    (output/'comparison.json').write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps(data, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('nk', type=int)
    parser.add_argument('--atom', type=int, default=3)
    parser.add_argument('--mpi', type=int, default=1)
    args = parser.parse_args()
    run(args.nk, args.atom, args.mpi)
