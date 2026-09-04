"""Repeat only PYATB polarization with its lossless writer; preserve old files."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from zstar.pyatb_compat import _OPTICAL_BLOCK_RE
from zstar.workflow import prepare_pyatb_assets
from zstar.shared_abacus import collect_shared_abacus

ROOT = Path('/home/zhuxd/abacus/agent-runs/20260904-shared-response-benchmark')
PYTHON = '/home/zhuxd/Software/anaconda3/envs/icu_copy/bin/python'


def refine(case, scheme):
    root = ROOT/case/scheme
    meta = json.loads((root/'shared_response.json').read_text())
    for name in ['0.no-move'] + [s['name'] for s in meta['stages']]:
        stage = root/name
        old = stage/'pyatb/Out/Polarization'
        dest = stage/'pyatb-precision'
        marker = dest/'precision_timing.json'
        if not marker.is_file():
            dest.mkdir(exist_ok=True)
            source_input = stage/'pyatb/Input'
            if not source_input.is_file():
                raise FileNotFoundError(source_input)
            text = _OPTICAL_BLOCK_RE.sub('', source_input.read_text())
            (dest/'Input').write_text(text)
            shutil.copy2(stage/'STRU', dest/'STRU')
            prepare_pyatb_assets(stage, dest)
            env = os.environ.copy()
            env['OMP_NUM_THREADS'] = env['MKL_NUM_THREADS'] = '40'
            env['OPENBLAS_NUM_THREADS'] = '1'
            start = time.monotonic()
            with (dest/'run.log').open('w') as log:
                subprocess.run(['mpirun', '-np', '1', PYTHON, '-m', 'zstar.pyatb_precision'],
                               cwd=dest, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            elapsed = time.monotonic()-start
            marker.write_text(json.dumps({'wall_seconds': elapsed, 'reserved_core_hours': elapsed*40/3600,
                                         'purpose': 'output-precision diagnostic; no new SCF',
                                         'original_input_sha256': hashlib.sha256(source_input.read_bytes()).hexdigest()}, indent=2)+'\n')
        precise = dest/'Out/Polarization'
        if not (precise/'zstar_precision.json').is_file():
            raise RuntimeError(f'Missing full precision writer record: {dest}')
        backup = stage/'pyatb/Out/Polarization.before_precision'
        if not backup.exists():
            shutil.copytree(old, backup)
        for filename in ('polarization.dat', 'polarization.rounded.dat', 'zstar_precision.json'):
            shutil.copy2(precise/filename, old/filename)
        print(case, scheme, name, 'precision collected', flush=True)
    for filename in ('shared_response_result.json', 'Z-BORN-all.out', 'Z-BORN-symm.out', 'BORN'):
        path = root/filename
        if path.is_file() and not path.with_suffix(path.suffix+'.before_precision').exists():
            shutil.copy2(path, path.with_suffix(path.suffix+'.before_precision'))
    collect_shared_abacus(root)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('case', choices=['sic', 'in2se3', 'hfo2', 'remaining-small'])
    parser.add_argument('scheme', choices=['shared', 'cartesian', 'shared-half', 'both'])
    args = parser.parse_args()
    if args.case == 'remaining-small':
        for case, scheme in [('in2se3', 'shared-half'), ('sic', 'shared'), ('sic', 'cartesian')]:
            refine(case, scheme)
    else:
        for scheme in (['shared', 'cartesian'] if args.scheme == 'both' else [args.scheme]):
            refine(args.case, scheme)
