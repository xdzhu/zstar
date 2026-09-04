"""Run the matched Cartesian central control without site-specific paths."""

import argparse
from pathlib import Path
import shutil
import subprocess
import sys

from zstar.shared_abacus import MANIFEST, load_manifest, prepare_shared_abacus
from zstar.shared_response import DEFAULT_DISTANCE


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', choices=['SiC', 't_HfO2', 'alpha_In2Se3', 'hBN', 'MoS2', 'H2O', 'CH4'])
    parser.add_argument('--work', type=Path)
    parser.add_argument('--half-step', action='store_true')
    parser.add_argument('--mp-density', type=float, default=.08)
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    case = Path(__file__).resolve().parent/args.case
    work = args.work.resolve() if args.work else case/('work-cartesian-half' if args.half_step else 'work-cartesian')
    dim = 0 if args.case in ('H2O', 'CH4') else 2 if args.case in ('alpha_In2Se3', 'hBN', 'MoS2') else 3
    if not work.exists():
        shutil.copytree(case/'run', work)
    if not (work/MANIFEST).is_file():
        prepare_shared_abacus(work/'STRU', root=work, scf_input=work/'INPUT', dimension=dim,
                              method='central', displacement_scheme='cartesian-control',
                              displacement_angstrom=DEFAULT_DISTANCE/(2 if args.half_step else 1))
    else:
        meta = load_manifest(work)
        distance = DEFAULT_DISTANCE/(2 if args.half_step else 1)
        if meta['displacement_scheme'] != 'cartesian-control' or meta['dimension'] != dim or abs(meta['nominal_distance_A']-distance) > 1e-12:
            raise ValueError('Existing control uses different settings; choose a fresh --work directory')
    if args.prepare_only:
        print(work)
        return
    commands = [['bec', 'run', '--mp-density', str(args.mp_density)], ['bec', 'stat'], ['bec', 'post']]
    if dim != 0:
        commands.append(['dielectric', 'static', '--dim', str(dim)])
    for command in commands:
        subprocess.run([sys.executable, '-m', 'zstar', *command], cwd=work, check=True)
    if dim == 0:
        script = Path(__file__).resolve().parents[2] / 'tools/shared_response/molecular_validation.py'
        subprocess.run([sys.executable, str(script), str(work), '--output',
                        str(work/'molecular_internal_response.json')], check=True)


if __name__ == '__main__':
    main()
