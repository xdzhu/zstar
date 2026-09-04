"""Execute isolated shell/Slurm/Torque header smoke tests without a scheduler."""

import json
from pathlib import Path
import subprocess
import tempfile

from zstar.workflow import generate_backend_script
from zstar.spectra_frontend import _write_driver


def run():
    results = {}
    with tempfile.TemporaryDirectory(prefix='zstar-header-smoke-') as folder:
        for system, directive in [('shell', '# local'), ('slurm', '#SBATCH --ntasks=1'),
                                  ('torque', '#PBS -l nodes=1:ppn=1')]:
            root = Path(folder) / system
            (root / '0.no-move').mkdir(parents=True)
            (root / '0.no-move/INPUT-scf').write_text('INPUT_PARAMETERS\n')
            header = root / 'header.sh'
            header.write_text('#!/usr/bin/env bash\n'+directive+'\n'
                              'export ZSTAR_HEADER_PROBE=ready\n'
                              'zstar() { test "$ZSTAR_HEADER_PROBE" = ready; printf "PROBE %s\\n" "$*"; }\n')
            script = generate_backend_script(root, backend=system, header_file=header, tasks=1, cpus_per_task=1)
            subprocess.run(['bash', '-n', str(script)], check=True)
            execution = subprocess.run(['bash', str(script)], check=True, capture_output=True, text=True)
            assert 'PROBE workflow run' in execution.stdout
            # A setup failure must prevent the mocked calculator command.
            header.write_text(directive+'\nexit 9\n')
            failed = generate_backend_script(root, backend=system, header_file=header, tasks=1, cpus_per_task=1)
            stopped = subprocess.run(['bash', str(failed)], capture_output=True, text=True)
            assert stopped.returncode == 9 and 'PROBE' not in stopped.stdout
            header.write_text(directive+'\nexport ZSTAR_HEADER_PROBE=ready\n'
                              'zstar() { test "$ZSTAR_HEADER_PROBE" = ready; printf "PROBE %s\\n" "$*"; }\n')
            spectra = {}
            for calculator in ['abacus', 'vasp', 'cp2k', 'qe']:
                script = _write_driver(str(root), system=system, output=None,
                    job_name='probe', nodes=1, tasks=1, cpus_per_task=1, walltime='00:01:00',
                    queue=None, account=None, env_script=None, dry_run=True,
                    header_file=str(header), calculator=calculator)
                subprocess.run(['bash', '-n', str(script)], check=True)
                execution = subprocess.run(['bash', str(script)], check=True, capture_output=True, text=True)
                assert 'PROBE spectra run' in execution.stdout and 'PROBE spectra post' not in execution.stdout
                spectra[calculator] = 'syntax and mocked execution passed'
            results[system] = {'bash_syntax': 'passed', 'embedded_environment': 'passed',
                               'setup_failure_stops_execution': 'passed',
                               'scheduler_submission': 'not performed', 'DFT_execution': 'mocked',
                               'canonical_spectra': spectra}
    print(json.dumps(results, indent=2))
    return results


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    results = run()
    if args.output:
        args.output.write_text(json.dumps(results, indent=2)+'\n')
