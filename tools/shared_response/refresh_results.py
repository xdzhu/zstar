"""Refresh completed derived records with the current code; no solver launch."""

import argparse
from pathlib import Path

from zstar.shared_abacus import collect_shared_abacus


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, default=Path('/home/zhuxd/abacus/agent-runs/20260904-shared-response-benchmark'))
    args = parser.parse_args()
    for case in ('sic', 'hfo2', 'in2se3'):
        for scheme in ('shared', 'cartesian', 'shared-half', 'cartesian-half'):
            path = args.root/case/scheme
            if (path/'shared_response_result.json').is_file():
                collect_shared_abacus(path)
