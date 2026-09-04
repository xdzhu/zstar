"""Export a completed Berry-only refinement without replacing baseline results."""

import argparse
import json
from pathlib import Path
import tarfile

from export_examples import copy, export_response


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mesh', type=int, required=True)
    parser.add_argument('--tar', type=Path, required=True)
    args = parser.parse_args()
    for scheme in ['unified', 'cartesian']:
        name = f'{scheme}-mesh{args.mesh}'
        root = args.source / name
        if not (root / 'shared_response_result.json').is_file():
            raise ValueError(f'Incomplete refinement: {root}')
        export_response(root, args.output / name)
        meta = json.loads((root / 'shared_response.json').read_text())
        for stage in ['0.no-move'] + [s['name'] for s in meta['stages']]:
            for item in ['STRU', 'INPUT-scf', 'INPUT', 'KPT']:
                copy(root / stage / item, args.output / name / stage / item)
    copy(args.source / f'mesh{args.mesh}-comparison.json',
         args.output / f'mesh{args.mesh}-comparison.json')
    with tarfile.open(args.tar, 'w:gz') as archive:
        for path in sorted(args.output.iterdir()):
            archive.add(path, arcname=path.name)
