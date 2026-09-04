"""Recover retained molecular APT evidence without changing source calculations."""

import argparse
import hashlib
import json
from pathlib import Path
import tarfile


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    checksums = {}
    with tarfile.open(args.output, 'w:gz') as archive:
        for case in ['H2O', 'CH4']:
            root = args.root / case
            selected = set()
            for pattern in ['*.json', 'Z-BORN*.out', '0.no-move/STRU', '0.no-move/KPT',
                            '**/INPUT', '**/running_scf.log', '**/time.json']:
                selected.update(p for p in root.glob(pattern) if p.is_file())
            for path in sorted(selected):
                name = (Path(case) / path.relative_to(root)).as_posix()
                archive.add(path, arcname=name)
                checksums[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    args.output.with_suffix('.sha256.json').write_text(json.dumps(checksums, indent=2) + '\n')
    print(f'Exported {len(checksums)} inputs, logs and response records; cubes remain at source.')
