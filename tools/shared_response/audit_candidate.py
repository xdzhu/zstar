"""Check portable candidate inputs and freeze additive archive checksums."""

import argparse
import json
from pathlib import Path
import tempfile

from zstar.abacus_assets import prepare_stru_assets, sha256_file


def audit(root, write_checksums=False):
    records = json.loads((root/'examples/manifest.json').read_text(encoding='utf-8'))['cases']
    selected = [r for r in records if r['id'].startswith('unified_')
                or r['id'] == 'in2se3_pbesol_bec']
    result = {}
    for record in selected:
        case = root/'examples'/record['path']
        for name in ('README.md', 'README.zh-CN.md', 'run.sh', 'run/STRU', 'run/INPUT'):
            if not (case/name).is_file():
                raise ValueError(f'Missing {case/name}')
        with tempfile.TemporaryDirectory(prefix='zstar-asset-audit-') as temporary:
            assets = case/'run/assets'
            if not assets.is_dir():
                assets = case/'run'
            prepared = prepare_stru_assets(case/'run/STRU', pp_dir=assets,
                orb_dir=assets, output_dir=temporary)
        manifest = case/'checksums.json'
        existing = json.loads(manifest.read_text()) if manifest.is_file() else {}
        for name, digest in existing.items():
            if sha256_file(case/name) != digest:
                raise ValueError(f'Previously frozen evidence changed: {case/name}')
        checksums = {}
        for path in sorted(case.rglob('*')):
            if path.is_symlink():
                raise ValueError(f'Archive must not contain symlinks: {path}')
            if not path.is_file() or '__pycache__' in path.parts or path == manifest:
                continue
            if path.stat().st_size >= 100*1024**2:
                raise ValueError(f'GitHub file-size limit: {path}')
            checksums[path.relative_to(case).as_posix()] = sha256_file(path)
        if write_checksums:
            manifest.write_text(json.dumps(checksums, indent=2)+'\n', encoding='utf-8')
        result[record['id']] = {'status': 'passed', 'portable_assets': len(prepared.assets),
            'previous_checksums_verified': len(existing), 'files': len(checksums)}
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, default=Path.cwd())
    parser.add_argument('--write-checksums', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = audit(args.root, args.write_checksums)
    encoded = json.dumps(report, indent=2)+'\n'
    if args.output:
        args.output.write_text(encoded, encoding='utf-8')
    print(encoded)
