"""Audit paired production runs and export independently reprocessable evidence."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tarfile

from export_examples import export_response, copy
from report_benchmark import result, compare


def solver_cost(root):
    records = [json.loads(line) for line in (root / 'component_times.jsonl').read_text().splitlines()]
    groups = {'ABACUS': 0., 'PYATB': 0., 'preparation': 0.}
    calls = {key: 0 for key in groups}
    for record in records:
        if not record.get('success', True):
            raise ValueError(f'Failed timing record needs an explicit production selection: {root}')
        command = record['command']
        kind = ('ABACUS' if command.rstrip().endswith('/abacus') else
                'preparation' if 'pyatb_input' in command else 'PYATB')
        groups[kind] += record['wall_seconds'] * record.get('allocated_cores', 40) / 3600
        calls[kind] += 1
    return {'solver_core_hours': groups['ABACUS'] + groups['PYATB'],
            'component_core_hours': groups, 'calls': calls,
            'definition': 'successful ABACUS + PYATB calls, including reference and band gate; input preparation excluded'}


def audit(root, output, export_cases=None):
    output.mkdir(parents=True, exist_ok=True)
    report = {'cases': {}, 'basis': 'matched seeds and 40-core allocation within each pair'}
    for case in ['hBN', 'MoS2', 'H2O', 'CH4']:
        source = root / case
        if not all((source / scheme / 'completed.json').is_file() for scheme in ['unified', 'cartesian']):
            report['cases'][case] = {'status': 'pending'}
            continue
        a, b = result(source / 'unified'), result(source / 'cartesian')
        for name, data in [('unified', a), ('cartesian', b)]:
            data['solver_cost'] = solver_cost(source / name)
            if data['solver_cost']['calls']['ABACUS'] != data['total_stages']:
                raise ValueError(f'{case}/{name}: SCF count does not match stage count')
        comparison = compare(a, b)
        comparison.pop('measured_shared_over_cartesian_core_hours', None)
        comparison['solver_speedup'] = b['solver_cost']['solver_core_hours'] / a['solver_cost']['solver_core_hours']
        report['cases'][case] = {'status': 'paired_result_audited', 'unified': a, 'cartesian': b,
                                  'comparison': comparison,
                                  'plan': json.loads((source / 'plan.json').read_text())}
        print(case, json.dumps(comparison), flush=True)
        if export_cases and case not in export_cases:
            continue
        dest = output / case
        shutil.copytree(source / 'seed', dest / 'run', dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns('.zstar-assets', 'INPUT.source', 'STRU.direct'))
        copy(source / 'plan.json', dest / 'results/plan.json')
        relaxation = root / 'relaxation' / case
        if (relaxation / 'relaxation.json').is_file():
            copy(relaxation / 'relaxation.json', dest / 'results/relaxation.json')
            for pattern in ['INPUT', 'STRU', 'KPT', 'OUT.*/running_relax.log', 'OUT.*/STRU_ION_D']:
                for path in relaxation.glob(pattern):
                    copy(path, dest / 'results/relaxation' / path.relative_to(relaxation))
        for name in ['unified', 'cartesian']:
            src, target = source / name, dest / 'results' / name
            export_response(src, target)
            for item in ['worker.json', 'completed.json']:
                copy(src / item, target / item)
            meta = json.loads((src / 'shared_response.json').read_text())
            for stage in ['0.no-move'] + [s['name'] for s in meta['stages']]:
                for item in ['STRU', 'INPUT-scf', 'INPUT', 'KPT']:
                    copy(src / stage / item, target / stage / item)
            for path in (src / '.zstar/stages').glob('*.json'):
                copy(path, target / '.zstar/stages' / path.name)
        (dest / 'results/comparison.json').write_text(json.dumps(comparison, indent=2) + '\n')
        hashes = {p.relative_to(dest).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in dest.rglob('*') if p.is_file() and p.name != 'checksums.json'}
        (dest / 'checksums.json').write_text(json.dumps(hashes, indent=2) + '\n')
    (output / 'four_new_benchmarks.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--tar', type=Path)
    parser.add_argument('--export-case', action='append')
    args = parser.parse_args()
    audit(args.root, args.output, args.export_case)
    if args.tar:
        with tarfile.open(args.tar, 'w:gz') as archive:
            for path in sorted(args.output.iterdir()):
                archive.add(path, arcname=path.name)
