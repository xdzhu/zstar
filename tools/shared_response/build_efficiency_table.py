"""Build the eight-system table from successful solver timing ledgers."""

import argparse
import hashlib
import json
from pathlib import Path

from report_eight_systems import solver_cost


CASES = [('cubic_BaTiO3', r'cubic BaTiO$_3$'), ('SiC', '3C-SiC'),
         ('t_HfO2', r't-HfO$_2$'), ('alpha_In2Se3', r'$\alpha$-In$_2$Se$_3$'),
         ('hBN', 'hBN'), ('MoS2', r'MoS$_2$'), ('H2O', r'H$_2$O'), ('CH4', r'CH$_4$')]


def entry(root, case):
    base = root / case / 'results'
    if case in ('H2O', 'CH4') and not (base / 'relaxation.json').is_file():
        return {'status': 'pending_relaxed_production_pair'}
    if case == 'cubic_BaTiO3':
        u, c, f = [solver_cost(base / path) for path in ['unified', 'legacy_bec', 'legacy_phonon']]
        costs = {'Cartesian': c['solver_core_hours']+f['solver_core_hours'], 'Unified': u['solver_core_hours']}
        counts = {'Cartesian': [9, 13], 'Unified': [3, 4]}
        files = [base / path / 'component_times.jsonl' for path in ['unified', 'legacy_bec', 'legacy_phonon']]
        protocol = 'forward Cartesian BEC plus independent phonons'
    else:
        old = case in ('SiC', 't_HfO2', 'alpha_In2Se3')
        paths = {'Cartesian': base / ('controls/cartesian' if old else 'cartesian'),
                 'Unified': base / ('shared' if old else 'unified')}
        costs, counts, files = {}, {}, []
        for name, path in paths.items():
            if not (path/'shared_response_result.json').is_file():
                return {'status': 'pending_production_pair'}
            meta = json.loads((path/'shared_response.json').read_text())
            cost = solver_cost(path)
            costs[name] = cost['solver_core_hours']
            counts[name] = [len(meta['stages']), cost['calls']['ABACUS']]
            if counts[name][1] != counts[name][0]+1:
                raise ValueError(f'{case}/{name}: unmatched displacement and SCF counts')
            files.append(path/'component_times.jsonl')
        protocol = 'central Cartesian joint-response control'
    return {'status': 'complete', 'protocol': protocol, 'counts': counts, 'solver_core_hours': costs,
            'speedup': costs['Cartesian']/costs['Unified'],
            'timing_sha256': {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, default=Path('examples/Shared_Response'))
    parser.add_argument('--output', type=Path, default=Path('docs/research/eight_system_efficiency.json'))
    parser.add_argument('--require-complete', action='store_true')
    args = parser.parse_args()
    report = {'definition': 'ABACUS + PYATB successful solver calls, including reference and band gate, excluding input preparation, relaxation, refinements and diagnostics',
              'cases': {case: entry(args.root, case) for case, _ in CASES}}
    if args.require_complete and any(data['status'] != 'complete' for data in report['cases'].values()):
        raise ValueError('Eight final production pairs are required before freezing the table')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    rows = []
    for case, label in CASES:
        data = report['cases'][case]
        if rows:
            rows.append(r'\addlinespace[2pt]')
        for index, scheme in enumerate(['Cartesian', 'Unified']):
            material = rf'\multirow{{2}}{{*}}{{{label}}}' if index == 0 else ''
            if data['status'] == 'complete':
                count, scf = data['counts'][scheme]
                cost = f"{data['solver_core_hours'][scheme]:.3f}"
                speedup = '1.00' if scheme == 'Cartesian' else f"{data['speedup']:.2f}"
            else:
                count = scf = cost = speedup = '--'
            rows.append(f'{material} & {scheme} & {count} & {scf} & {cost} & {speedup} '+r'\\')
    args.output.with_suffix('.tex').write_text('\n'.join(rows)+'\n')
    print('\n'.join(rows))
