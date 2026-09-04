"""Data-only response sensitivity checks for the CPC revision."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from analyze_archive import static_response
from zstar.spectra import load_gamma_modes
from scipy.constants import atomic_mass, elementary_charge, epsilon_0, speed_of_light


def decompose(pair):
    a, b = pair.get('unified', pair.get('shared')), pair['cartesian']
    geometry = {'cell_A': a['cell_A'], 'masses_amu': a['masses_amu']}
    results = {}
    for h_name, h in [('U', a), ('C', b)]:
        for z_name, z in [('U', a), ('C', b)]:
            value = static_response(np.array(h['hessian_raw_eV_A2']), np.array(z['born_raw_e']), geometry, a['dimension'])
            if value['status'] != 'computed':
                return {'status': value['status']}
            results[h_name + z_name] = np.array(value['tensor'])
    base = results['CC']
    return {'status': 'computed', 'index_order': 'Hessian then BEC: U=Unified, C=Cartesian',
            'tensors': {k: v.tolist() for k, v in results.items()},
            'relative_changes_vs_CC': {k: float(np.linalg.norm(v-base)/np.linalg.norm(base)) for k, v in results.items()},
            'signed_diagonal_changes_vs_CC': {k: np.diag(v-base).tolist() for k, v in results.items()},
            'unit': value['unit']}


def hfo2(root):
    base = root / 'examples/IR_Raman_Spectra/Bulk_HfO2/results'
    modes = load_gamma_modes(base / 'qpoints.yaml')
    volume = modes.volume_angstrom3 if hasattr(modes, 'volume_angstrom3') else abs(np.linalg.det(modes.cell))
    prefactor = elementary_charge**2 / (epsilon_0 * volume * 1e-30 * atomic_mass * (2*np.pi*speed_of_light*100)**2)
    rows = list(csv.DictReader((base/'ir/ir_modes.csv').open()))
    groups = []
    for row in rows:
        frequency = float(row['frequency_cm-1'])
        charge = np.array([float(row[f'Zmode_{axis}']) for axis in 'xyz'])
        tensor = prefactor*np.outer(charge, charge)/frequency**2
        if np.linalg.norm(tensor) < 1e-8:
            continue
        if groups and abs(frequency - groups[-1]['frequency_cm-1']) < .001:
            groups[-1]['contribution'] += tensor
            groups[-1]['multiplicity'] += 1
        else:
            groups.append({'frequency_cm-1': frequency, 'multiplicity': 1, 'contribution': tensor})
    total_ph = sum(g['contribution'] for g in groups)
    total = np.array(json.loads((base/'ir/static_response.json').read_text())['tensor'])
    soft = groups[0]['contribution']
    sensitivity = total - soft + soft*(groups[0]['frequency_cm-1']/129)**2
    return {'volume_A3': volume, 'epsilon_ph': total_ph.tolist(), 'epsilon_infinity': (total-total_ph).tolist(),
            'groups': [{**g, 'contribution': g['contribution'].tolist()} for g in groups],
            'soft_pair_fraction_inplane_ph': float(soft[0,0]/total_ph[0,0]),
            'soft_pair_fraction_trace_ph': float(np.trace(soft)/np.trace(total_ph)),
            'original_mean_epsilon': float(np.trace(total)/3),
            'fixed_charge_129cm1_sensitivity_mean_epsilon': float(np.trace(sensitivity)/3),
            'sensitivity_only': 'Replace only the 96.1 cm-1 pair with 129 cm-1 at fixed oscillator strengths; not a corrected prediction.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--new-results', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    data = json.loads((args.root/'examples/Shared_Response/benchmark_summary.json').read_text())
    result = {'In2Se3': {'baseline': decompose(data['cases']['in2se3']),
                        'refined': decompose(data['cases']['in2se3']['dense_mesh'])},
              'HfO2': hfo2(args.root)}
    if args.new_results:
        new = json.loads(args.new_results.read_text())
        result['MoS2'] = decompose(new['cases']['MoS2'])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
