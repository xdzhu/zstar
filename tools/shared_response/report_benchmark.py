"""Summarize measured shared/control runs; never substitute absent results."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from analyze_archive import static_response
from zstar.shared_response import make_phonopy, read_structure, reconstruct_responses, symmetry_operations
from zstar.pyatb_compat import read_static_dielectric

ROOT = Path('/home/zhuxd/abacus/agent-runs/20260904-shared-response-benchmark')


def result(root):
    metadata = root / 'shared_response.json'
    if not metadata.is_file():
        return {'status': 'not_prepared'}
    meta = json.loads(metadata.read_text())
    states = [json.loads(p.read_text()) for p in sorted((root / '.zstar/stages').glob('*.json'))]
    summary = {'status': 'running_or_pending', 'displacements': len(meta['stages']),
               'total_stages': len(meta['stages'])+1,
               'completed_stages': sum(s['status'] == 'completed' for s in states),
               'failed_stages': [s for s in states if s['status'] == 'failed']}
    summary['full_precision_polarization_stages'] = sum(
        (root/name/'pyatb/Out/Polarization/zstar_precision.json').is_file()
        for name in ['0.no-move'] + [s['name'] for s in meta['stages']])
    summary['output_precision_diagnostic_core_hours'] = sum(
        json.loads(path.read_text())['reserved_core_hours']
        for path in root.glob('*/pyatb-precision/precision_timing.json'))
    rounded_checks = []
    for name in ['0.no-move']+[s['name'] for s in meta['stages']]:
        base = root/name/'pyatb/Out'
        before, repeated = base/'Polarization.before_precision/polarization.dat', base/'Polarization/polarization.rounded.dat'
        if before.is_file() and repeated.is_file():
            rounded_checks.append(before.read_bytes() == repeated.read_bytes())
    summary['precision_rerun_legacy_text_check'] = {'checked_stages': len(rounded_checks), 'identical_stages': sum(rounded_checks)}
    if (root/'refinement_cost.json').is_file():
        summary['berry_refinement_cost'] = json.loads((root/'refinement_cost.json').read_text())
        summary['berry_refinement'] = json.loads((root/'refinement.json').read_text())
    timing_path = root / 'component_times.jsonl'
    times = [json.loads(line) for line in timing_path.read_text().splitlines()] if timing_path.is_file() else []
    groups = {}
    for item in times:
        command = item['command']
        kind = ('SCF' if command.endswith('/abacus') else
                'input_preparation' if 'pyatb_input' in command else
                'band_gate' if str(item['cwd']).endswith('pyatb-band') else 'polarization_and_electronic_response')
        group = groups.setdefault(kind, {'wall_seconds': 0., 'reserved_core_hours': 0., 'calls': 0})
        group['wall_seconds'] += item['wall_seconds']
        # The 40-core worker remains reserved even during its input preparation.
        group['reserved_core_hours'] += item['wall_seconds'] * 40 / 3600
        group['calls'] += 1
    summary['timing_completed_components'] = groups
    summary['reserved_core_hours_completed_components'] = sum(g['reserved_core_hours'] for g in groups.values())
    path = root / 'shared_response_result.json'
    if not path.is_file():
        return summary
    data = json.loads(path.read_text())
    atoms = read_structure(root / '0.no-move/STRU')
    p = make_phonopy(atoms, symprec=meta['symprec_A'])
    raw = reconstruct_responses(len(atoms), data['observations'], symmetry_operations(p, dimension=meta['dimension']),
                                reference_forces=data['reference_forces_eV_A'])
    hessian = raw.force_constants.transpose(0,2,1,3).reshape(len(atoms)*3, len(atoms)*3)
    geometry = {'cell_A': atoms.cell, 'masses_amu': atoms.masses}
    if meta['dimension'] == 0:
        from molecular_validation import molecular_response
        static = molecular_response(hessian, raw.born, atoms.positions, atoms.masses,
                                    data['reference_forces_eV_A'])
    else:
        static = static_response(hessian, raw.born, geometry, meta['dimension'])
    p.dataset = {'natom': len(atoms), 'first_atoms': [
        {'number': s['atom'], 'displacement': np.array(s['displacement_A']),
         'forces': np.array(s['forces_eV_A'])-data['reference_forces_eV_A']}
        for s in data['observations']]}
    p.produce_force_constants(fc_calculator='traditional')
    independent_fc_error = float(np.max(np.abs(p.force_constants.transpose(1,0,3,2)-raw.force_constants)))
    epsilon, _ = read_static_dielectric(root / '0.no-move/pyatb')
    summary.update({'status': 'joint_result_available', 'dimension': meta['dimension'],
                    'diagnostics': data['diagnostics'], 'frequencies_cm1': (np.array(data['frequencies_THz']) * 33.3564095198152).tolist(),
                    'born_raw_e': data['born_raw_e'], 'born_projected_e': data['born_projected_e'],
                    'phonon_static_response': static, 'epsilon_infinity': epsilon.tolist(),
                    'hessian_raw_eV_A2': hessian.tolist(), 'cell_A': atoms.cell.tolist(),
                    'masses_amu': atoms.masses.tolist()})
    summary['independent_phonopy_raw_force_constants_max_difference_eV_A2'] = independent_fc_error
    if meta['dimension'] == 0:
        summary['molecular_internal_validation'] = static
    if static['status'] == 'computed' and meta['dimension'] != 0:
        from zstar.spectra import load_gamma_modes, read_born_data, BornData, calculate_ir_spectrum, mode_effective_charges
        born = read_born_data(root/'BORN', natoms=len(atoms))
        phonon_only = BornData(tensors=born.tensors, electronic_dielectric=None, source=born.source)
        modes = load_gamma_modes(root/'qpoints.yaml')
        spectrum = calculate_ir_spectrum(modes, phonon_only,
            dimensionality=meta['dimension'], acoustic_cutoff_cm1=.001,
            imaginary_tolerance_cm1=.001, points=3)
        tensor = spectrum.response_real[0]
        if meta['dimension'] == 3:
            tensor = tensor - np.eye(3)
        else:
            tensor = tensor/(4*np.pi)
        summary['independent_mode_sum_phonon_static_tensor'] = tensor.tolist()
        summary['mode_sum_vs_Hessian_static_relative_difference'] = float(
            np.linalg.norm(tensor-static['tensor'])/np.linalg.norm(static['tensor']))
        charges = mode_effective_charges(modes, born.tensors)
        groups = []
        for index, frequency in enumerate(modes.frequencies_cm1):
            if frequency <= .001:
                continue
            oscillator = np.outer(charges[index], charges[index].conj()).real
            if groups and abs(frequency-groups[-1]['last_frequency_cm1']) < .01:
                groups[-1]['multiplicity'] += 1
                groups[-1]['oscillator_e2_amu'] += oscillator
                groups[-1]['last_frequency_cm1'] = float(frequency)
            else:
                groups.append({'last_frequency_cm1': float(frequency), 'multiplicity': 1,
                               'oscillator_e2_amu': oscillator})
        for group in groups:
            group['oscillator_e2_amu'] = group['oscillator_e2_amu'].tolist()
        summary['degenerate_mode_oscillator_groups'] = groups
    if static['status'] == 'computed' and meta['dimension'] == 3:
        summary['epsilon_total_static'] = (epsilon + np.array(static['tensor'])).tolist()
    return summary


def compare(a, b):
    if any(r['status'] != 'joint_result_available' for r in (a, b)):
        return {'status': 'waiting_for_both_joint_results'}
    if not np.allclose(a['cell_A'], b['cell_A'], rtol=0, atol=1e-8):
        return {'status': 'geometry_mismatch'}
    dz = np.array(a['born_projected_e']) - b['born_projected_e']
    dh = np.array(a['hessian_raw_eV_A2']) - b['hessian_raw_eV_A2']
    df = np.array(a['frequencies_cm1']) - b['frequencies_cm1']
    output = {'status': 'compared', 'max_BEC_difference_raw_e': float(np.max(np.abs(np.array(a['born_raw_e'])-b['born_raw_e']))),
              'max_BEC_difference_projected_e': float(np.max(np.abs(dz))),
              'max_Hessian_difference_eV_A2': float(np.max(np.abs(dh))),
              'relative_Hessian_Frobenius_difference': float(np.linalg.norm(dh)/np.linalg.norm(b['hessian_raw_eV_A2'])),
              'max_frequency_difference_cm1': float(np.max(np.abs(df))),
              'measured_shared_over_cartesian_core_hours': a['reserved_core_hours_completed_components']/b['reserved_core_hours_completed_components']}
    output['both_have_full_precision_polarization'] = all(
        r['full_precision_polarization_stages'] == r['total_stages'] for r in (a, b))
    ma, mb = a.get('molecular_internal_validation'), b.get('molecular_internal_validation')
    if ma and mb:
        output['molecular_internal_mode_count'] = ma['internal_mode_count']
        output['max_internal_frequency_difference_cm1'] = float(np.max(np.abs(
            np.array(ma['frequencies_cm1']) - mb['frequencies_cm1']), initial=0))
    ga, gb = a.get('degenerate_mode_oscillator_groups', []), b.get('degenerate_mode_oscillator_groups', [])
    if ga and [g['multiplicity'] for g in ga] == [g['multiplicity'] for g in gb]:
        oa, ob = (np.array([g['oscillator_e2_amu'] for g in groups]) for groups in (ga, gb))
        output['max_group_oscillator_difference_e2_amu'] = float(np.max(np.abs(oa-ob)))
        output['relative_group_oscillator_Frobenius_difference'] = float(np.linalg.norm(oa-ob)/max(np.linalg.norm(ob), 1e-30))
        output['mode_group_pairing'] = 'ascending frequency with equal degeneracies; within-group oscillator tensors are summed'
    sa, sb = a['phonon_static_response'], b['phonon_static_response']
    if sa['status'] == sb['status'] == 'computed':
        delta = np.array(sa['tensor']) - sb['tensor']
        output['max_phonon_static_difference'] = float(np.max(np.abs(delta)))
        output['relative_phonon_static_Frobenius_difference'] = float(np.linalg.norm(delta)/np.linalg.norm(sb['tensor']))
        output['static_unit'] = sa['unit']
    else:
        output['static_comparison'] = {'status': 'not_computed',
                                       'unified_reason': sa['status'],
                                       'cartesian_reason': sb['status']}
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT/'benchmark_summary.json')
    args = parser.parse_args()
    report = {'timing_basis': 'monotonic completed command durations x 40 reserved CPU cores; relaxation and contended attempts reported separately', 'cases': {}}
    for case in ('sic', 'hfo2', 'in2se3'):
        a, b = result(ROOT/case/'shared'), result(ROOT/case/'cartesian')
        report['cases'][case] = {'shared': a, 'cartesian': b, 'comparison': compare(a, b)}
        if (ROOT/case/'shared-mesh88/shared_response_result.json').is_file() and (ROOT/case/'cartesian-mesh88/shared_response_result.json').is_file():
            dense_a, dense_b = result(ROOT/case/'shared-mesh88'), result(ROOT/case/'cartesian-mesh88')
            comparison = compare(dense_a, dense_b)
            comparison.pop('measured_shared_over_cartesian_core_hours', None)
            report['cases'][case]['dense_mesh'] = {'shared': dense_a, 'cartesian': dense_b, 'comparison': comparison}
            print('  dense-mesh comparison:', comparison)
        if (ROOT/case/'shared-half/shared_response_result.json').is_file():
            half = result(ROOT/case/'shared-half')
            report['cases'][case]['half_step'] = half
            report['cases'][case]['step_comparison'] = compare(a, half)
            # This ratio is a step-convergence cost, not a Cartesian reduction.
            report['cases'][case]['step_comparison'].pop('measured_shared_over_cartesian_core_hours', None)
            print('  half-step comparison:', report['cases'][case]['step_comparison'])
            if (ROOT/case/'cartesian-half/shared_response_result.json').is_file():
                control_half = result(ROOT/case/'cartesian-half')
                report['cases'][case]['cartesian_half'] = control_half
                report['cases'][case]['matched_half_step_comparison'] = compare(half, control_half)
                report['cases'][case]['cartesian_step_comparison'] = compare(b, control_half)
                print('  matched half-step:', report['cases'][case]['matched_half_step_comparison'])
        print(case, a['status'], a.get('completed_stages'), '/', a.get('total_stages'),
              b['status'], b.get('completed_stages'), '/', b.get('total_stages'))
        if a['status'] == 'joint_result_available':
            print('  shared BEC diagonal:', np.diagonal(a['born_projected_e'], axis1=1, axis2=2))
            print('  frequencies:', np.round(a['frequencies_cm1'], 3))
            print('  static:', a['phonon_static_response'])
        print('  comparison:', report['cases'][case]['comparison'])
        print('  reserved core-hours:', a.get('reserved_core_hours_completed_components'),
              b.get('reserved_core_hours_completed_components'))
    report['in2se3_Berry_mesh_diagnostics'] = []
    for nk in (44,66,88):
        path = ROOT/'in2se3'/f'berry-mesh-{nk}'/'atom-3/comparison.json'
        if path.is_file():
            report['in2se3_Berry_mesh_diagnostics'].append(json.loads(path.read_text()))
    args.output.write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
