"""Reprocess each complete shared archive in a temporary directory, offline."""

import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import sys

import numpy as np

from zstar.shared_abacus import collect_shared_abacus, prepare_shared_abacus
from zstar.spectra import BornData, calculate_ir_spectrum, load_gamma_modes, read_born_data

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools/shared_response'))


def verify(case, scheme='shared'):
    archive = case/'results'/scheme
    if not (archive/'shared_response_result.json').is_file():
        print(case.name, 'PENDING: no complete archived shared result')
        return False
    checksums = case/'checksums.json'
    if checksums.is_file():
        for name, expected in json.loads(checksums.read_text()).items():
            digest = hashlib.sha256()
            with (case/name).open('rb') as handle:
                for block in iter(lambda: handle.read(1024*1024), b''):
                    digest.update(block)
            if digest.hexdigest() != expected:
                raise ValueError(f'Archive checksum mismatch: {case.name}/{name}')
    expected = json.loads((archive/'shared_response_result.json').read_text())
    with tempfile.TemporaryDirectory(prefix='zstar-shared-check-') as temporary:
        copied = Path(temporary)/'shared'
        shutil.copytree(archive, copied)
        result = collect_shared_abacus(copied)
        np.testing.assert_allclose(result['born_raw_e'], expected['born_raw_e'], atol=1e-8, rtol=0)
        np.testing.assert_allclose(result['born_projected_e'], expected['born_projected_e'], atol=1e-8, rtol=0)
        # Different LAPACK backends can assign tiny signed values to the three
        # translation modes; no broad cutoff is applied to optical modes.
        actual_f = np.asarray(result['frequencies_THz'])
        expected_f = np.asarray(expected['frequencies_THz'])
        # Compare eigenvalues, not their ill-conditioned square roots, at zero.
        zero = np.maximum(np.abs(actual_f), np.abs(expected_f)) < 1e-5
        np.testing.assert_allclose(actual_f[zero]*np.abs(actual_f[zero]),
                                   expected_f[zero]*np.abs(expected_f[zero]), atol=1e-10, rtol=0)
        np.testing.assert_allclose(actual_f[~zero], expected_f[~zero], atol=1e-6, rtol=1e-7)
        if case.name in ('SiC', 't_HfO2', 'alpha_In2Se3'):
            benchmark = json.loads((case.parent/'benchmark_summary.json').read_text())
            key = {'SiC': 'sic', 't_HfO2': 'hfo2', 'alpha_In2Se3': 'in2se3'}[case.name]
            entry = benchmark['cases'][key]
            if scheme == 'shared-mesh88':
                entry = entry['dense_mesh']
            reference = entry['shared']['phonon_static_response']
        elif '-mesh112' in scheme:
            reference = json.loads((case/'results/mesh112-comparison.json').read_text())[scheme.split('-')[0]]['phonon_static_response']
        else:
            entry = json.loads((case.parent/'four_new_benchmarks.json').read_text())['cases'][case.name]
            reference = entry[scheme]['phonon_static_response']
        assert reference['status'] == 'computed'
        if result['dimension'] == 0:
            from report_benchmark import result as audit_result
            checked = audit_result(copied)['phonon_static_response']
            assert checked['status'] == 'computed'
            tensor = checked['tensor']
        else:
            modes = load_gamma_modes(copied/'qpoints.yaml')
            born = read_born_data(copied/'BORN', natoms=len(modes.masses_amu))
            response = calculate_ir_spectrum(modes, BornData(born.tensors, None, born.source),
                dimensionality=result['dimension'], acoustic_cutoff_cm1=.001,
                imaginary_tolerance_cm1=.001, points=3)
            tensor = response.response_real[0]
            tensor = tensor-np.eye(3) if result['dimension'] == 3 else tensor/(4*np.pi)
        np.testing.assert_allclose(tensor, reference['tensor'], atol=2e-6, rtol=1e-6)
        prepare_shared_abacus(case/'run/STRU', root=Path(temporary)/'new-preparation',
                              scf_input=case/'run/INPUT', dimension=result['dimension'])
    print(case.name, scheme, 'PASS: checksums, raw/projected BEC, Gamma frequencies, and static response')
    return True


if __name__ == '__main__':
    root = Path(__file__).resolve().parent
    completed = [verify(root/name) for name in ('SiC', 't_HfO2', 'alpha_In2Se3')]
    if (root/'alpha_In2Se3/results/shared-mesh88/shared_response_result.json').is_file():
        completed.append(verify(root/'alpha_In2Se3', 'shared-mesh88'))
    for name in ('hBN', 'MoS2', 'H2O', 'CH4'):
        for scheme in ('unified', 'cartesian'):
            completed.append(verify(root/name, scheme))
    for scheme in ('unified-mesh112', 'cartesian-mesh112'):
        completed.append(verify(root/'MoS2', scheme))
    if not all(completed):
        raise SystemExit('Some shared example archives are still pending.')
