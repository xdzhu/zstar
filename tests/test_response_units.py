"""Field boundaries, tensor coupling and explicit Raman normalization."""

import numpy as np
import pytest

from zstar.response_units import raman_convention, slab_effective_dielectric


def supercell_response(layer, fraction):
    # Solve interface continuity for three independent (E_x, E_y, D_z) fields.
    averages_e, averages_d = [], []
    for source in np.eye(3):
        layer_e = source.copy()
        layer_e[2] = (source[2] - layer[2, :2] @ source[:2]) / layer[2, 2]
        layer_d = layer @ layer_e
        averages_e.append(fraction * layer_e + (1-fraction) * source)
        averages_d.append(fraction * layer_d + (1-fraction) * source)
    return np.array(averages_d).T @ np.linalg.inv(np.array(averages_e).T)


def test_diagonal_slab_parallel_and_inverse_perpendicular():
    sc = np.diag([2.25, 2., 1.2])
    np.testing.assert_allclose(slab_effective_dielectric(sc, 20., 5.), np.diag([6., 5., 3.]))


def test_coupled_slab_tensor_and_rotated_normal():
    layer = np.array([[6., .5, .3], [.5, 5., .2], [.3, .2, 3.]])
    sc = supercell_response(layer, .25)
    np.testing.assert_allclose(slab_effective_dielectric(sc, 20., 5.), layer, atol=1e-14)
    angle = .63
    r = np.array([[np.cos(angle), 0, np.sin(angle)], [0, 1, 0],
                  [-np.sin(angle), 0, np.cos(angle)]])
    np.testing.assert_allclose(
        slab_effective_dielectric(r @ sc @ r.T, 20., 5., normal=r[:, 2]),
        r @ layer @ r.T, atol=1e-14,
    )


def test_frequency_dependent_complex_slab():
    layers = np.array([np.diag([4+1j, 5+2j, 2+.1j]), np.diag([3+.2j, 4+.3j, 2+.2j])])
    sc = np.array([np.diag([1+.2*(e[0,0]-1), 1+.2*(e[1,1]-1),
                           1/(.8+.2/e[2,2])]) for e in layers])
    np.testing.assert_allclose(slab_effective_dielectric(sc, 25., 5.), layers)


@pytest.mark.parametrize('thickness', [0, -1, 21, np.nan, np.inf])
def test_invalid_slab_thickness(thickness):
    with pytest.raises(ValueError, match='thickness'):
        slab_effective_dielectric(np.eye(3), 20., thickness)


@pytest.mark.parametrize('dim,unit,factor', [
    (0, 'angstrom^3', '4*pi'), (1, 'angstrom^2', '4*pi'),
    (2, 'angstrom', 'L_perp'), (3, '1', 'epsilon_r'),
])
def test_raman_units_are_explicit(dim, unit, factor):
    result = raman_convention(dim)
    assert result['response_unit'] == unit
    assert factor in result['response_definition']
    assert result['normal_coordinate_unit'] == 'angstrom*sqrt(amu)'
    assert not result['local_field_correction_added']
def test_slab_spectrum_preserves_default_and_converts_total_tensor():
    from zstar.spectra import GammaModes, BornData, calculate_ir_spectrum
    from zstar.response_units import slab_effective_dielectric
    modes = GammaModes(frequencies_thz=np.array([8., 9., 10.]),
                       eigenvectors=np.eye(3)[:, None, :], masses_amu=np.ones(1),
                       lattice_angstrom=np.diag([10., 10., 20.]), symbols=('X',),
                       positions_fractional=np.zeros((1, 3)))
    born = BornData(tensors=np.eye(3)[None, :, :]*.1,
                    electronic_dielectric=np.diag([1.4, 1.5, 1.1]), source='test')
    sheet = calculate_ir_spectrum(modes, born, dimensionality=2, points=11)
    with pytest.warns(UserWarning):
        normalized = calculate_ir_spectrum(modes, born, dimensionality=2, thickness_angstrom=5., points=11)
    np.testing.assert_allclose(normalized.response_real, np.eye(3)+sheet.response_real/5.)
    converted = calculate_ir_spectrum(modes, born, dimensionality=2, thickness_angstrom=5.,
                                     slab_boundary='macroscopic', points=11)
    expected = slab_effective_dielectric(np.eye(3)+(sheet.response_real+1j*sheet.response_imag)/20., 20., 5.)
    np.testing.assert_allclose(converted.response_real+1j*converted.response_imag, expected)
