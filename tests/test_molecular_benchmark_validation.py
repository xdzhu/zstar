import numpy as np
import pytest
from scipy.constants import elementary_charge, epsilon_0

from tools.shared_response.molecular_validation import internal_basis, molecular_response


@pytest.mark.parametrize('positions,masses,rank', [
    ([[0, 0, 0]], [4], 3),
    ([[-1, 0, 0], [1, 0, 0]], [1, 2], 5),
    ([[0, 0, 0], [0, 1, 1], [0, -1, 1]], [16, 1, 1], 6),
])
def test_rigid_subspace_rank(positions, masses, rank):
    internal, rigid = internal_basis(positions, masses)
    assert rigid.shape[1] == rank
    assert internal.shape[1] == 3*len(masses)-rank
    np.testing.assert_allclose(internal.T @ rigid, 0, atol=1e-14)
    np.testing.assert_allclose(internal @ internal.T + rigid @ rigid.T,
                               np.eye(3*len(masses)), atol=1e-14)


def diatomic_data():
    h = np.zeros((6, 6))
    h[np.ix_([0, 3], [0, 3])] = [[2, -2], [-2, 2]]
    z = np.zeros((2, 3, 3))
    z[:, 0, 0] = [1, -1]
    return h, z, np.array([[-1., 0, 0], [1, 0, 0]])


def test_diatomic_static_units_and_mass_invariance():
    h, z, positions = diatomic_data()
    a = molecular_response(h, z, positions, [1, 1])
    b = molecular_response(h, z, positions, [2, 2])
    expected = elementary_charge / (4*np.pi*epsilon_0) * 1e10 / 2
    assert a['tensor'][0][0] == pytest.approx(expected)
    np.testing.assert_allclose(a['tensor'], b['tensor'], atol=1e-12)
    assert b['frequencies_cm1'][0] == pytest.approx(a['frequencies_cm1'][0]/np.sqrt(2))
    assert a['mode_sum_closure_relative'] < 1e-14


def test_covariance_and_translation_invariance():
    h, z, positions = diatomic_data()
    rotation = np.array([[0., -1, 0], [1, 0, 0], [0, 0, 1]])
    r = np.kron(np.eye(2), rotation)
    a = molecular_response(h, z, positions, [1, 2])
    b = molecular_response(r @ h @ r.T, rotation @ z @ rotation.T,
                           positions @ rotation.T + 14, [1, 2])
    np.testing.assert_allclose(b['tensor'], rotation @ np.array(a['tensor']) @ rotation.T, atol=1e-12)
    np.testing.assert_allclose(a['frequencies_cm1'], b['frequencies_cm1'])


def test_no_projection_can_hide_internal_instability_or_unrelaxed_geometry():
    h, z, positions = diatomic_data()
    assert molecular_response(-h, z, positions, [1, 1])['status'] == 'singular_or_unstable_internal_subspace'
    result = molecular_response(h, z, positions, [1, 1], [[.1, 0, 0], [-.1, 0, 0]])
    assert result['status'] == 'nonstationary_reference'
    assert 'tensor' not in result


def test_invalid_masses_rejected():
    with pytest.raises(ValueError, match='positive'):
        internal_basis([[0, 0, 0]], [0])
