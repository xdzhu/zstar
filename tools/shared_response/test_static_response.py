"""Do not hide an unresolved optical zero mode in a pseudoinverse."""

import numpy as np
from analyze_archive import static_response


def model():
    k = np.diag([1., 2., 3.])
    phi = np.block([[k, -k], [-k, k]])
    z = np.array([np.eye(3), -np.eye(3)])
    geom = {'masses_amu': [1., 2.], 'cell_A': np.diag([5., 5., 5.])}
    return phi, z, geom


def test_extra_optical_zero_is_not_silently_dropped():
    phi, z, geom = model()
    phi[[0, 3], :] = 0
    phi[:, [0, 3]] = 0
    assert static_response(phi, z, geom, 3)['status'] == 'singular_or_unstable_optical_subspace'


def test_sheet_and_wire_normalization():
    phi, z, geom = model()
    bulk = np.array(static_response(phi, z, geom, 3)['tensor'])
    sheet = np.array(static_response(phi, z, geom, 2)['tensor'])
    wire = np.array(static_response(phi, z, geom, 1)['tensor'])
    np.testing.assert_allclose(sheet, bulk * 5 / (4*np.pi))
    np.testing.assert_allclose(wire, bulk * 25 / (4*np.pi))


def test_molecule_requires_rotational_subspace_check():
    assert static_response(*model(), 0)['status'] == 'molecular_rotational_subspace_not_validated'
