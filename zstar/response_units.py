"""Explicit response conventions and a screened slab capacitor conversion."""

from __future__ import annotations

import numpy as np


def raman_convention(dimension: int) -> dict:
    """Describe the existing derivatives without silently changing their scale."""
    definitions = {
        0: ('V*(epsilon_r-I)/(4*pi)', 'angstrom^3', 'Gaussian polarizability volume'),
        1: ('A_perp*(epsilon_r-I)/(4*pi)', 'angstrom^2', 'Gaussian line polarizability'),
        2: ('L_perp*(epsilon_r-I)', 'angstrom', 'SI alpha_2D/epsilon_0, source-field normalized'),
        3: ('epsilon_r', '1', 'relative dielectric tensor'),
    }
    if dimension not in definitions:
        raise ValueError('dimensionality must be 0, 1, 2, or 3')
    definition, unit, convention = definitions[dimension]
    return {'response_definition': definition, 'response_unit': unit,
            'response_convention': convention,
            'normal_coordinate': 'q', 'normal_coordinate_unit': 'angstrom*sqrt(amu)',
            'tensor_unit': f'{unit}/(angstrom*sqrt(amu))',
            'field_convention': 'inherited from the electronic-response source',
            'local_field_correction_added': False}


def slab_effective_dielectric(epsilon_supercell, cell_height: float,
                              thickness: float, *, normal=(0., 0., 1.)):
    """Remove vacuum from a *screened macroscopic* slab dielectric tensor.

    The input must describe response to the supercell-averaged electric field.
    This is not a local-field correction to an independent-particle spectrum.
    Tangential E and normal D are continuous at the vacuum/slab interface.
    Off-diagonal couplings are retained rather than inverted elementwise.
    """
    if not np.isfinite(cell_height) or not np.isfinite(thickness) or not 0 < thickness <= cell_height:
        raise ValueError('Slab thickness must be positive and no larger than the cell height')
    epsilon = np.asarray(epsilon_supercell)
    if epsilon.shape[-2:] != (3, 3) or not np.all(np.isfinite(epsilon)):
        raise ValueError('Expected finite (..., 3, 3) dielectric tensors')
    n = np.asarray(normal, dtype=float)
    if n.shape != (3,) or not np.all(np.isfinite(n)) or np.linalg.norm(n) < 1e-12:
        raise ValueError('Slab normal must be a finite nonzero three-vector')
    n = n / np.linalg.norm(n)
    tangent = np.eye(3)[np.argmin(np.abs(n))]
    tangent -= np.dot(tangent, n) * n
    tangent /= np.linalg.norm(tangent)
    basis = np.column_stack((tangent, np.cross(n, tangent), n))
    local = np.einsum('ia,...ij,jb->...ab', basis, epsilon, basis)
    nn = local[..., 2, 2]
    if np.any(np.abs(nn) < 1e-12):
        raise ValueError('Singular perpendicular supercell permittivity')
    mixed = np.empty_like(local, dtype=np.result_type(local, float))
    mixed[..., :2, :2] = local[..., :2, :2] - (
        local[..., :2, 2, None] * local[..., None, 2, :2] / nn[..., None, None])
    mixed[..., :2, 2] = local[..., :2, 2] / nn[..., None]
    mixed[..., 2, :2] = -local[..., 2, :2] / nn[..., None]
    mixed[..., 2, 2] = 1 / nn
    layer = np.eye(3) + (cell_height / thickness) * (mixed - np.eye(3))
    ll = layer[..., 2, 2]
    if np.any(np.abs(ll) < 1e-12):
        raise ValueError('Singular effective perpendicular permittivity at the supplied thickness')
    result = np.empty_like(layer)
    result[..., 2, 2] = 1 / ll
    result[..., :2, 2] = layer[..., :2, 2] / ll[..., None]
    result[..., 2, :2] = -layer[..., 2, :2] / ll[..., None]
    result[..., :2, :2] = layer[..., :2, :2] - (
        layer[..., :2, 2, None] * layer[..., None, 2, :2] / ll[..., None, None])
    return np.einsum('ia,...ab,jb->...ij', basis, result, basis)
