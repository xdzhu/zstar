"""Analytic safeguards for the research-only retrospective analysis."""

import numpy as np
import pytest
import spglib

from analyze_archive import (
    expand_charges, fit_charge, fit_force_constants, frequencies, operations,
    orbit_matrix, project_hessian, select_stages, static_response,
)
from collect_archive import read_electronic_dielectric


def crystal_groups():
    groups = {}
    for hall in range(1, 531):
        symbol = spglib.get_spacegroup_type(hall).pointgroup_international
        if symbol not in groups:
            rotations = np.unique(spglib.get_symmetry_from_database(hall)["rotations"], axis=0)
            metric = sum(w.T @ w for w in rotations) / len(rotations)
            s = np.linalg.cholesky(metric).T
            groups[symbol] = [(s @ w @ np.linalg.inv(s), np.array([0])) for w in rotations]
    return groups


GROUPS = crystal_groups()


@pytest.mark.parametrize("symbol", list(GROUPS))
def test_all_32_groups_recover_general_allowed_tensor(symbol):
    ops = GROUPS[symbol]
    rng = np.random.default_rng(917)
    raw = rng.normal(size=(3, 3))
    z = sum(r @ raw @ r.T for r, _ in ops) / len(ops)
    stages = []
    for axis in range(3):
        for sign in (1, -1):
            u = np.eye(3)[axis] * 0.01 * sign
            stages.append({"name": "xyz"[axis] + ("+" if sign > 0 else "-"), "displacement_A": u, "dipole_change_e_A": z @ u})
    chosen, _ = select_stages(stages, ops)
    assert np.linalg.matrix_rank(orbit_matrix(chosen, ops), tol=1e-7) == 3
    np.testing.assert_allclose(fit_charge(chosen, ops), z, atol=1e-12)


def test_rank_one_svd_is_not_rank_three():
    with pytest.raises(ValueError, match="span three"):
        select_stages([{"name": "x+", "displacement_A": [0.01, 0, 0]}], GROUPS["1"])


def test_no_symmetry_needs_all_axes_and_real_minus_for_central():
    data = [{"name": a + s, "displacement_A": np.eye(3)[i] * (0.01 if s == "+" else -0.01)} for i, a in enumerate("xyz") for s in ("+", "-")]
    assert len(select_stages(data, GROUPS["1"])[0]) == 6


def test_charge_tensor_must_not_be_forced_symmetric():
    z = np.array([[1, 2, 3], [-4, 5, 6], [7, -8, 9.]])
    data = [{"name": "xyz"[i] + "+", "displacement_A": u, "dipole_change_e_A": z @ u} for i, u in enumerate(np.eye(3) * 0.01)]
    np.testing.assert_allclose(fit_charge(data, GROUPS["1"]), z)


def test_molecular_force_permutation_and_common_ensemble():
    positions = np.array([[0, 0, 0], [1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1.]])
    geom = {"cell_A": (np.eye(3) * 20).tolist(), "fractional_positions": (positions / 20 + .5).tolist(), "symbols": ["C"] + ["H"] * 4, "masses_amu": [12.] + [1.] * 4}
    ops, _ = operations(geom, 0)
    assert len(ops) == 24
    rng = np.random.default_rng(82)
    h = rng.normal(size=(15, 15))
    h = h.T @ h
    phi = np.zeros_like(h)
    for r, perm in ops:
        u = np.zeros_like(h)
        for i, j in enumerate(perm):
            u[3*j:3*j+3, 3*i:3*i+3] = r
        phi += u @ h @ u.T / len(ops)
    phi = project_hessian(phi)
    data = []
    for atom in (0, 1):
        for axis in range(3):
            for sign in (1, -1):
                d = np.eye(3)[axis] * 0.01 * sign
                u = np.zeros((5, 3)); u[atom] = d
                data.append({"name": f"{atom}/" + "xyz"[axis] + ("+" if sign == 1 else "-"), "atom": atom, "displacement_A": d, "forces_eV_A": (-phi @ u.ravel()).reshape(5, 3)})
    selected = []
    for atom in (0, 1):
        selected += select_stages([s for s in data if s["atom"] == atom], [(r, p) for r, p in ops if p[atom] == atom])[0]
    fit, _ = fit_force_constants(selected, ops, 5)
    np.testing.assert_allclose(fit, phi, atol=1e-11)
    np.testing.assert_allclose(frequencies(fit, geom["masses_amu"])[3:], frequencies(phi, geom["masses_amu"])[3:], atol=1e-8)


def test_neutrality_uses_all_atoms_not_representative_count():
    ops = [(np.eye(3), np.array([0, 1, 2])), (np.eye(3), np.array([0, 2, 1]))]
    full = expand_charges({0: np.eye(3) * 2, 1: -np.eye(3)}, ops, 3)
    np.testing.assert_allclose(full.sum(axis=0), 0)


def test_static_pseudoinverse_matches_mass_weighted_mode_sum():
    rng = np.random.default_rng(71)
    x = rng.normal(size=(9, 9))
    phi = project_hessian(x.T @ x)
    z = rng.normal(size=(3, 3, 3)); z -= z.mean(axis=0)
    masses = np.array([1., 12., 16.])
    m = np.repeat(np.sqrt(masses), 3)
    w, v = np.linalg.eigh(phi / m[:, None] / m[None, :])
    b = np.concatenate(z, axis=1)
    active = w > 1e-8
    mode_charge = (b / m) @ v[:, active]
    mode_sum = (mode_charge / w[active]) @ mode_charge.T
    np.testing.assert_allclose(b @ np.linalg.pinv(phi, rcond=1e-7) @ b.T, mode_sum, atol=1e-9)


def test_imaginary_optical_mode_blocks_static_response():
    phi = -project_hessian(np.eye(6))
    geom = {"cell_A": np.eye(3) * 10, "masses_amu": [1, 1]}
    result = static_response(phi, np.array([np.eye(3), -np.eye(3)]), geom, 3)
    assert result["status"] == "unstable_harmonic_reference"


@pytest.mark.parametrize("header", ["# xx xy xz yx yy yz zx zy zz\n", "14.3996\n"])
def test_dielectric_is_not_first_atom_born_tensor(tmp_path, header):
    path = tmp_path / "BORN"
    path.write_text(header + "6 0 0 0 6 0 0 0 6\n2.7 0 0 0 2.7 0 0 0 2.7\n")
    np.testing.assert_allclose(read_electronic_dielectric(path), np.eye(3) * 6)


def test_static_units_against_diatomic_closed_form():
    from scipy.constants import elementary_charge, epsilon_0
    k, charge, volume = 20., 2.7, 25.
    phi = np.kron(np.array([[1., -1.], [-1., 1.]]), np.eye(3) * k)
    z = np.array([np.eye(3), -np.eye(3)]) * charge
    geom = {"cell_A": np.eye(3) * volume**(1/3), "masses_amu": [12., 28.]}
    response = static_response(phi, z, geom, 3)
    expected = elementary_charge / epsilon_0 * 1e10 / volume * charge**2 / k
    np.testing.assert_allclose(response["tensor"], np.eye(3) * expected, atol=1e-12)
