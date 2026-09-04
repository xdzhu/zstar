"""Independent, fixed-orientation molecular response audit for the benchmarks.

Raw production Hessians are never overwritten. Rigid-motion projection follows
the mass-weighted vibrational analysis described at gaussian.com/vib/.
"""

from __future__ import annotations

import numpy as np
from scipy.constants import atomic_mass, elementary_charge, epsilon_0, speed_of_light


def internal_basis(positions, masses):
    positions = np.asarray(positions, dtype=float)
    masses = np.asarray(masses, dtype=float)
    if (positions.shape != (len(masses), 3) or not len(masses)
            or not np.isfinite(positions).all() or not np.isfinite(masses).all()
            or np.any(masses <= 0)):
        raise ValueError('Finite Cartesian coordinates and positive atomic masses are required')
    centered = positions - np.average(positions, axis=0, weights=masses)
    weights = np.sqrt(masses)[:, None]
    translations = np.tile(np.eye(3), (len(masses), 1)) * weights.repeat(3, axis=0)
    rotations = np.column_stack([
        (np.cross(axis, centered) * weights).ravel() for axis in np.eye(3)
    ])
    rigid = np.column_stack([translations, rotations])
    norms = np.linalg.norm(rigid, axis=0)
    rigid = rigid[:, norms > 1e-12] / norms[norms > 1e-12]
    vectors, singular, _ = np.linalg.svd(rigid, full_matrices=True)
    rank = int(np.count_nonzero(singular > singular[0] * 1e-10))
    return vectors[:, rank:], vectors[:, :rank]


def molecular_response(hessian, born, positions, masses, reference_forces=None):
    """Return internal modes and Gaussian vibrational polarizability in A^3.

    BEC/APT arrays use [atom, dipole direction, displacement direction].
    Coordinates must describe one contiguous neutral molecule, not wrapped atoms.
    This is a fixed-orientation response; free-rotor orientational polarization
    is not included. Large residual forces invalidate the equilibrium audit.
    """
    masses = np.asarray(masses, dtype=float)
    internal, rigid = internal_basis(positions, masses)
    n = len(masses)
    hessian = np.asarray(hessian, dtype=float)
    born = np.asarray(born, dtype=float)
    if (hessian.shape != (3*n, 3*n) or born.shape != (n, 3, 3)
            or not np.isfinite(hessian).all() or not np.isfinite(born).all()):
        raise ValueError('Incompatible or non-finite Hessian/APT arrays')
    mass = np.repeat(np.sqrt(masses), 3)
    dynamical = (hessian + hessian.T) / 2 / mass[:, None] / mass[None, :]
    reduced = internal.T @ dynamical @ internal
    values, vectors = np.linalg.eigh(reduced)
    conversion = np.sqrt(elementary_charge / (atomic_mass * 1e-20)) / (2*np.pi*speed_of_light*100)
    frequencies = np.sign(values) * np.sqrt(np.abs(values)) * conversion
    data = {
        'status': 'computed', 'unit': 'Angstrom^3 Gaussian vibrational polarizability',
        'rigid_motion_rank': rigid.shape[1], 'internal_mode_count': len(values),
        'frequencies_cm1': frequencies.tolist(),
        'raw_rigid_coupling_relative': float(np.linalg.norm(dynamical @ rigid)
                                            / max(np.linalg.norm(dynamical), 1e-30)),
        'coordinate_convention': 'fixed orientation, mass-weighted COM Eckart subspace',
        'raw_Hessian_modified': False,
    }
    if reference_forces is not None:
        forces = np.asarray(reference_forces, dtype=float)
        if forces.shape != (n, 3) or not np.isfinite(forces).all():
            raise ValueError('Invalid reference force array')
        data['reference_max_force_eV_A'] = float(np.linalg.norm(forces, axis=1).max())
        if data['reference_max_force_eV_A'] > .005:
            data['status'] = 'nonstationary_reference'
            return data
    if len(values) and values.min() <= max(1e-12, np.abs(values).max() * 1e-8):
        data['status'] = 'singular_or_unstable_internal_subspace'
        return data
    neutral = born - born.mean(axis=0)
    coupling = (np.concatenate(neutral, axis=1) / mass) @ internal
    prefactor = elementary_charge / (4*np.pi*epsilon_0) * 1e10
    tensor = prefactor * coupling @ np.linalg.solve(reduced, coupling.T)
    charges = coupling @ vectors
    mode_sum = prefactor * (charges / values) @ charges.T
    data.update(tensor=tensor.tolist(), mode_charges_e_sqrt_amu=charges.T.tolist(),
                mode_sum_closure_relative=float(np.linalg.norm(tensor-mode_sum)
                                               / max(np.linalg.norm(tensor), 1e-30)))
    return data


if __name__ == '__main__':
    import argparse
    import json
    from pathlib import Path
    from zstar.shared_response import read_structure, make_phonopy, reconstruct_responses, symmetry_operations

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    meta = json.loads((args.root / 'shared_response.json').read_text())
    if meta['dimension'] != 0:
        raise ValueError('This audit requires a molecular ensemble (dim=0)')
    data = json.loads((args.root / 'shared_response_result.json').read_text())
    atoms = read_structure(args.root / '0.no-move/STRU')
    phonon = make_phonopy(atoms, symprec=meta['symprec_A'])
    raw = reconstruct_responses(len(atoms), data['observations'],
        symmetry_operations(phonon, dimension=0), reference_forces=data['reference_forces_eV_A'])
    h = raw.force_constants.transpose(0, 2, 1, 3).reshape(3*len(atoms), 3*len(atoms))
    report = molecular_response(h, raw.born, atoms.positions, atoms.masses, data['reference_forces_eV_A'])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))
    if report['status'] != 'computed':
        raise SystemExit('Molecular equilibrium/internal-mode validation failed')
