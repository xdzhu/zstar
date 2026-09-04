"""Geometry and tensor algebra for shared Gamma displacement responses.

Lengths are Angstrom, forces eV/Angstrom, and dipole changes e*Angstrom.
Born tensors use Z[atom, polarization, displacement], not its transpose.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

from .stru_analyzer import stru_analyzer

BOHR_ANGSTROM = 0.529177210903
DEFAULT_DISTANCE = 0.02 * BOHR_ANGSTROM


def read_structure(path: str | Path) -> PhonopyAtoms:
    values = stru_analyzer(str(path))
    a0, lattice, elements, counts, mode, coords, _, moments, masses, _, _ = values
    if any(np.any(np.asarray(m) != 0) for s in elements for m in moments[s]):
        raise ValueError("Shared displacements require a nonmagnetic reference; magnetic symmetry is not implemented.")
    cell = np.asarray(lattice) * a0 * BOHR_ANGSTROM
    positions = np.asarray([p for s in elements for p in coords[s]])
    if mode == "Cartesian":
        positions = positions * a0 * BOHR_ANGSTROM @ np.linalg.inv(cell)
    return PhonopyAtoms(
        symbols=[s for s in elements for _ in range(counts[s])],
        masses=[float(m) for s, m in zip(elements, masses) for _ in range(counts[s])],
        cell=cell, scaled_positions=positions,
    )


def write_structure(template: str | Path, destination: str | Path, atoms: PhonopyAtoms) -> None:
    """Preserve species/basis declarations and write unambiguous Direct coordinates."""
    source = Path(template)
    values = stru_analyzer(str(source))
    elements, counts, moves = values[2], values[3], values[6]
    expected = [s for s in elements for _ in range(counts[s])]
    if expected != atoms.symbols:
        raise ValueError("Atom ordering differs from the STRU template")
    prefix = source.read_text(encoding="utf-8").split("LATTICE_CONSTANT", 1)[0]
    lines = [prefix.rstrip(), "", "LATTICE_CONSTANT", f"{1 / BOHR_ANGSTROM:.16g}",
             "", "LATTICE_VECTORS"]
    lines.extend(" ".join(f"{x:.16g}" for x in row) for row in atoms.cell)
    lines.extend(["", "ATOMIC_POSITIONS", "Direct"])
    index = 0
    for symbol in elements:
        lines.extend(["", symbol, "0", str(counts[symbol])])
        for movement in moves[symbol]:
            lines.append(" ".join(f"{x:.16g}" for x in atoms.scaled_positions[index])
                         + " m " + " ".join(str(int(x)) for x in movement))
            index += 1
    Path(destination).write_text("\n".join(lines) + "\n", encoding="utf-8")


def actual_displacement(reference: PhonopyAtoms, displaced: PhonopyAtoms,
                        *, tolerance: float = 1e-8) -> tuple[int, np.ndarray]:
    if reference.symbols != displaced.symbols or not np.allclose(
            reference.cell, displaced.cell, atol=tolerance, rtol=0):
        raise ValueError("Shared responses require identical lattices and atom ordering")
    delta = displaced.scaled_positions - reference.scaled_positions
    delta -= np.rint(delta)
    cart = delta @ reference.cell
    changed = np.flatnonzero(np.linalg.norm(cart, axis=1) > tolerance)
    if len(changed) != 1:
        raise ValueError(f"Expected one displaced atom, found {len(changed)}")
    atom = int(changed[0])
    return atom, cart[atom]


def make_phonopy(atoms: PhonopyAtoms, *, symprec: float = 1e-5) -> Phonopy:
    # No ABACUS native-unit interface: geometry and force constants remain in A.
    return Phonopy(atoms, supercell_matrix=np.eye(3, dtype=int),
                   primitive_matrix=np.eye(3), symprec=symprec)


def symmetry_operations(phonon: Phonopy, *, dimension: int = 3):
    """Return Cartesian rotations and source-to-target permutations."""
    atoms = phonon.supercell
    cell, pos = atoms.cell, atoms.scaled_positions
    labels = np.asarray(atoms.numbers)
    ops = []
    tolerance = phonon.symmetry.tolerance
    for w, t in zip(phonon.symmetry.symmetry_operations["rotations"],
                    phonon.symmetry.symmetry_operations["translations"]):
        r = cell.T @ w @ np.linalg.inv(cell.T)
        if not np.allclose(r.T @ r, np.eye(3), atol=1e-7, rtol=0):
            raise ValueError("Approximate lattice symmetry is nonorthogonal; refine the structure or lower --symmprec")
        # Periodic and open directions have different electrostatic boundaries.
        if dimension in (1, 2):
            normal = np.cross(cell[0], cell[1])
            normal /= np.linalg.norm(normal)
            if not np.isclose(abs(normal @ r @ normal), 1, atol=1e-8):
                raise ValueError("Detected symmetry mixes periodic and open directions")
        delta = (pos @ w.T + t)[:, None, :] - pos[None, :, :]
        delta -= np.rint(delta)
        distances = np.linalg.norm(delta @ cell, axis=2)
        distances[labels[:, None] != labels[None, :]] = np.inf
        perm = distances.argmin(axis=1)
        if len(set(perm)) != len(pos) or np.max(distances[np.arange(len(pos)), perm]) > tolerance * 1.01:
            raise ValueError("Ambiguous symmetry atom mapping")
        ops.append((r, perm))
    return ops


@dataclass
class SharedResponse:
    born: np.ndarray
    force_constants: np.ndarray
    diagnostics: dict


def reconstruct_responses(natoms: int, stages: list[dict], operations,
                          *, reference_forces=None) -> SharedResponse:
    """Fit site stabilizer orbits, then transform responses to equivalent atoms.

    All supplied operations must preserve the full reference Hamiltonian and
    electrostatic boundary conditions, not just the atomic coordinates.
    """
    if not stages or not operations:
        raise ValueError("No shared response observations or symmetry operations")
    f0 = np.zeros((natoms, 3)) if reference_forces is None else np.asarray(reference_forces, dtype=float)
    if f0.shape != (natoms, 3) or not np.all(np.isfinite(f0)):
        raise ValueError("Invalid reference forces")
    born = np.zeros((natoms, 3, 3))
    fc = np.zeros((natoms, natoms, 3, 3))
    counts = np.zeros(natoms)
    reports = []
    for atom in sorted({int(s["atom"]) for s in stages}):
        if atom < 0 or atom >= natoms:
            raise ValueError("Displaced atom index is outside the structure")
        stabilizer = [(r, p) for r, p in operations if p[atom] == atom]
        x, y = [], []
        for stage in (s for s in stages if s["atom"] == atom):
            u = np.asarray(stage["displacement_A"], dtype=float)
            dipole = np.asarray(stage["dipole_change_e_A"], dtype=float)
            force = np.asarray(stage["forces_eV_A"], dtype=float) - f0
            if u.shape != (3,) or dipole.shape != (3,) or force.shape != (natoms, 3):
                raise ValueError("Invalid displacement, dipole, or force shape")
            if not all(np.all(np.isfinite(a)) for a in (u, dipole, force)) or np.linalg.norm(u) < 1e-10:
                raise ValueError("Nonfinite response or zero displacement")
            for r, perm in stabilizer:
                rotated = np.empty_like(force)
                rotated[perm] = force @ r.T
                x.append(r @ u)
                y.append(np.r_[r @ dipole, -rotated.ravel()])
        x, y = np.asarray(x), np.asarray(y)
        coef, _, rank, singular = np.linalg.lstsq(x, y, rcond=1e-10)
        if rank != 3:
            raise ValueError(f"Atom {atom + 1}: displacement-orbit rank {rank}/3; more independent directions are required")
        z = coef[:, :3].T
        column = coef[:, 3:].T.reshape(natoms, 3, 3)
        residual = y - x @ coef
        reports.append({"atom": atom, "rank": int(rank),
                        "condition_number": float(singular[0] / singular[-1]),
                        "dipole_fit_max_e_A": float(np.max(np.abs(residual[:, :3]))),
                        "force_fit_max_eV_A": float(np.max(np.abs(residual[:, 3:]))),
                        "observations": len(x)})
        for r, perm in operations:
            target = perm[atom]
            born[target] += r @ z @ r.T
            fc[perm, target] += np.einsum("ab,jbc,dc->jad", r, column, r)
            counts[target] += 1
    if np.any(counts == 0):
        raise ValueError("No response observations for one or more inequivalent atoms")
    born /= counts[:, None, None]
    fc /= counts[None, :, None, None]
    hessian = fc.transpose(0, 2, 1, 3).reshape(3 * natoms, 3 * natoms)
    diagnostics = {"site_fits": reports,
                   "born_asr_max_e": float(np.max(np.abs(born.sum(axis=0)))),
                   "hessian_reciprocity_max_eV_A2": float(np.max(np.abs(hessian - hessian.T))),
                   "hessian_reciprocity_relative_Frobenius": float(np.linalg.norm(hessian-hessian.T)/max(np.linalg.norm(hessian), 1e-30)),
                   "force_asr_max_eV_A2": float(np.max(np.abs(fc.sum(axis=1)))),
                   "total_force_derivative_max_eV_A2": float(np.max(np.abs(fc.sum(axis=0)))),
                   "reference_force_max_eV_A": float(np.max(np.abs(f0))),
                   "constraints_applied": False}
    return SharedResponse(born, fc, diagnostics)


def project_response(response: SharedResponse) -> SharedResponse:
    """Separate, explicit reciprocity/translation projection; retain raw results."""
    n = len(response.born)
    h = response.force_constants.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n)
    t = np.tile(np.eye(3), (n, 1)) / np.sqrt(n)
    q = np.eye(3 * n) - t @ t.T
    h = q @ ((h + h.T) / 2) @ q
    return SharedResponse(
        response.born - response.born.mean(axis=0),
        h.reshape(n, 3, n, 3).transpose(0, 2, 1, 3),
        {**response.diagnostics, "constraints_applied": True,
         "born_projection_max_e": float(np.max(np.abs(response.born.mean(axis=0)))),
         "hessian_projection_relative_Frobenius": float(np.linalg.norm(
             h-response.force_constants.transpose(0,2,1,3).reshape(3*n,3*n))/max(np.linalg.norm(h), 1e-30))},
    )
