"""Compare the ABACUS/PYATB and VASP GaAs-nanowire response tensors."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np
from ase.io import read
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from plot_spectroscopy_across_dimensions import read_abacus_stru
from zstar.vasp_bec import parse_vasp_outcar


def source_record(path: Path) -> dict[str, object]:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {
        "name": path.name,
        "bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def read_zborn(path: Path) -> tuple[list[str], np.ndarray]:
    labels: list[str] = []
    tensors: list[np.ndarray] = []
    pattern = re.compile(r"^\s*\*?\s*\d+\s+([A-Za-z][A-Za-z0-9]*)\s+(.+)$")
    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if match is None:
            continue
        values = [float(value) for value in match.group(2).split()[:9]]
        if len(values) != 9:
            continue
        labels.append(match.group(1))
        tensors.append(np.asarray(values, dtype=float).reshape(3, 3))
    if not tensors:
        raise ValueError(f"No BEC tensors found in {path}")
    return labels, np.asarray(tensors)


def circular_center(values: np.ndarray) -> float:
    resultant = np.mean(np.exp(2j * np.pi * np.asarray(values, dtype=float)))
    if abs(resultant) < 1.0e-12:
        raise ValueError("Cannot determine a circular center for an open direction")
    return float((np.angle(resultant) / (2.0 * np.pi)) % 1.0)


def localized_coordinates(lattice: np.ndarray, fractional: np.ndarray) -> np.ndarray:
    lengths = np.linalg.norm(lattice, axis=1)
    if not np.allclose(lattice, np.diag(np.diag(lattice)), atol=1.0e-8):
        raise ValueError("GaAs benchmark comparison requires orthogonal cells")
    result = np.zeros_like(fractional, dtype=float)
    for axis in (0, 1):
        center = circular_center(fractional[:, axis])
        result[:, axis] = ((fractional[:, axis] - center + 0.5) % 1.0 - 0.5) * lengths[axis]
    result[:, 2] = fractional[:, 2] * lengths[2]
    return result


def match_atoms(
    first_labels: list[str],
    first_coordinates: np.ndarray,
    second_labels: list[str],
    second_coordinates: np.ndarray,
    period_z: float,
) -> np.ndarray:
    mapping = np.full(len(first_labels), -1, dtype=int)
    for label in sorted(set(first_labels)):
        first = np.asarray([i for i, item in enumerate(first_labels) if item == label])
        second = np.asarray([i for i, item in enumerate(second_labels) if item == label])
        if len(first) != len(second):
            raise ValueError(f"Species count differs for {label}: {len(first)} != {len(second)}")
        delta = first_coordinates[first, None, :] - second_coordinates[None, second, :]
        delta[:, :, 2] -= np.rint(delta[:, :, 2] / period_z) * period_z
        cost = np.linalg.norm(delta, axis=2)
        row, column = linear_sum_assignment(cost)
        if float(np.max(cost[row, column])) > 1.0e-5:
            raise ValueError(
                f"Coordinate match for {label} exceeds tolerance: {np.max(cost[row, column])}"
            )
        mapping[first[row]] = second[column]
    if np.any(mapping < 0):
        raise ValueError("Atom mapping is incomplete")
    return mapping


def line_polarizability(epsilon: np.ndarray, transverse_area_angstrom2: float) -> np.ndarray:
    return transverse_area_angstrom2 * (epsilon - np.eye(3)) / (4.0 * np.pi)


def response_quantity(path: Path, name: str) -> np.ndarray:
    data = json.loads(path.read_text(encoding="utf-8"))
    for quantity in data.get("quantities", []):
        if quantity.get("name") == name:
            return np.asarray(quantity["values"], dtype=float)
    raise ValueError(f"Quantity {name!r} is absent from {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--abacus-stru", type=Path, required=True)
    parser.add_argument("--abacus-zborn", type=Path, required=True)
    parser.add_argument("--abacus-response", type=Path, required=True)
    parser.add_argument("--vasp-poscar", type=Path, required=True)
    parser.add_argument("--vasp-outcar", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    abacus_lattice, abacus_labels, abacus_fractional = read_abacus_stru(args.abacus_stru)
    zborn_labels, abacus_bec = read_zborn(args.abacus_zborn)
    if zborn_labels != abacus_labels:
        raise ValueError("ABACUS STRU and Z-BORN-symm.out atom orders differ")
    abacus_coordinates = localized_coordinates(abacus_lattice, abacus_fractional)

    vasp_atoms = read(args.vasp_poscar, format="vasp")
    vasp_lattice = np.asarray(vasp_atoms.cell.array, dtype=float)
    vasp_fractional = np.asarray(vasp_atoms.get_scaled_positions(wrap=True), dtype=float)
    vasp_labels = list(vasp_atoms.get_chemical_symbols())
    vasp_coordinates = localized_coordinates(vasp_lattice, vasp_fractional)
    vasp_epsilon, vasp_bec = parse_vasp_outcar(args.vasp_outcar)

    mapping = match_atoms(
        abacus_labels,
        abacus_coordinates,
        vasp_labels,
        vasp_coordinates,
        float(np.linalg.norm(abacus_lattice[2])),
    )
    mapped_vasp_bec = vasp_bec[mapping]
    delta = abacus_bec - mapped_vasp_bec

    abacus_alpha = response_quantity(args.abacus_response, "line_polarizability")
    vasp_area = float(np.linalg.norm(np.cross(vasp_lattice[0], vasp_lattice[1])))
    vasp_alpha = line_polarizability(vasp_epsilon, vasp_area)
    alpha_delta = abacus_alpha - vasp_alpha
    longitudinal_delta = float(alpha_delta[2, 2])
    longitudinal_relative = abs(longitudinal_delta) / max(
        abs(float(vasp_alpha[2, 2])), np.finfo(float).tiny
    )

    species_metrics = {}
    for label in sorted(set(abacus_labels)):
        indices = [i for i, item in enumerate(abacus_labels) if item == label]
        block = delta[indices]
        species_metrics[label] = {
            "atoms": len(indices),
            "max_abs_e": float(np.max(np.abs(block))),
            "rms_e": float(np.sqrt(np.mean(block**2))),
        }

    report = {
        "schema_version": 1,
        "tensor_convention": "rows=displacement/force; columns=polarization/electric field",
        "atom_mapping_abacus_to_vasp_one_based": (mapping + 1).tolist(),
        "coordinate_match_max_angstrom": float(
            max(
                np.linalg.norm(
                    abacus_coordinates[i, :2] - vasp_coordinates[j, :2]
                )
                for i, j in enumerate(mapping)
            )
        ),
        "bec": {
            "max_abs_delta_e": float(np.max(np.abs(delta))),
            "rms_delta_e": float(np.sqrt(np.mean(delta**2))),
            "mean_abs_delta_e": float(np.mean(np.abs(delta))),
            "species": species_metrics,
            "abacus_acoustic_sum_e": np.sum(abacus_bec, axis=0).tolist(),
            "vasp_acoustic_sum_e": np.sum(mapped_vasp_bec, axis=0).tolist(),
            "representative_abacus_atom_1": {
                "label": abacus_labels[0],
                "vasp_atom_one_based": int(mapping[0] + 1),
                "abacus": abacus_bec[0].tolist(),
                "vasp": mapped_vasp_bec[0].tolist(),
                "delta": delta[0].tolist(),
            },
        },
        "electronic_response": {
            "comparison_scope": (
                "The periodic z component is the like-for-like 1D validation. "
                "VASP LEPSILON includes DFT local-field effects, whereas the "
                "PYATB Kubo response is independent-particle; transverse x/y "
                "differences are retained as a method-convention diagnostic."
            ),
            "abacus_line_polarizability_angstrom2": abacus_alpha.tolist(),
            "vasp_line_polarizability_angstrom2": vasp_alpha.tolist(),
            "delta_angstrom2": alpha_delta.tolist(),
            "longitudinal_z_abs_delta_angstrom2": abs(longitudinal_delta),
            "longitudinal_z_relative_delta": longitudinal_relative,
            "transverse_xy_max_abs_delta_angstrom2": float(
                np.max(np.abs(alpha_delta[:2, :2]))
            ),
            "max_abs_delta_angstrom2": float(np.max(np.abs(alpha_delta))),
        },
        "sources": {
            "abacus_stru": source_record(args.abacus_stru),
            "abacus_zborn": source_record(args.abacus_zborn),
            "abacus_response": source_record(args.abacus_response),
            "vasp_poscar": source_record(args.vasp_poscar),
            "vasp_outcar": source_record(args.vasp_outcar),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
