"""Molecular atomic polar tensors from ABACUS/PYATB finite differences."""

from __future__ import annotations

import json
from pathlib import Path
import re
from typing import Sequence

import numpy as np

from .polarization_2d import read_pyatb_lattice
from .spectra import ELEMENTARY_CHARGE, read_pyatb_polarization
from .stru_analyzer import stru_analyzer


_ATOM_DIR_RE = re.compile(r"^([1-9]\d*)\.([^.]+)$")


def _direction_directory(atom_dir: Path, direction: str, sign: str) -> Path:
    candidate = atom_dir / f"{direction}{sign}"
    if candidate.is_dir():
        return candidate
    if sign == "+" and (atom_dir / direction).is_dir():
        return atom_dir / direction
    raise FileNotFoundError(
        f"Displacement directory not found for {direction}{sign} under {atom_dir}"
    )


def _stage_polarization(stage: Path) -> tuple[np.ndarray, np.ndarray, Path]:
    errors: list[str] = []
    for candidate in (stage / "pyatb", stage / "pyatb-polar", stage):
        try:
            return read_pyatb_polarization(candidate)
        except FileNotFoundError as exc:
            errors.append(str(exc))
    raise FileNotFoundError("; ".join(errors))


def _nearest_branch(delta: np.ndarray, quanta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    adjusted = np.asarray(delta, dtype=float).copy()
    shifts = np.zeros(3, dtype=int)
    valid = np.abs(quanta) > np.finfo(float).tiny
    shifts[valid] = np.rint(adjusted[valid] / quanta[valid]).astype(int)
    adjusted[valid] -= shifts[valid] * quanta[valid]
    return adjusted, shifts


def calculate_molecular_apt(
    atom_directory: str | Path,
    reference_directory: str | Path,
    *,
    method: str = "central",
    displacement_angstrom: float = 0.01,
    directions: Sequence[str] = ("x", "y", "z"),
) -> tuple[np.ndarray, dict]:
    """Calculate one molecular APT in displacement-by-dipole convention."""

    atom_dir = Path(atom_directory).resolve()
    reference = Path(reference_directory).resolve()
    method_key = method.lower()
    if method_key not in {"forward", "central"}:
        raise ValueError("method must be forward or central")
    displacement_m = float(displacement_angstrom) * 1.0e-10
    if displacement_m <= 0.0:
        raise ValueError("displacement_angstrom must be positive")

    reference_polarization, reference_quanta, reference_source = _stage_polarization(
        reference
    )
    lattice_path = reference / "pyatb" / "Out" / "input.json"
    lattice = read_pyatb_lattice(lattice_path)
    volume_m3 = abs(float(np.linalg.det(lattice))) * 1.0e-30
    unit_vectors = lattice / np.linalg.norm(lattice, axis=1, keepdims=True)

    tensor = np.zeros((3, 3), dtype=float)
    diagnostics = {
        "method": method_key,
        "displacement_angstrom": float(displacement_angstrom),
        "cell_volume_angstrom3": abs(float(np.linalg.det(lattice))),
        "reference_polarization": str(reference_source),
        "directions": {},
    }
    normalized_directions = [str(direction).lower() for direction in directions]
    if set(normalized_directions) - {"x", "y", "z"}:
        raise ValueError("directions must contain only x, y, and z")

    for direction in normalized_directions:
        row = "xyz".index(direction)
        plus, plus_quanta, plus_source = _stage_polarization(
            _direction_directory(atom_dir, direction, "+")
        )
        if method_key == "central":
            minus, minus_quanta, minus_source = _stage_polarization(
                _direction_directory(atom_dir, direction, "-")
            )
            delta = plus - minus
            quanta = 0.5 * (np.abs(plus_quanta) + np.abs(minus_quanta))
            denominator_m = 2.0 * displacement_m
        else:
            minus = reference_polarization
            minus_quanta = reference_quanta
            minus_source = reference_source
            delta = plus - reference_polarization
            quanta = 0.5 * (np.abs(plus_quanta) + np.abs(reference_quanta))
            denominator_m = displacement_m

        delta, shifts = _nearest_branch(delta, quanta)
        delta_cartesian = delta @ unit_vectors
        tensor[row, :] = delta_cartesian * volume_m3 / (
            denominator_m * ELEMENTARY_CHARGE
        )
        diagnostics["directions"][direction] = {
            "plus": str(plus_source),
            "minus": str(minus_source),
            "delta_polarization_C_m2": delta.tolist(),
            "branch_shifts": shifts.tolist(),
            "plus_quanta_C_m2": plus_quanta.tolist(),
            "minus_quanta_C_m2": minus_quanta.tolist(),
        }

    return tensor, diagnostics


def _write_zborn(path: Path, atoms: list[dict]) -> None:
    names = ("xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz")
    lines = [f"{'No. Atom': <8} " + " ".join(f"{name:>14}" for name in names) + "\n"]
    for atom in atoms:
        values = " ".join(f"{value: 14.8f}" for value in np.asarray(atom["tensor"]).reshape(9))
        lines.append(f"*{int(atom['index']):>4} {str(atom['label']):<3} {values}\n")
    path.write_text("".join(lines), encoding="utf-8", newline="\n")


def _read_zborn(path: Path) -> list[dict]:
    atoms: list[dict] = []
    for raw in path.read_text(encoding="utf-8").splitlines()[1:]:
        fields = raw.replace("*", " ", 1).split()
        if len(fields) < 11:
            continue
        atoms.append(
            {
                "index": int(fields[0]),
                "label": fields[1],
                "tensor": np.asarray([float(value) for value in fields[2:11]])
                .reshape(3, 3)
                .tolist(),
            }
        )
    if not atoms:
        raise ValueError(f"No molecular APT tensors found in {path}")
    return atoms


def _natoms_from_stru(path: Path) -> int:
    structure = stru_analyzer(str(path))
    coordinates = structure[5]
    return sum(len(values) for values in coordinates.values())


def collect_molecular_apts(
    root: str | Path = ".",
    *,
    method: str = "central",
    displacement_angstrom: float = 0.01,
    symprec: float = 1.0e-3,
    response_output: str | Path = "zstar_response.json",
) -> dict:
    """Collect molecular APTs, reconstruct symmetry, and enforce translation invariance."""

    root_path = Path(root).resolve()
    reference = root_path / "0.no-move"
    stru = reference / "STRU"
    if not stru.is_file():
        raise FileNotFoundError(stru)

    atoms: list[dict] = []
    diagnostics: list[dict] = []
    for atom_dir in sorted(root_path.iterdir(), key=lambda path: path.name):
        if not atom_dir.is_dir():
            continue
        match = _ATOM_DIR_RE.match(atom_dir.name)
        if match is None:
            continue
        tensor, atom_diagnostics = calculate_molecular_apt(
            atom_dir,
            reference,
            method=method,
            displacement_angstrom=displacement_angstrom,
        )
        atom = {
            "index": int(match.group(1)),
            "label": match.group(2),
            "tensor": tensor.tolist(),
        }
        atoms.append(atom)
        diagnostics.append({**atom_diagnostics, "index": atom["index"], "label": atom["label"]})
    atoms.sort(key=lambda item: item["index"])
    if not atoms:
        raise ValueError(f"No molecular atom-displacement folders found under {root_path}")

    raw_path = root_path / "Z-BORN-reduced.out"
    _write_zborn(raw_path, atoms)
    natoms_total = _natoms_from_stru(stru)
    all_path = root_path / "Z-BORN-all.out"
    if len(atoms) == natoms_total:
        _write_zborn(all_path, atoms)

    from .verify_born_symmetry import run_symcheck

    kwargs = {
        "stru": str(stru),
        "reduced": str(raw_path),
        "symprec": float(symprec),
        "out": str(root_path / "molecular_apt_symmetry_report.txt"),
        "json_path": str(root_path / "molecular_apt_symmetry_report.json"),
        "csv_path": None,
        "symm_out": str(root_path / "Z-BORN-symm.out"),
    }
    if all_path.is_file():
        kwargs["all"] = str(all_path)
    run_symcheck(**kwargs)
    symm_path = root_path / "Z-BORN-symm.out"
    corrected_atoms = _read_zborn(symm_path)
    symmetry_report_path = root_path / "molecular_apt_symmetry_report.json"
    symmetry_report = json.loads(symmetry_report_path.read_text(encoding="utf-8"))
    symmetry_born = symmetry_report["symmetry_born"]
    raw_expanded = symmetry_born["Z_symmetry_mean"]
    corrected = symmetry_born["Z_corrected"]
    for atom in corrected_atoms:
        atom["tensor"] = corrected[str(atom["index"])]
        atom["gapt"] = float(np.trace(np.asarray(atom["tensor"])) / 3.0)
    for atom in atoms:
        atom["gapt"] = float(np.trace(np.asarray(atom["tensor"])) / 3.0)
    raw_reduced_sum = np.sum(
        [np.asarray(atom["tensor"]) for atom in atoms], axis=0
    )
    raw_sum = np.sum(
        [np.asarray(raw_expanded[key]) for key in sorted(raw_expanded, key=int)],
        axis=0,
    )
    corrected_sum = np.sum(
        [np.asarray(atom["tensor"]) for atom in corrected_atoms], axis=0
    )

    result = {
        "schema_version": 1,
        "backend": "abacus-pyatb",
        "dimensionality": 0,
        "quantity": "atomic_polar_tensor",
        "method": method.lower(),
        "displacement_angstrom": float(displacement_angstrom),
        "natoms_total": natoms_total,
        "natoms_calculated": len(atoms),
        "tensor_convention": (
            "rows=atomic displacement; columns=molecular dipole; "
            "molecular atomic polar tensor in units of e"
        ),
        "raw_atoms": atoms,
        "atoms": corrected_atoms,
        "raw_reduced_sum_tensor": raw_reduced_sum.tolist(),
        "raw_acoustic_sum_tensor": raw_sum.tolist(),
        "acoustic_sum_tensor": corrected_sum.tolist(),
        "diagnostics": diagnostics,
        "files": {
            "raw_reduced": str(raw_path),
            "raw_all": str(all_path) if all_path.is_file() else None,
            "symmetry_corrected": str(symm_path),
            "symmetry_report": str(symmetry_report_path),
        },
    }
    json_path = root_path / "molecular_apt.json"
    json_path.write_text(json.dumps(result, indent=2), encoding="utf-8", newline="\n")

    from .response_schema import response_record_from_bec_result

    response_path = Path(response_output)
    if not response_path.is_absolute():
        response_path = root_path / response_path
    response_record_from_bec_result(
        result,
        dimensionality=0,
        provenance={
            "collector": "zstar.molecular_bec.collect_molecular_apts",
            "source": str(root_path),
            "legacy_result": str(json_path),
        },
    ).write(response_path)
    result["json_output"] = str(json_path)
    result["response_output"] = str(response_path)
    return result
