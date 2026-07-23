"""Hybrid 2D polarization and Born-charge analysis.

In-plane polarization is evaluated with the Berry-phase result from PYATB.
Out-of-plane polarization is evaluated as the real-space slab dipole from an
ABACUS charge-density cube, including the ionic contribution.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import re
from typing import Optional, Sequence

import numpy as np


BOHR_M = 5.29177210903e-11
ELEMENTARY_CHARGE = 1.602176634e-19
_FLOAT_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


@dataclass(frozen=True)
class SlabDipole:
    cube: str
    dipole_e_bohr: float
    polarization_C_per_m: float
    area_bohr2: float
    height_bohr: float
    electron_count_raw: float
    ionic_charge: float
    neutrality_scale: float
    normal: tuple[float, float, float]

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class Hybrid2DBorn:
    tensor: np.ndarray
    in_plane_tensor: np.ndarray
    out_of_plane_row: np.ndarray
    method: str
    displacement_angstrom: float
    diagnostics: dict


def _find_geometry_block(lines: Sequence[str]) -> int:
    for index in range(min(len(lines) - 3, 40)):
        rows = [lines[index + offset].split() for offset in range(4)]
        if any(len(row) < 4 for row in rows):
            continue
        try:
            int(rows[0][0])
            for row in rows[1:]:
                int(row[0])
            for row in rows:
                [float(value) for value in row[1:4]]
        except ValueError:
            continue
        if all(abs(int(row[0])) > 0 for row in rows[1:]):
            return index
    raise ValueError("Could not locate cube geometry header")


def find_charge_cube(directory: str | Path) -> Path:
    root = Path(directory)
    candidates: list[Path] = []
    patterns = (
        "OUT.*/SPIN*_CHG*.cube",
        "OUT.*/*CHG*.cube",
        "SPIN*_CHG*.cube",
        "*CHG*.cube",
    )
    for pattern in patterns:
        candidates.extend(root.glob(pattern))
    unique = sorted(
        {
            path.resolve()
            for path in candidates
            if path.is_file()
            and path.stat().st_size > 0
            and "elecstaticpot" not in path.name.lower()
        }
    )
    if not unique:
        raise FileNotFoundError(f"No ABACUS charge-density cube found under {root}")
    return unique[0]


def _periodic_unwrap(values: np.ndarray, center: float, period: float) -> np.ndarray:
    return center + (values - center) - np.round((values - center) / period) * period


def integrate_slab_dipole(
    cube_path: str | Path,
    *,
    neutrality_tolerance: float = 0.05,
) -> SlabDipole:
    """Integrate the total out-of-plane dipole of a neutral periodic slab."""

    path = Path(cube_path).resolve()
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    geometry_index = _find_geometry_block(lines)
    geometry = lines[geometry_index].split()
    natoms = abs(int(geometry[0]))
    origin = np.asarray([float(value) for value in geometry[1:4]], dtype=float)

    dimensions: list[int] = []
    steps: list[list[float]] = []
    for line in lines[geometry_index + 1 : geometry_index + 4]:
        fields = line.split()
        dimensions.append(abs(int(fields[0])))
        steps.append([float(value) for value in fields[1:4]])
    nx, ny, nz = dimensions
    step_vectors = np.asarray(steps, dtype=float)
    cell = step_vectors * np.asarray(dimensions, dtype=float)[:, None]
    area = float(np.linalg.norm(np.cross(cell[0], cell[1])))
    volume = float(abs(np.linalg.det(cell)))
    if area <= 0.0 or volume <= 0.0:
        raise ValueError("Cube cell has zero area or volume")
    normal_vector = np.cross(cell[0], cell[1])
    normal = normal_vector / np.linalg.norm(normal_vector)
    if float(np.dot(normal, cell[2])) < 0.0:
        normal *= -1.0
    height = volume / area

    atom_start = geometry_index + 4
    atom_stop = atom_start + natoms
    ionic_charges: list[float] = []
    atom_positions: list[list[float]] = []
    for line in lines[atom_start:atom_stop]:
        fields = line.split()
        if len(fields) < 5:
            raise ValueError(f"Malformed cube atom line: {line}")
        ionic_charges.append(float(fields[1]))
        atom_positions.append([float(value) for value in fields[2:5]])
    charges = np.asarray(ionic_charges, dtype=float)
    positions = np.asarray(atom_positions, dtype=float)
    ionic_total = float(np.sum(charges))
    if ionic_total <= 0.0:
        raise ValueError("Cube header does not contain positive ionic valence charges")

    values = np.fromiter(
        (
            float(token)
            for line in lines[atom_stop:]
            for token in line.split()
        ),
        dtype=float,
    )
    expected = nx * ny * nz
    if len(values) != expected:
        raise ValueError(
            f"Cube data size mismatch: expected {expected}, found {len(values)}"
        )
    density = values.reshape(nx, ny, nz)
    voxel_volume = float(abs(np.linalg.det(step_vectors)))
    electron_count = float(np.sum(density) * voxel_volume)
    if electron_count <= 0.0:
        raise ValueError("Integrated cube electron density is not positive")
    relative_error = abs(electron_count - ionic_total) / ionic_total
    if relative_error > float(neutrality_tolerance):
        raise ValueError(
            "Charge cube fails neutrality check: "
            f"electrons={electron_count:.8f}, ionic={ionic_total:.8f}, "
            f"relative error={relative_error:.3%}"
        )
    neutrality_scale = ionic_total / electron_count

    atom_normal = positions @ normal
    center = float(np.sum(charges * atom_normal) / ionic_total)
    origin_normal = float(origin @ normal)
    center = origin_normal + ((center - origin_normal) % height)
    atom_unwrapped = _periodic_unwrap(atom_normal, center, height)

    grid_normal = (
        origin_normal
        + np.arange(nz, dtype=float) * float(step_vectors[2] @ normal)
    )
    grid_unwrapped = _periodic_unwrap(grid_normal, center, height)
    density_by_plane = np.sum(density, axis=(0, 1))

    ionic_dipole = float(np.sum(charges * atom_unwrapped))
    electronic_dipole = -float(
        np.sum(density_by_plane * grid_unwrapped)
        * voxel_volume
        * neutrality_scale
    )
    dipole = ionic_dipole + electronic_dipole
    polarization = dipole / area * ELEMENTARY_CHARGE / BOHR_M
    return SlabDipole(
        cube=str(path),
        dipole_e_bohr=dipole,
        polarization_C_per_m=polarization,
        area_bohr2=area,
        height_bohr=height,
        electron_count_raw=electron_count,
        ionic_charge=ionic_total,
        neutrality_scale=neutrality_scale,
        normal=tuple(float(value) for value in normal),
    )


def _parse_pyatb_polarization(path: Path) -> tuple[np.ndarray, np.ndarray]:
    pattern = re.compile(
        rf"direction is in\s+([abc]),\s*P\s*=\s*({_FLOAT_RE})\s*"
        rf"\(mod\s*({_FLOAT_RE})\)\s*C/m\^2",
        re.IGNORECASE,
    )
    values = {"a": None, "b": None, "c": None}
    quanta = {"a": None, "b": None, "c": None}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        match = pattern.search(line)
        if match:
            values[match.group(1).lower()] = float(match.group(2))
            quanta[match.group(1).lower()] = float(match.group(3))
    if any(value is None for value in values.values()):
        raise ValueError(f"Incomplete PYATB polarization output: {path}")
    return (
        np.asarray([values[key] for key in ("a", "b", "c")], dtype=float),
        np.asarray([quanta[key] for key in ("a", "b", "c")], dtype=float),
    )


def _read_pyatb_cell(input_json: Path) -> tuple[np.ndarray, float]:
    data = json.loads(input_json.read_text(encoding="utf-8"))
    candidates = (
        data.get("lattice_vector"),
        (data.get("LATTICE") or {}).get("lattice_vector"),
        (data.get("lattice") or {}).get("lattice_vector"),
    )
    lattice = next((item for item in candidates if item is not None), None)
    if lattice is None:
        def search(value):
            if isinstance(value, dict):
                for key, child in value.items():
                    if key.lower() in {"lattice_vector", "lattice_vectors"}:
                        return child
                    found = search(child)
                    if found is not None:
                        return found
            return None

        lattice = search(data)
    matrix_angstrom = np.asarray(lattice, dtype=float)
    if matrix_angstrom.shape != (3, 3):
        raise ValueError(f"Invalid lattice in {input_json}")
    unit_vectors = matrix_angstrom / np.linalg.norm(
        matrix_angstrom, axis=1, keepdims=True
    )
    volume_m3 = float(abs(np.linalg.det(matrix_angstrom))) * 1.0e-30
    return unit_vectors, volume_m3


def _nearest_branch(target: np.ndarray, reference: np.ndarray, quanta: np.ndarray) -> np.ndarray:
    result = target - reference
    for index, quantum in enumerate(quanta):
        if abs(float(quantum)) > np.finfo(float).tiny:
            result[index] -= round(result[index] / quantum) * quantum
    return result


def _direction_directory(atom_dir: Path, direction: str, sign: str) -> Path:
    candidate = atom_dir / f"{direction}{sign}"
    if candidate.is_dir():
        return candidate
    if sign == "+" and (atom_dir / direction).is_dir():
        return atom_dir / direction
    raise FileNotFoundError(
        f"Displacement directory not found for {direction}{sign} under {atom_dir}"
    )


def _pyatb_polar_file(directory: Path) -> Path:
    path = directory / "pyatb" / "Out" / "Polarization" / "polarization.dat"
    if not path.is_file():
        raise FileNotFoundError(f"PYATB polarization output not found: {path}")
    return path


def calculate_hybrid_2d_born(
    atom_directory: str | Path,
    reference_directory: str | Path,
    *,
    method: str = "central",
    displacement_angstrom: float = 0.01,
    directions: Sequence[str] = ("x", "y", "z"),
    neutrality_tolerance: float = 0.05,
) -> Hybrid2DBorn:
    """Combine in-plane Berry derivatives with out-of-plane cube derivatives."""

    atom_dir = Path(atom_directory).resolve()
    reference = Path(reference_directory).resolve()
    method_key = method.lower()
    if method_key not in {"forward", "central"}:
        raise ValueError("method must be forward or central")

    reference_polar, reference_quanta = _parse_pyatb_polarization(
        _pyatb_polar_file(reference)
    )
    cell_json = reference / "pyatb" / "Out" / "input.json"
    unit_vectors, volume_m3 = _read_pyatb_cell(cell_json)
    reference_dipole = integrate_slab_dipole(
        find_charge_cube(reference), neutrality_tolerance=neutrality_tolerance
    )
    normal = np.asarray(reference_dipole.normal, dtype=float)
    if abs(abs(float(normal[2])) - 1.0) > 1.0e-6:
        raise ValueError(
            "Hybrid 2D Born tensors currently require the slab normal to be "
            "aligned with Cartesian z"
        )

    displacement_m = float(displacement_angstrom) * 1.0e-10
    displacement_bohr = float(displacement_angstrom) * 1.0e-10 / BOHR_M
    tensor = np.zeros((3, 3), dtype=float)
    diagnostics: dict = {"reference_dipole": reference_dipole.to_dict(), "directions": {}}

    direction_columns = {"x": 0, "y": 1, "z": 2}
    normalized_directions = [str(direction).lower() for direction in directions]
    invalid = [item for item in normalized_directions if item not in direction_columns]
    if invalid:
        raise ValueError(f"Unsupported displacement directions: {invalid}")
    if len(set(normalized_directions)) != len(normalized_directions):
        raise ValueError("Displacement directions must not contain duplicates")

    for direction in normalized_directions:
        beta = direction_columns[direction]
        plus_dir = _direction_directory(atom_dir, direction, "+")
        plus_polar, plus_quanta = _parse_pyatb_polarization(
            _pyatb_polar_file(plus_dir)
        )
        plus_dipole = integrate_slab_dipole(
            find_charge_cube(plus_dir), neutrality_tolerance=neutrality_tolerance
        )
        if method_key == "central":
            minus_dir = _direction_directory(atom_dir, direction, "-")
            minus_polar, minus_quanta = _parse_pyatb_polarization(
                _pyatb_polar_file(minus_dir)
            )
            minus_dipole = integrate_slab_dipole(
                find_charge_cube(minus_dir), neutrality_tolerance=neutrality_tolerance
            )
            quanta = np.where(
                np.abs(plus_quanta) > 0.0, plus_quanta, minus_quanta
            )
            delta_lattice = _nearest_branch(plus_polar, minus_polar, quanta)
            polar_derivative = (delta_lattice @ unit_vectors) / (
                2.0 * displacement_m
            )
            dipole_delta = plus_dipole.dipole_e_bohr - minus_dipole.dipole_e_bohr
            dipole_delta -= round(
                dipole_delta / reference_dipole.height_bohr
            ) * reference_dipole.height_bohr
            z_born = dipole_delta / (2.0 * displacement_bohr)
            minus_report: Optional[dict] = minus_dipole.to_dict()
        else:
            delta_lattice = _nearest_branch(
                plus_polar, reference_polar, reference_quanta
            )
            polar_derivative = (delta_lattice @ unit_vectors) / displacement_m
            dipole_delta = (
                plus_dipole.dipole_e_bohr - reference_dipole.dipole_e_bohr
            )
            dipole_delta -= round(
                dipole_delta / reference_dipole.height_bohr
            ) * reference_dipole.height_bohr
            z_born = dipole_delta / displacement_bohr
            minus_report = None

        # PYATB P is a 3D polarization. Multiplication by cell volume gives the
        # vacuum-independent in-plane Born response.
        in_plane = volume_m3 / ELEMENTARY_CHARGE * polar_derivative[:2]
        tensor[:2, beta] = in_plane
        tensor[2, beta] = z_born
        diagnostics["directions"][direction] = {
            "plus": plus_dipole.to_dict(),
            "minus": minus_report,
            "berry_delta_C_per_m2": delta_lattice.tolist(),
            "in_plane_column": in_plane.tolist(),
            "out_of_plane_component": float(z_born),
        }

    return Hybrid2DBorn(
        tensor=tensor,
        in_plane_tensor=tensor[:2].copy(),
        out_of_plane_row=tensor[2].copy(),
        method=method_key,
        displacement_angstrom=float(displacement_angstrom),
        diagnostics=diagnostics,
    )


def write_hybrid_2d_report(
    path: str | Path,
    result: Hybrid2DBorn,
) -> None:
    data = {
        "method": result.method,
        "displacement_angstrom": result.displacement_angstrom,
        "tensor": result.tensor.tolist(),
        "in_plane_tensor": result.in_plane_tensor.tolist(),
        "out_of_plane_row": result.out_of_plane_row.tolist(),
        "diagnostics": result.diagnostics,
    }
    Path(path).write_text(json.dumps(data, indent=2), encoding="utf-8")
