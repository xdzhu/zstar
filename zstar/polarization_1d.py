"""Hybrid one-dimensional polarization and Born-charge analysis.

For a wire aligned with Cartesian ``z``, the periodic polarization component
is evaluated from the PYATB Berry phase.  The two transverse components are
finite real-space dipoles integrated from neutral ABACUS charge-density cubes.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from .polarization_2d import (
    BOHR_M,
    ELEMENTARY_CHARGE,
    _direction_directory,
    _nearest_branch,
    _parse_pyatb_polarization,
    _pyatb_polar_file,
    _read_cube_profile_data,
    _read_pyatb_cell,
    _periodic_unwrap,
    _periodic_weighted_center,
    find_charge_cube,
    read_pyatb_lattice,
)


@dataclass(frozen=True)
class TransverseDipole:
    cube: str
    axis: str
    dipole_e_bohr: float
    period_bohr: float
    electron_count_raw: float
    ionic_charge: float
    neutrality_scale: float
    unwrap_center_bohr: float

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class Hybrid1DBorn:
    tensor: np.ndarray
    transverse_tensor: np.ndarray
    periodic_column: np.ndarray
    method: str
    displacement_angstrom: float
    periodic_axis: str
    diagnostics: dict


def integrate_transverse_dipole(
    cube_path: str | Path,
    axis: str,
    *,
    neutrality_tolerance: float = 0.05,
    ionic_valence_charges: Optional[Sequence[float]] = None,
    unwrap_center_bohr: Optional[float] = None,
) -> TransverseDipole:
    """Integrate one localized transverse dipole from a neutral cube."""

    axis_key = str(axis).lower()
    if axis_key not in {"x", "y"}:
        raise ValueError("transverse axis must be x or y for a z-periodic wire")
    axis_index = "xyz".index(axis_key)
    cube = _read_cube_profile_data(
        cube_path,
        neutrality_tolerance=neutrality_tolerance,
        ionic_valence_charges=ionic_valence_charges,
    )
    cell_vector = cube.cell[axis_index]
    period = float(np.linalg.norm(cell_vector))
    unit = cell_vector / period
    expected = np.eye(3)[axis_index]
    if not np.allclose(unit, expected, atol=1.0e-8):
        raise ValueError(
            "Hybrid 1D transverse dipoles currently require the nonperiodic "
            "cell vectors to align with positive Cartesian x and y"
        )
    for index, vector in enumerate(cube.cell):
        if index != axis_index and abs(float(vector @ unit)) > 1.0e-8 * period:
            raise ValueError("Hybrid 1D transverse dipoles require orthogonal cell axes")

    ionic_total = float(np.sum(cube.charges))
    atom_coordinate = cube.positions @ unit
    origin_coordinate = float(cube.origin @ unit)
    if unwrap_center_bohr is None:
        center = _periodic_weighted_center(
            atom_coordinate,
            cube.charges,
            origin_coordinate,
            period,
        )
    else:
        center = float(unwrap_center_bohr)
    center = origin_coordinate + ((center - origin_coordinate) % period)
    atom_unwrapped = _periodic_unwrap(atom_coordinate, center, period)

    grid_coordinate = (
        origin_coordinate
        + np.arange(cube.dimensions[axis_index], dtype=float)
        * float(cube.step_vectors[axis_index] @ unit)
    )
    grid_unwrapped = _periodic_unwrap(grid_coordinate, center, period)
    summed_axes = tuple(index for index in range(3) if index != axis_index)
    density_by_plane = np.sum(cube.density, axis=summed_axes)

    ionic_dipole = float(np.sum(cube.charges * atom_unwrapped))
    electronic_dipole = -float(
        np.sum(density_by_plane * grid_unwrapped)
        * cube.voxel_volume_bohr3
        * cube.neutrality_scale
    )
    return TransverseDipole(
        cube=str(cube.path),
        axis=axis_key,
        dipole_e_bohr=ionic_dipole + electronic_dipole,
        period_bohr=period,
        electron_count_raw=cube.electron_count_raw,
        ionic_charge=ionic_total,
        neutrality_scale=cube.neutrality_scale,
        unwrap_center_bohr=center,
    )


def _branch_delta(value: float, reference: float, period: float) -> float:
    delta = float(value) - float(reference)
    return delta - round(delta / float(period)) * float(period)


def calculate_hybrid_1d_born(
    atom_directory: str | Path,
    reference_directory: str | Path,
    *,
    method: str = "central",
    displacement_angstrom: float = 0.01,
    directions: Sequence[str] = ("x", "y", "z"),
    neutrality_tolerance: float = 0.05,
) -> Hybrid1DBorn:
    """Combine transverse cube dipoles with the periodic Berry derivative."""

    atom_dir = Path(atom_directory).resolve()
    reference = Path(reference_directory).resolve()
    method_key = method.lower()
    if method_key not in {"forward", "central"}:
        raise ValueError("method must be forward or central")

    reference_polar, reference_quanta = _parse_pyatb_polarization(
        _pyatb_polar_file(reference)
    )
    unit_vectors, volume_m3 = _read_pyatb_cell(
        reference / "pyatb" / "Out" / "input.json"
    )
    if not np.allclose(unit_vectors, np.eye(3), atol=1.0e-8):
        raise ValueError(
            "Hybrid 1D Born tensors currently require an orthogonal cell with "
            "the wire aligned with Cartesian z"
        )

    reference_dipoles = {
        axis: integrate_transverse_dipole(
            find_charge_cube(reference),
            axis,
            neutrality_tolerance=neutrality_tolerance,
        )
        for axis in ("x", "y")
    }
    displacement_m = float(displacement_angstrom) * 1.0e-10
    displacement_bohr = displacement_m / BOHR_M
    if not math.isfinite(displacement_bohr) or displacement_bohr <= 0.0:
        raise ValueError("displacement_angstrom must be finite and positive")

    direction_columns = {"x": 0, "y": 1, "z": 2}
    normalized_directions = [str(direction).lower() for direction in directions]
    invalid = [item for item in normalized_directions if item not in direction_columns]
    if invalid:
        raise ValueError(f"Unsupported displacement directions: {invalid}")
    if len(set(normalized_directions)) != len(normalized_directions):
        raise ValueError("Displacement directions must not contain duplicates")

    tensor = np.zeros((3, 3), dtype=float)
    diagnostics: dict = {
        "periodic_axis": "z",
        "transverse_axes": ["x", "y"],
        "reference_dipoles": {
            axis: value.to_dict() for axis, value in reference_dipoles.items()
        },
        "directions": {},
    }

    for direction in normalized_directions:
        beta = direction_columns[direction]
        plus_dir = _direction_directory(atom_dir, direction, "+")
        plus_polar, plus_quanta = _parse_pyatb_polarization(
            _pyatb_polar_file(plus_dir)
        )
        plus_dipoles = {
            axis: integrate_transverse_dipole(
                find_charge_cube(plus_dir),
                axis,
                neutrality_tolerance=neutrality_tolerance,
                unwrap_center_bohr=reference_dipoles[axis].unwrap_center_bohr,
            )
            for axis in ("x", "y")
        }

        if method_key == "central":
            minus_dir = _direction_directory(atom_dir, direction, "-")
            minus_polar, minus_quanta = _parse_pyatb_polarization(
                _pyatb_polar_file(minus_dir)
            )
            minus_dipoles = {
                axis: integrate_transverse_dipole(
                    find_charge_cube(minus_dir),
                    axis,
                    neutrality_tolerance=neutrality_tolerance,
                    unwrap_center_bohr=reference_dipoles[axis].unwrap_center_bohr,
                )
                for axis in ("x", "y")
            }
            quanta = np.where(np.abs(plus_quanta) > 0.0, plus_quanta, minus_quanta)
            delta_lattice = _nearest_branch(plus_polar, minus_polar, quanta)
            polar_derivative = (delta_lattice @ unit_vectors) / (
                2.0 * displacement_m
            )
            denominator = 2.0 * displacement_bohr
        else:
            minus_dipoles = None
            delta_lattice = _nearest_branch(
                plus_polar, reference_polar, reference_quanta
            )
            polar_derivative = (delta_lattice @ unit_vectors) / displacement_m
            denominator = displacement_bohr

        for alpha, axis in enumerate(("x", "y")):
            baseline = (
                minus_dipoles[axis].dipole_e_bohr
                if minus_dipoles is not None
                else reference_dipoles[axis].dipole_e_bohr
            )
            delta = _branch_delta(
                plus_dipoles[axis].dipole_e_bohr,
                baseline,
                plus_dipoles[axis].period_bohr,
            )
            tensor[beta, alpha] = delta / denominator

        tensor[beta, 2] = (
            volume_m3 / ELEMENTARY_CHARGE * polar_derivative[2]
        )
        diagnostics["directions"][direction] = {
            "berry_delta_C_per_m2": delta_lattice.tolist(),
            "periodic_component": float(tensor[beta, 2]),
            "plus_dipoles": {
                axis: value.to_dict() for axis, value in plus_dipoles.items()
            },
            "minus_dipoles": (
                None
                if minus_dipoles is None
                else {axis: value.to_dict() for axis, value in minus_dipoles.items()}
            ),
        }

    return Hybrid1DBorn(
        tensor=tensor,
        transverse_tensor=tensor[:, :2].copy(),
        periodic_column=tensor[:, 2].copy(),
        method=method_key,
        displacement_angstrom=float(displacement_angstrom),
        periodic_axis="z",
        diagnostics=diagnostics,
    )


def write_hybrid_1d_report(path: str | Path, result: Hybrid1DBorn) -> None:
    data = {
        "method": result.method,
        "displacement_angstrom": result.displacement_angstrom,
        "periodic_axis": result.periodic_axis,
        "tensor_convention": (
            "rows=atomic displacement/force; "
            "columns=polarization/electric field"
        ),
        "tensor": result.tensor.tolist(),
        "transverse_tensor": result.transverse_tensor.tolist(),
        "periodic_column": result.periodic_column.tolist(),
        "diagnostics": result.diagnostics,
    }
    Path(path).write_text(json.dumps(data, indent=2), encoding="utf-8")


def write_reference_1d_polarization(
    directory: str | Path = ".",
    output: str | Path = "zstar_1d_polarization.json",
    *,
    neutrality_tolerance: float = 0.05,
) -> Path:
    """Write transverse dipoles and periodic line polarization for one wire."""

    root = Path(directory).resolve()
    polar, quanta = _parse_pyatb_polarization(_pyatb_polar_file(root))
    lattice = read_pyatb_lattice(root / "pyatb" / "Out" / "input.json")
    length_angstrom = float(np.linalg.norm(lattice[2]))
    cross_section_angstrom2 = float(abs(np.linalg.det(lattice))) / length_angstrom
    cross_section_m2 = cross_section_angstrom2 * 1.0e-20
    dipoles = {
        axis: integrate_transverse_dipole(
            find_charge_cube(root),
            axis,
            neutrality_tolerance=neutrality_tolerance,
        )
        for axis in ("x", "y")
    }
    target = Path(output)
    if not target.is_absolute():
        target = root / target
    target.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "dimensionality": 1,
                "periodic_axis": "z",
                "transverse_dipole_e_angstrom": {
                    axis: value.dipole_e_bohr * BOHR_M / 1.0e-10
                    for axis, value in dipoles.items()
                },
                "periodic_line_polarization_C": float(polar[2] * cross_section_m2),
                "periodic_line_polarization_quantum_C": float(
                    quanta[2] * cross_section_m2
                ),
                "supercell_polarization_C_per_m2": polar.tolist(),
                "supercell_polarization_quanta_C_per_m2": quanta.tolist(),
                "cross_section_angstrom2": cross_section_angstrom2,
                "sources": {
                    "polarization": str(_pyatb_polar_file(root)),
                    "charge_cube": str(find_charge_cube(root)),
                },
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return target
