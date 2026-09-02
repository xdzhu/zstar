"""Hybrid 2D polarization and Born-charge analysis.

In-plane polarization is evaluated with the Berry-phase result from PYATB.
Out-of-plane polarization is evaluated as the real-space slab dipole from an
ABACUS charge-density cube, including the ionic contribution.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import csv
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
class MolecularDipole:
    """Three-dimensional dipole reconstructed from a neutral charge cube."""

    cube: str
    dipole_e_bohr: tuple[float, float, float]
    dipole_debye: tuple[float, float, float]
    electron_count_raw: float
    ionic_charge: float
    neutrality_scale: float
    center_fractional: tuple[float, float, float]
    diagnostics: dict

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class Hybrid2DBorn:
    tensor: np.ndarray
    in_plane_tensor: np.ndarray
    out_of_plane_column: np.ndarray
    method: str
    displacement_angstrom: float
    diagnostics: dict


@dataclass(frozen=True)
class SlabChargeDifference:
    reference_cube: str
    displaced_cube: str
    coordinate_angstrom: np.ndarray
    reference_electron_line_density_e_per_angstrom: np.ndarray
    displaced_electron_line_density_e_per_angstrom: np.ndarray
    electron_charge_difference_e_per_angstrom: np.ndarray
    cumulative_charge_difference_e: np.ndarray
    cumulative_dipole_difference_e_angstrom: np.ndarray
    reference_ion_positions_angstrom: np.ndarray
    displaced_ion_positions_angstrom: np.ndarray
    total_dipole_change_e_angstrom: float
    sheet_polarization_change_C_per_m: float
    displacement_angstrom: Optional[float]
    effective_charge_e: Optional[float]
    diagnostics: dict


@dataclass(frozen=True)
class _CubeProfileData:
    path: Path
    dimensions: tuple[int, int, int]
    origin: np.ndarray
    step_vectors: np.ndarray
    cell: np.ndarray
    normal: np.ndarray
    area_bohr2: float
    height_bohr: float
    charges: np.ndarray
    positions: np.ndarray
    density: np.ndarray
    voxel_volume_bohr3: float
    electron_count_raw: float
    neutrality_scale: float


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
        "*charge*density*.cube",
        "*ELECTRON_DENSITY*.cube",
        "*.cube",
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
            and "potential" not in path.name.lower()
        }
    )
    if not unique:
        raise FileNotFoundError(f"No ABACUS charge-density cube found under {root}")
    return unique[0]


def _periodic_unwrap(values: np.ndarray, center: float, period: float) -> np.ndarray:
    return center + (values - center) - np.round((values - center) / period) * period


def _periodic_weighted_center(
    values: np.ndarray,
    weights: np.ndarray,
    origin: float,
    period: float,
) -> float:
    """Return a weighted circular center for coordinates in one periodic cell."""

    coordinates = np.asarray(values, dtype=float)
    charges = np.asarray(weights, dtype=float)
    if coordinates.shape != charges.shape or coordinates.ndim != 1:
        raise ValueError("periodic center values and weights must be matching 1D arrays")
    total = float(np.sum(charges))
    if total <= 0.0 or period <= 0.0:
        raise ValueError("periodic center weights and period must be positive")
    phases = 2.0 * math.pi * (coordinates - float(origin)) / float(period)
    resultant = np.sum(charges * np.exp(1j * phases)) / total
    if abs(resultant) < 1.0e-10:
        # A delocalized or exactly inversion-balanced distribution has no
        # unique circular center. Preserve the historical arithmetic fallback.
        center = float(np.sum(charges * coordinates) / total)
    else:
        phase = float(np.angle(resultant)) % (2.0 * math.pi)
        center = float(origin) + phase * float(period) / (2.0 * math.pi)
    return float(origin) + ((center - float(origin)) % float(period))


def _cube_sidecar_charges(path: Path, natoms: int) -> np.ndarray | None:
    sidecar = path.with_suffix(path.suffix + ".zstar.json")
    if not sidecar.is_file():
        return None
    data = json.loads(sidecar.read_text(encoding="utf-8"))
    if data.get("schema") != "zstar-cube":
        raise ValueError(f"Invalid ZStar cube sidecar schema: {sidecar}")
    if data.get("density_sign", "positive_electron_density") != "positive_electron_density":
        raise ValueError(f"Unsupported cube density sign in {sidecar}")
    charges = np.asarray(data.get("ionic_valence_charges", ()), dtype=float)
    if charges.shape != (natoms,) or np.any(charges <= 0.0):
        raise ValueError(
            f"Cube sidecar ionic charges must have shape {(natoms,)}: {sidecar}"
        )
    return charges


def _read_cube_profile_data(
    cube_path: str | Path,
    *,
    neutrality_tolerance: float,
    ionic_valence_charges: Optional[Sequence[float]] = None,
) -> _CubeProfileData:
    if (
        not math.isfinite(float(neutrality_tolerance))
        or float(neutrality_tolerance) < 0.0
    ):
        raise ValueError("neutrality_tolerance must be finite and non-negative")
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
    normal = np.cross(cell[0], cell[1])
    normal /= np.linalg.norm(normal)
    if float(np.dot(normal, cell[2])) < 0.0:
        normal *= -1.0
    height = volume / area

    atom_start = geometry_index + 4
    atom_stop = atom_start + natoms
    charges: list[float] = []
    positions: list[list[float]] = []
    for line in lines[atom_start:atom_stop]:
        fields = line.split()
        if len(fields) < 5:
            raise ValueError(f"Malformed cube atom line: {line}")
        charges.append(float(fields[1]))
        positions.append([float(value) for value in fields[2:5]])
    ionic_charges = np.asarray(charges, dtype=float)
    sidecar_charges = _cube_sidecar_charges(path, natoms)
    if ionic_valence_charges is not None:
        supplied = np.asarray(ionic_valence_charges, dtype=float)
        if supplied.shape != (natoms,) or np.any(supplied <= 0.0):
            raise ValueError(
                f"ionic_valence_charges must have shape {(natoms,)} and be positive"
            )
        ionic_charges = supplied
    elif sidecar_charges is not None:
        ionic_charges = sidecar_charges
    ionic_total = float(np.sum(ionic_charges))
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

    return _CubeProfileData(
        path=path,
        dimensions=(nx, ny, nz),
        origin=origin,
        step_vectors=step_vectors,
        cell=cell,
        normal=normal,
        area_bohr2=area,
        height_bohr=height,
        charges=ionic_charges,
        positions=np.asarray(positions, dtype=float),
        density=density,
        voxel_volume_bohr3=voxel_volume,
        electron_count_raw=electron_count,
        neutrality_scale=ionic_total / electron_count,
    )


def integrate_slab_dipole(
    cube_path: str | Path,
    *,
    neutrality_tolerance: float = 0.05,
    ionic_valence_charges: Optional[Sequence[float]] = None,
) -> SlabDipole:
    """Integrate a neutral slab dipole from any standard electron-density cube."""

    cube = _read_cube_profile_data(
        cube_path,
        neutrality_tolerance=neutrality_tolerance,
        ionic_valence_charges=ionic_valence_charges,
    )
    charges = cube.charges
    ionic_total = float(np.sum(charges))
    atom_normal = cube.positions @ cube.normal
    origin_normal = float(cube.origin @ cube.normal)
    center = _periodic_weighted_center(
        atom_normal,
        charges,
        origin_normal,
        cube.height_bohr,
    )
    atom_unwrapped = _periodic_unwrap(
        atom_normal, center, cube.height_bohr
    )

    grid_normal = (
        origin_normal
        + np.arange(cube.dimensions[2], dtype=float)
        * float(cube.step_vectors[2] @ cube.normal)
    )
    grid_unwrapped = _periodic_unwrap(
        grid_normal, center, cube.height_bohr
    )
    density_by_plane = np.sum(cube.density, axis=(0, 1))

    ionic_dipole = float(np.sum(charges * atom_unwrapped))
    electronic_dipole = -float(
        np.sum(density_by_plane * grid_unwrapped)
        * cube.voxel_volume_bohr3
        * cube.neutrality_scale
    )
    dipole = ionic_dipole + electronic_dipole
    polarization = dipole / cube.area_bohr2 * ELEMENTARY_CHARGE / BOHR_M
    return SlabDipole(
        cube=str(cube.path),
        dipole_e_bohr=dipole,
        polarization_C_per_m=polarization,
        area_bohr2=cube.area_bohr2,
        height_bohr=cube.height_bohr,
        electron_count_raw=cube.electron_count_raw,
        ionic_charge=ionic_total,
        neutrality_scale=cube.neutrality_scale,
        normal=tuple(float(value) for value in cube.normal),
    )


def integrate_molecular_dipole(
    cube_path: str | Path,
    *,
    neutrality_tolerance: float = 0.05,
    ionic_valence_charges: Optional[Sequence[float]] = None,
) -> MolecularDipole:
    """Integrate a neutral molecular dipole from a 3D electron-density cube.

    The cube is interpreted as a periodic numerical box with localized molecular
    charge. Atomic positions and the density grid are unwrapped around a
    charge-weighted molecular center before integrating, so a molecule crossing
    a box boundary does not create a spurious dipole jump. The returned dipole
    is in ``e bohr`` and Debye; it is origin independent when the neutrality
    check passes.
    """

    cube = _read_cube_profile_data(
        cube_path,
        neutrality_tolerance=neutrality_tolerance,
        ionic_valence_charges=ionic_valence_charges,
    )
    cell = np.asarray(cube.cell, dtype=float)
    inverse_cell = np.linalg.inv(cell)
    origin = np.asarray(cube.origin, dtype=float)
    atom_fractional = (cube.positions - origin) @ inverse_cell
    center = np.asarray(
        [
            _periodic_weighted_center(
                atom_fractional[:, axis],
                cube.charges,
                0.0,
                1.0,
            )
            for axis in range(3)
        ],
        dtype=float,
    )
    atom_unwrapped_fractional = np.column_stack(
        [
            _periodic_unwrap(atom_fractional[:, axis], center[axis], 1.0)
            for axis in range(3)
        ]
    )

    # The density array can be large (ABACUS HSE cubes are often 243^3).
    # Compute its first moments using one-dimensional coordinate vectors rather
    # than materializing three full Cartesian coordinate arrays.
    fractional_axes = [
        _periodic_unwrap(
            np.arange(size, dtype=float) / float(size), center[axis], 1.0
        )
        for axis, size in enumerate(cube.dimensions)
    ]
    density = cube.density * cube.voxel_volume_bohr3 * cube.neutrality_scale
    electronic_fractional_moment = np.asarray(
        [
            float(
                np.sum(
                    density
                    * fractional_axes[axis].reshape(
                        (-1, 1, 1) if axis == 0 else
                        (1, -1, 1) if axis == 1 else (1, 1, -1)
                    )
                )
            )
            for axis in range(3)
        ],
        dtype=float,
    )
    ionic_fractional_moment = np.sum(
        cube.charges[:, None] * atom_unwrapped_fractional, axis=0
    )
    ionic_dipole = ionic_fractional_moment @ cell
    electronic_dipole = -electronic_fractional_moment @ cell
    dipole = ionic_dipole + electronic_dipole
    dipole_debye = dipole * BOHR_M * ELEMENTARY_CHARGE / 3.33564e-30
    return MolecularDipole(
        cube=str(cube.path),
        dipole_e_bohr=tuple(float(value) for value in dipole),
        dipole_debye=tuple(float(value) for value in dipole_debye),
        electron_count_raw=cube.electron_count_raw,
        ionic_charge=float(np.sum(cube.charges)),
        neutrality_scale=cube.neutrality_scale,
        center_fractional=tuple(float(value) for value in center),
        diagnostics={
            "cell_bohr": cell.tolist(),
            "dimensions": list(cube.dimensions),
            "ionic_dipole_e_bohr": ionic_dipole.tolist(),
            "electronic_dipole_e_bohr": electronic_dipole.tolist(),
            "neutrality_error_relative": abs(
                cube.electron_count_raw - float(np.sum(cube.charges))
            ) / float(np.sum(cube.charges)),
        },
    )


def compare_slab_charge_profiles(
    reference_cube: str | Path,
    displaced_cube: str | Path,
    *,
    displacement_angstrom: Optional[float] = None,
    neutrality_tolerance: float = 0.05,
    ionic_valence_charges: Optional[Sequence[float]] = None,
) -> SlabChargeDifference:
    """Resolve a slab dipole change into planar electronic and ionic terms."""

    reference = _read_cube_profile_data(
        reference_cube,
        neutrality_tolerance=neutrality_tolerance,
        ionic_valence_charges=ionic_valence_charges,
    )
    displaced = _read_cube_profile_data(
        displaced_cube,
        neutrality_tolerance=neutrality_tolerance,
        ionic_valence_charges=ionic_valence_charges,
    )
    if reference.dimensions != displaced.dimensions:
        raise ValueError("Reference and displaced cube grids have different dimensions")
    if not np.allclose(reference.origin, displaced.origin, atol=1.0e-8):
        raise ValueError("Reference and displaced cube grids have different origins")
    if not np.allclose(reference.step_vectors, displaced.step_vectors, atol=1.0e-8):
        raise ValueError("Reference and displaced cube grids have different steps")
    if reference.charges.shape != displaced.charges.shape or not np.allclose(
        reference.charges, displaced.charges, atol=1.0e-8
    ):
        raise ValueError("Reference and displaced cubes have different ionic charges")

    normal = reference.normal
    ionic_total = float(np.sum(reference.charges))
    atom_normal = reference.positions @ normal
    origin_normal = float(reference.origin @ normal)
    center = _periodic_weighted_center(
        atom_normal,
        reference.charges,
        origin_normal,
        reference.height_bohr,
    )

    grid_normal = (
        origin_normal
        + np.arange(reference.dimensions[2], dtype=float)
        * float(reference.step_vectors[2] @ normal)
    )
    grid_unwrapped = _periodic_unwrap(
        grid_normal, center, reference.height_bohr
    )
    order = np.argsort(grid_unwrapped)
    coordinate_angstrom = (
        grid_unwrapped[order] - center
    ) * BOHR_M / 1.0e-10

    reference_plane_electrons = (
        np.sum(reference.density, axis=(0, 1))
        * reference.voxel_volume_bohr3
        * reference.neutrality_scale
    )
    displaced_plane_electrons = (
        np.sum(displaced.density, axis=(0, 1))
        * displaced.voxel_volume_bohr3
        * displaced.neutrality_scale
    )
    reference_plane_electrons = reference_plane_electrons[order]
    displaced_plane_electrons = displaced_plane_electrons[order]
    dz_angstrom = (
        abs(float(reference.step_vectors[2] @ normal))
        * BOHR_M
        / 1.0e-10
    )
    reference_line_density = reference_plane_electrons / dz_angstrom
    displaced_line_density = displaced_plane_electrons / dz_angstrom
    electron_charge_difference = -(
        displaced_line_density - reference_line_density
    )
    electron_plane_charge = electron_charge_difference * dz_angstrom

    reference_ion_positions = (
        _periodic_unwrap(
            reference.positions @ normal, center, reference.height_bohr
        )
        - center
    ) * BOHR_M / 1.0e-10
    displaced_ion_positions = (
        _periodic_unwrap(
            displaced.positions @ normal, center, reference.height_bohr
        )
        - center
    ) * BOHR_M / 1.0e-10

    cumulative_electron_charge = np.cumsum(electron_plane_charge)
    cumulative_electron_dipole = np.cumsum(
        electron_plane_charge * coordinate_angstrom
    )
    cumulative_ionic_charge = np.zeros_like(coordinate_angstrom)
    cumulative_ionic_dipole = np.zeros_like(coordinate_angstrom)
    thresholds = coordinate_angstrom + 0.5 * dz_angstrom
    thresholds[-1] = np.inf
    for index, threshold in enumerate(thresholds):
        displaced_mask = displaced_ion_positions <= threshold
        reference_mask = reference_ion_positions <= threshold
        cumulative_ionic_charge[index] = float(
            np.sum(reference.charges[displaced_mask])
            - np.sum(reference.charges[reference_mask])
        )
        cumulative_ionic_dipole[index] = float(
            np.sum(
                reference.charges[displaced_mask]
                * displaced_ion_positions[displaced_mask]
            )
            - np.sum(
                reference.charges[reference_mask]
                * reference_ion_positions[reference_mask]
            )
        )

    cumulative_charge = cumulative_electron_charge + cumulative_ionic_charge
    cumulative_dipole = cumulative_electron_dipole + cumulative_ionic_dipole
    direct_reference = integrate_slab_dipole(
        reference.path,
        neutrality_tolerance=neutrality_tolerance,
        ionic_valence_charges=ionic_valence_charges,
    )
    direct_displaced = integrate_slab_dipole(
        displaced.path,
        neutrality_tolerance=neutrality_tolerance,
        ionic_valence_charges=ionic_valence_charges,
    )
    direct_delta_bohr = (
        direct_displaced.dipole_e_bohr - direct_reference.dipole_e_bohr
    )
    direct_delta_bohr -= round(
        direct_delta_bohr / reference.height_bohr
    ) * reference.height_bohr
    direct_delta_e_angstrom = direct_delta_bohr * BOHR_M / 1.0e-10
    closure_error = float(cumulative_dipole[-1] - direct_delta_e_angstrom)

    area_angstrom2 = (
        reference.area_bohr2 * (BOHR_M / 1.0e-10) ** 2
    )
    sheet_polarization = (
        direct_delta_e_angstrom
        / area_angstrom2
        * ELEMENTARY_CHARGE
        / 1.0e-10
    )
    displacement = (
        None if displacement_angstrom is None else float(displacement_angstrom)
    )
    if displacement is not None and (
        not math.isfinite(displacement) or displacement == 0.0
    ):
        raise ValueError("displacement_angstrom must be finite and non-zero")
    effective_charge = (
        None
        if displacement is None
        else direct_delta_e_angstrom / displacement
    )

    return SlabChargeDifference(
        reference_cube=str(reference.path),
        displaced_cube=str(displaced.path),
        coordinate_angstrom=coordinate_angstrom,
        reference_electron_line_density_e_per_angstrom=reference_line_density,
        displaced_electron_line_density_e_per_angstrom=displaced_line_density,
        electron_charge_difference_e_per_angstrom=electron_charge_difference,
        cumulative_charge_difference_e=cumulative_charge,
        cumulative_dipole_difference_e_angstrom=cumulative_dipole,
        reference_ion_positions_angstrom=reference_ion_positions,
        displaced_ion_positions_angstrom=displaced_ion_positions,
        total_dipole_change_e_angstrom=float(direct_delta_e_angstrom),
        sheet_polarization_change_C_per_m=float(sheet_polarization),
        displacement_angstrom=displacement,
        effective_charge_e=(
            None if effective_charge is None else float(effective_charge)
        ),
        diagnostics={
            "reference_electron_count_raw": reference.electron_count_raw,
            "displaced_electron_count_raw": displaced.electron_count_raw,
            "reference_neutrality_scale": reference.neutrality_scale,
            "displaced_neutrality_scale": displaced.neutrality_scale,
            "reference_dipole_e_angstrom": (
                direct_reference.dipole_e_bohr * BOHR_M / 1.0e-10
            ),
            "displaced_dipole_e_angstrom": (
                direct_displaced.dipole_e_bohr * BOHR_M / 1.0e-10
            ),
            "area_angstrom2": area_angstrom2,
            "height_angstrom": reference.height_bohr * BOHR_M / 1.0e-10,
            "profile_closure_error_e_angstrom": closure_error,
        },
    )


def write_slab_charge_difference(
    outdir: str | Path,
    result: SlabChargeDifference,
    *,
    plot: bool = True,
) -> dict:
    """Write source data and a compact visualization of a slab dipole change."""

    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "slab_charge_profile.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "z_angstrom",
                "reference_electron_line_density_e_per_angstrom",
                "displaced_electron_line_density_e_per_angstrom",
                "electron_charge_difference_e_per_angstrom",
                "cumulative_charge_difference_e",
                "cumulative_dipole_difference_e_angstrom",
            ]
        )
        for row in zip(
            result.coordinate_angstrom,
            result.reference_electron_line_density_e_per_angstrom,
            result.displaced_electron_line_density_e_per_angstrom,
            result.electron_charge_difference_e_per_angstrom,
            result.cumulative_charge_difference_e,
            result.cumulative_dipole_difference_e_angstrom,
        ):
            writer.writerow([float(value) for value in row])

    summary_path = output / "slab_dipole_summary.json"
    summary = {
        "reference_cube": Path(result.reference_cube).name,
        "displaced_cube": Path(result.displaced_cube).name,
        "total_dipole_change_e_angstrom": result.total_dipole_change_e_angstrom,
        "sheet_polarization_change_C_per_m": (
            result.sheet_polarization_change_C_per_m
        ),
        "displacement_angstrom": result.displacement_angstrom,
        "effective_charge_e": result.effective_charge_e,
        "reference_ion_positions_angstrom": (
            result.reference_ion_positions_angstrom.tolist()
        ),
        "displaced_ion_positions_angstrom": (
            result.displaced_ion_positions_angstrom.tolist()
        ),
        "diagnostics": result.diagnostics,
        "files": {"profile": csv_path.name},
    }

    if plot:
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        with mpl.rc_context(
            {
                "font.family": "sans-serif",
                "font.sans-serif": ["Arial", "DejaVu Sans"],
                "font.size": 8,
                "axes.spines.right": False,
                "axes.spines.top": False,
                "pdf.fonttype": 42,
                "svg.fonttype": "none",
            }
        ):
            fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.8))
            axes[0].plot(
                result.coordinate_angstrom,
                result.electron_charge_difference_e_per_angstrom,
                color="#2f6b9a",
                linewidth=1.4,
            )
            axes[0].axhline(0.0, color="#6b7280", linewidth=0.6)
            for position in result.reference_ion_positions_angstrom:
                axes[0].axvline(position, color="#d1d5db", linewidth=0.5)
            axes[0].set_xlabel(r"Slab-normal coordinate ($\AA$)")
            axes[0].set_ylabel(r"$\Delta\lambda_e(z)$ ($e$ $\AA^{-1}$)")
            if result.displacement_angstrom is None:
                finite_difference_x = np.asarray([0.0, 1.0])
                axes[1].set_xticks(
                    finite_difference_x, ["reference", "displaced"]
                )
                axes[1].set_xlabel("Charge-density cube")
            else:
                finite_difference_x = np.asarray(
                    [0.0, result.displacement_angstrom]
                )
                axes[1].set_xlabel(r"Ionic displacement $\Delta u_z$ ($\AA$)")
            axes[1].plot(
                finite_difference_x,
                [0.0, result.total_dipole_change_e_angstrom],
                color="#b2472f",
                marker="o",
                linewidth=1.6,
            )
            axes[1].axhline(0.0, color="#6b7280", linewidth=0.6)
            axes[1].set_ylabel(
                r"$\mu_z-\mu_z^{\mathrm{ref}}$ ($e\AA$)"
            )
            annotation = (
                rf"$\Delta\mu_z={result.total_dipole_change_e_angstrom:.4g}"
                r"\ e\AA$"
            )
            if result.effective_charge_e is not None:
                annotation += (
                    "\n"
                    rf"$Z^*_{{zz}}={result.effective_charge_e:.4g}\ e$"
                )
            axes[1].text(
                0.04,
                0.92,
                annotation,
                transform=axes[1].transAxes,
                ha="left",
                va="top",
            )
            axes[1].grid(axis="y", color="#d9dde1", linewidth=0.5)
            fig.tight_layout()
            plot_base = output / "slab_dipole_profile"
            fig.savefig(plot_base.with_suffix(".png"), dpi=300)
            fig.savefig(plot_base.with_suffix(".pdf"))
            fig.savefig(plot_base.with_suffix(".svg"))
            plt.close(fig)
        summary["files"].update(
            {
                "plot": plot_base.with_suffix(".png").name,
                "plot_pdf": plot_base.with_suffix(".pdf").name,
                "plot_svg": plot_base.with_suffix(".svg").name,
            }
        )

    summary["files"]["summary"] = summary_path.name
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


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


def read_pyatb_lattice(input_json: str | Path) -> np.ndarray:
    input_json = Path(input_json)
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
    return matrix_angstrom


def _read_pyatb_cell(input_json: Path) -> tuple[np.ndarray, float]:
    matrix_angstrom = read_pyatb_lattice(input_json)
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
        tensor[beta, :2] = in_plane
        tensor[beta, 2] = z_born
        diagnostics["directions"][direction] = {
            "plus": plus_dipole.to_dict(),
            "minus": minus_report,
            "berry_delta_C_per_m2": delta_lattice.tolist(),
            "in_plane_column": in_plane.tolist(),
            "out_of_plane_component": float(z_born),
        }

    return Hybrid2DBorn(
        tensor=tensor,
        in_plane_tensor=tensor[:, :2].copy(),
        out_of_plane_column=tensor[:, 2].copy(),
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
        "tensor_convention": (
            "rows=atomic displacement/force; "
            "columns=polarization/electric field"
        ),
        "tensor": result.tensor.tolist(),
        "in_plane_tensor": result.in_plane_tensor.tolist(),
        "out_of_plane_column": result.out_of_plane_column.tolist(),
        "diagnostics": result.diagnostics,
    }
    Path(path).write_text(json.dumps(data, indent=2), encoding="utf-8")
