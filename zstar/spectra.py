"""Infrared and Raman spectra for three-dimensional crystals and 2D slabs."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import csv
import json
import math
from pathlib import Path
import re
import shutil
from typing import Iterable, Optional, Sequence

import numpy as np
import yaml

from .pyatb_compat import read_static_dielectric
from .stru_analyzer import stru_analyzer


THZ_TO_CM1 = 33.35640951981521
ELEMENTARY_CHARGE = 1.602176634e-19
VACUUM_PERMITTIVITY = 8.8541878128e-12
ATOMIC_MASS_UNIT = 1.66053906660e-27
SPEED_OF_LIGHT = 299792458.0
BOLTZMANN = 1.380649e-23
PLANCK = 6.62607015e-34
BOHR_ANGSTROM = 0.529177210903
DEBYE_C_M = 3.33564e-30
_FLOAT_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?")
_PYATB_POLARIZATION_RE = re.compile(
    rf"direction is in\s+([abc]),\s*P\s*=\s*({_FLOAT_RE.pattern})\s*"
    rf"\(mod\s*({_FLOAT_RE.pattern})\)\s*C/m\^2",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class GammaModes:
    frequencies_thz: np.ndarray
    eigenvectors: np.ndarray
    masses_amu: np.ndarray
    lattice_angstrom: np.ndarray
    symbols: tuple[str, ...]
    positions_fractional: np.ndarray

    @property
    def frequencies_cm1(self) -> np.ndarray:
        return self.frequencies_thz * THZ_TO_CM1

    @property
    def volume_angstrom3(self) -> float:
        return float(abs(np.linalg.det(self.lattice_angstrom)))

    @property
    def area_angstrom2(self) -> float:
        return float(
            np.linalg.norm(
                np.cross(self.lattice_angstrom[0], self.lattice_angstrom[1])
            )
        )

    @property
    def cell_height_angstrom(self) -> float:
        return self.volume_angstrom3 / self.area_angstrom2


@dataclass(frozen=True)
class BornData:
    tensors: np.ndarray
    electronic_dielectric: Optional[np.ndarray]
    source: str


@dataclass(frozen=True)
class IRSpectrumResult:
    mode_numbers: np.ndarray
    frequencies_cm1: np.ndarray
    effective_charges: np.ndarray
    intensities: np.ndarray
    frequency_grid_cm1: np.ndarray
    spectrum: np.ndarray
    response_real: np.ndarray
    response_imag: np.ndarray
    dimensionality: int
    response_kind: str


@dataclass(frozen=True)
class MolecularIRSpectrumResult:
    mode_numbers: np.ndarray
    frequencies_cm1: np.ndarray
    dipole_derivatives: np.ndarray
    activities: np.ndarray
    normalized_activities: np.ndarray
    frequency_grid_cm1: np.ndarray
    spectrum: np.ndarray


@dataclass(frozen=True)
class RamanSpectrumResult:
    mode_numbers: np.ndarray
    frequencies_cm1: np.ndarray
    tensors: np.ndarray
    activities: np.ndarray
    depolarization_ratios: np.ndarray
    frequency_grid_cm1: np.ndarray
    spectrum: np.ndarray
    tensor_kind: str


@dataclass(frozen=True)
class NativeLineSpectrumResult:
    mode_numbers: np.ndarray
    frequencies_cm1: np.ndarray
    activities: np.ndarray
    frequency_grid_cm1: np.ndarray
    spectrum: np.ndarray
    activity_kind: str
    activity_unit: str


def _complex_component(value) -> complex:
    if isinstance(value, (list, tuple)):
        if len(value) == 2 and all(isinstance(item, (int, float)) for item in value):
            return complex(float(value[0]), float(value[1]))
        if len(value) == 1:
            return complex(float(value[0]), 0.0)
    return complex(float(value), 0.0)


def load_gamma_modes(path: str | Path = "qpoints.yaml") -> GammaModes:
    """Load the first Gamma-point eigensystem from a Phonopy YAML file."""

    yaml_path = Path(path)
    data = yaml.safe_load(yaml_path.read_text(encoding="utf-8"))
    phonons = data.get("phonon") or []
    if not phonons:
        raise ValueError(f"No phonon entries found in {yaml_path}")

    gamma = None
    for entry in phonons:
        q_position = np.asarray(entry.get("q-position", [0.0, 0.0, 0.0]), dtype=float)
        if np.linalg.norm(q_position) < 1.0e-7:
            gamma = entry
            break
    gamma = gamma or phonons[0]
    bands = gamma.get("band") or []
    if not bands:
        raise ValueError(f"No phonon bands found in {yaml_path}")

    companion = yaml_path.with_name("phonopy.yaml")
    metadata = None
    if companion.is_file() and companion.resolve() != yaml_path.resolve():
        metadata = yaml.safe_load(companion.read_text(encoding="utf-8"))

    cell = data.get("primitive_cell") or data.get("unit_cell") or {}
    if not cell.get("points"):
        if metadata:
            cell = (
                metadata.get("primitive_cell")
                or metadata.get("unit_cell")
                or {}
            )
    points = cell.get("points") or []
    if not points:
        raise ValueError(
            f"No primitive-cell points found in {yaml_path}; place the matching "
            "phonopy.yaml in the same directory."
        )
    masses = np.asarray([point["mass"] for point in points], dtype=float)
    symbols = tuple(str(point.get("symbol", "X")) for point in points)
    positions = np.asarray(
        [point.get("coordinates", [0.0, 0.0, 0.0]) for point in points],
        dtype=float,
    )
    lattice = np.asarray(cell.get("lattice"), dtype=float)
    if lattice.shape != (3, 3):
        raise ValueError(f"Invalid primitive-cell lattice in {yaml_path}")
    physical_unit = data.get("physical_unit") or {}
    if not physical_unit and metadata:
        physical_unit = metadata.get("physical_unit") or {}
    length_unit = str(physical_unit.get("length", "angstrom")).lower()
    if length_unit in {"au", "bohr", "a.u."}:
        lattice *= BOHR_ANGSTROM

    frequencies: list[float] = []
    eigenvectors: list[np.ndarray] = []
    for band in bands:
        frequencies.append(float(band["frequency"]))
        raw = band.get("eigenvector")
        if raw is None:
            raise ValueError("qpoints.yaml was written without eigenvectors")
        vector = np.empty((len(points), 3), dtype=complex)
        for atom_index, atom_vector in enumerate(raw):
            for direction in range(3):
                vector[atom_index, direction] = _complex_component(
                    atom_vector[direction]
                )
        eigenvectors.append(vector)

    result = GammaModes(
        frequencies_thz=np.asarray(frequencies, dtype=float),
        eigenvectors=np.asarray(eigenvectors, dtype=complex),
        masses_amu=masses,
        lattice_angstrom=lattice,
        symbols=symbols,
        positions_fractional=positions,
    )
    if result.eigenvectors.shape[1] != len(result.masses_amu):
        raise ValueError("Eigenvector atom count does not match primitive-cell masses")
    return result


def _numeric_lines(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        values = [float(item) for item in _FLOAT_RE.findall(line)]
        if len(values) >= 9:
            rows.append(values)
    return rows


def read_born_data(
    path: str | Path,
    *,
    natoms: Optional[int] = None,
    dielectric_path: Optional[str | Path] = None,
) -> BornData:
    """Read labeled Z-BORN output or a Phonopy-style BORN file."""

    born_path = Path(path)
    rows = _numeric_lines(born_path)
    if not rows:
        raise ValueError(f"No 3x3 tensors found in {born_path}")

    electronic: Optional[np.ndarray] = None
    tensors: list[np.ndarray] = []
    is_phonopy_born = born_path.name.upper() == "BORN" or all(
        len(row) == 9 for row in rows
    )

    if is_phonopy_born and len(rows) >= 2:
        electronic = np.asarray(rows[0][-9:], dtype=float).reshape(3, 3)
        tensors = [
            np.asarray(row[-9:], dtype=float).reshape(3, 3) for row in rows[1:]
        ]
    else:
        tensors = [
            np.asarray(row[-9:], dtype=float).reshape(3, 3) for row in rows
        ]

    if dielectric_path:
        dielectric_file = Path(dielectric_path)
        if dielectric_file.name.upper() == "BORN":
            dielectric_rows = _numeric_lines(dielectric_file)
            if dielectric_rows:
                electronic = np.asarray(
                    dielectric_rows[0][-9:], dtype=float
                ).reshape(3, 3)
        else:
            electronic, _ = read_static_dielectric(dielectric_file)

    tensor_source = born_path.resolve()
    if natoms is not None and len(tensors) != natoms:
        for candidate_name in ("Z-BORN-all.out", "Z-BORN-symm.out"):
            candidate = born_path.parent / candidate_name
            if not candidate.is_file():
                continue
            expanded_rows = _numeric_lines(candidate)
            if len(expanded_rows) == natoms:
                tensors = [
                    np.asarray(row[-9:], dtype=float).reshape(3, 3)
                    for row in expanded_rows
                ]
                tensor_source = candidate.resolve()
                break

    if natoms is not None and len(tensors) != natoms:
        raise ValueError(
            f"Born tensor count ({len(tensors)}) does not match phonon atoms "
            f"({natoms}). Provide a full per-atom Z-BORN-all.out or "
            "Z-BORN-symm.out in the same directory."
        )
    return BornData(
        tensors=np.asarray(tensors, dtype=float),
        electronic_dielectric=electronic,
        source=(
            str(born_path.resolve())
            if tensor_source == born_path.resolve()
            else f"{born_path.resolve()} + {tensor_source}"
        ),
    )


def _mode_phase_real(eigenvectors: np.ndarray) -> np.ndarray:
    """Choose a stable Gamma-mode phase and return its real representation."""

    output = np.empty(eigenvectors.shape, dtype=float)
    for mode_index, vector in enumerate(eigenvectors):
        flat = vector.reshape(-1)
        pivot = int(np.argmax(np.abs(flat)))
        phase = np.angle(flat[pivot])
        phased = vector * np.exp(-1j * phase)
        imaginary_norm = float(np.linalg.norm(phased.imag))
        real_norm = float(np.linalg.norm(phased.real))
        if imaginary_norm > max(1.0e-7, 1.0e-5 * real_norm):
            raise ValueError(
                f"Gamma eigenvector {mode_index + 1} cannot be represented as real "
                f"(imaginary norm {imaginary_norm:.3e})"
            )
        output[mode_index] = phased.real
    return output


def mode_effective_charges(
    modes: GammaModes,
    born_tensors: np.ndarray,
) -> np.ndarray:
    """Calculate mass-weighted mode effective-charge vectors."""

    born = np.asarray(born_tensors, dtype=float)
    if born.shape != (len(modes.masses_amu), 3, 3):
        raise ValueError(
            f"Born tensors must have shape {(len(modes.masses_amu), 3, 3)}, "
            f"got {born.shape}"
        )
    eigenvectors = _mode_phase_real(modes.eigenvectors)
    mass_weighted = eigenvectors / np.sqrt(modes.masses_amu)[None, :, None]
    return np.einsum("aij,maj->mi", born, mass_weighted)


def _selected_mode_indices(
    frequencies_cm1: np.ndarray,
    mode_numbers: Optional[Sequence[int]],
    acoustic_cutoff_cm1: float,
) -> np.ndarray:
    if mode_numbers:
        indices = np.asarray([int(number) - 1 for number in mode_numbers], dtype=int)
        if np.any(indices < 0) or np.any(indices >= len(frequencies_cm1)):
            raise IndexError("Requested mode number is outside qpoints.yaml")
        return indices
    return np.flatnonzero(frequencies_cm1 > float(acoustic_cutoff_cm1))


def _lorentzian(grid: np.ndarray, centers: np.ndarray, gamma: float) -> np.ndarray:
    half_width = max(float(gamma), 1.0e-12) / 2.0
    delta = grid[:, None] - centers[None, :]
    return (half_width / math.pi) / (delta * delta + half_width * half_width)


def calculate_native_line_spectrum(
    frequencies_cm1: Sequence[float],
    activities: Sequence[float],
    *,
    mode_numbers: Optional[Sequence[int]] = None,
    activity_kind: str,
    activity_unit: str,
    broadening_cm1: float = 8.0,
    max_frequency_cm1: Optional[float] = None,
    points: int = 2001,
) -> NativeLineSpectrumResult:
    """Broaden calculator-native IR or Raman activities without re-scaling them."""

    frequencies = np.asarray(frequencies_cm1, dtype=float)
    values = np.asarray(activities, dtype=float)
    if frequencies.ndim != 1 or values.shape != frequencies.shape:
        raise ValueError("frequencies and activities must be one-dimensional and aligned")
    numbers = (
        np.asarray(mode_numbers, dtype=int)
        if mode_numbers is not None
        else np.arange(1, len(frequencies) + 1, dtype=int)
    )
    if numbers.shape != frequencies.shape:
        raise ValueError("mode_numbers must align with frequencies")
    keep = frequencies > 0.0
    frequencies = frequencies[keep]
    values = np.maximum(values[keep], 0.0)
    numbers = numbers[keep]
    if not len(frequencies):
        raise ValueError("No positive-frequency modes are available")
    upper = max_frequency_cm1
    if upper is None:
        upper = max(100.0, float(np.max(frequencies)) + 5.0 * broadening_cm1)
    grid = np.linspace(0.0, float(upper), int(points))
    spectrum = _lorentzian(grid, frequencies, broadening_cm1) @ values
    return NativeLineSpectrumResult(
        mode_numbers=numbers,
        frequencies_cm1=frequencies,
        activities=values,
        frequency_grid_cm1=grid,
        spectrum=spectrum,
        activity_kind=activity_kind,
        activity_unit=activity_unit,
    )


def write_native_line_spectrum_outputs(
    outdir: str | Path,
    result: NativeLineSpectrumResult,
    *,
    stem: str,
    plot: bool = True,
) -> dict:
    """Write calculator-native line activities and a normalized display spectrum."""

    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    modes_path = output / f"{stem}_modes.csv"
    with modes_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["mode", "frequency_cm-1", result.activity_kind, "unit"])
        for number, frequency, activity in zip(
            result.mode_numbers, result.frequencies_cm1, result.activities
        ):
            writer.writerow(
                [int(number), float(frequency), float(activity), result.activity_unit]
            )
    spectrum_path = output / f"{stem}_spectrum.dat"
    np.savetxt(
        spectrum_path,
        np.column_stack([result.frequency_grid_cm1, result.spectrum]),
        header=f"frequency_cm-1 {result.activity_kind}",
    )

    plot_files: dict[str, str] = {}
    if plot:
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        with mpl.rc_context(
            {
                "font.family": "sans-serif",
                "font.sans-serif": ["Arial", "DejaVu Sans"],
                "font.size": 9,
                "axes.spines.right": True,
                "axes.spines.top": True,
                "axes.linewidth": 0.8,
                "xtick.direction": "in",
                "ytick.direction": "in",
                "xtick.top": True,
                "ytick.right": True,
                "pdf.fonttype": 42,
                "svg.fonttype": "none",
            }
        ):
            fig, ax = plt.subplots(figsize=(7.2, 4.5))
            scale = max(float(np.max(result.spectrum)), np.finfo(float).tiny)
            normalized = result.spectrum / scale
            color = "#b14b3c" if stem.lower().startswith("raman") else "#2f6b9a"
            ax.plot(result.frequency_grid_cm1, normalized, color=color, linewidth=1.5)
            activity_scale = max(float(np.max(result.activities)), np.finfo(float).tiny)
            ax.vlines(
                result.frequencies_cm1,
                0.0,
                0.12 * result.activities / activity_scale,
                color="#30343b",
                linewidth=0.9,
            )
            ax.set(
                xlabel=r"Wavenumber (cm$^{-1}$)",
                ylabel=f"Normalized {stem.upper()} intensity",
                xlim=(result.frequency_grid_cm1[0], result.frequency_grid_cm1[-1]),
                ylim=(0.0, 1.05),
            )
            fig.tight_layout()
            plot_files = _save_figure_bundle(fig, output, f"{stem}_spectrum")
            plt.close(fig)

    summary = {
        "activity_kind": result.activity_kind,
        "activity_unit": result.activity_unit,
        "modes": len(result.mode_numbers),
        "normalization": "plot maximum equals one; tabulated activities are unchanged",
        "files": {
            "modes": str(modes_path.resolve()),
            "spectrum": str(spectrum_path.resolve()),
            "plot": plot_files.get("plot"),
            "plot_pdf": plot_files.get("plot_pdf"),
            "plot_svg": plot_files.get("plot_svg"),
        },
    }
    (output / f"{stem}_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    return summary


def calculate_ir_spectrum(
    modes: GammaModes,
    born: BornData,
    *,
    dimensionality: int = 3,
    mode_numbers: Optional[Sequence[int]] = None,
    acoustic_cutoff_cm1: float = 5.0,
    broadening_cm1: float = 10.0,
    max_frequency_cm1: Optional[float] = None,
    points: int = 2001,
    thickness_angstrom: Optional[float] = None,
) -> IRSpectrumResult:
    """Calculate IR oscillator strengths and dielectric/sheet response."""

    if dimensionality not in (2, 3):
        raise ValueError("dimensionality must be 2 or 3")
    effective_all = mode_effective_charges(modes, born.tensors)
    frequencies_all = modes.frequencies_cm1
    indices = _selected_mode_indices(
        frequencies_all, mode_numbers, acoustic_cutoff_cm1
    )
    frequencies = frequencies_all[indices]
    effective = effective_all[indices]
    intensities = effective * effective

    upper = max_frequency_cm1
    if upper is None:
        upper = max(100.0, float(np.max(frequencies)) + 5.0 * broadening_cm1)
    grid = np.linspace(0.0, float(upper), int(points))
    line_shapes = _lorentzian(grid, frequencies, broadening_cm1)
    spectrum = line_shapes @ intensities

    oscillator_tensors = np.einsum("mi,mj->mij", effective, effective)
    denominator = (
        frequencies[None, :] ** 2
        - grid[:, None] ** 2
        - 1j * float(broadening_cm1) * grid[:, None]
    )
    common = (
        ELEMENTARY_CHARGE**2
        / (
            VACUUM_PERMITTIVITY
            * ATOMIC_MASS_UNIT
            * (2.0 * math.pi * SPEED_OF_LIGHT * 100.0) ** 2
        )
    )
    if dimensionality == 3:
        prefactor = common / (modes.volume_angstrom3 * 1.0e-30)
        response = prefactor * np.einsum(
            "mij,wm->wij", oscillator_tensors, 1.0 / denominator
        )
        if born.electronic_dielectric is None:
            response += np.eye(3)[None, :, :]
        else:
            response += born.electronic_dielectric[None, :, :]
        response_kind = "relative dielectric tensor"
    else:
        prefactor = common / (modes.area_angstrom2 * 1.0e-20)
        response = prefactor * np.einsum(
            "mij,wm->wij", oscillator_tensors, 1.0 / denominator
        )
        if born.electronic_dielectric is not None:
            electronic_sheet_m = (
                modes.cell_height_angstrom
                * 1.0e-10
                * (born.electronic_dielectric - np.eye(3))
            )
            response += electronic_sheet_m[None, :, :]
        if thickness_angstrom is not None:
            response = (
                np.eye(3)[None, :, :]
                + response / (float(thickness_angstrom) * 1.0e-10)
            )
            response_kind = "effective relative dielectric tensor"
        else:
            response *= 1.0e10
            response_kind = "2D sheet polarizability (Angstrom)"

    return IRSpectrumResult(
        mode_numbers=indices + 1,
        frequencies_cm1=frequencies,
        effective_charges=effective,
        intensities=intensities,
        frequency_grid_cm1=grid,
        spectrum=spectrum,
        response_real=response.real,
        response_imag=response.imag,
        dimensionality=dimensionality,
        response_kind=response_kind,
    )


def read_pyatb_polarization(
    path_or_directory: str | Path,
) -> tuple[np.ndarray, np.ndarray, Path]:
    """Read Berry polarization and branch quanta from a PYATB output."""

    root = Path(path_or_directory)
    candidates = [root] if root.is_file() else [
        root / "Out" / "Polarization" / "polarization.dat",
        root / "Polarization" / "polarization.dat",
        root / "polarization.dat",
    ]
    for candidate in candidates:
        if not candidate.is_file() or candidate.stat().st_size == 0:
            continue
        values: dict[str, float] = {}
        quanta: dict[str, float] = {}
        for line in candidate.read_text(
            encoding="utf-8", errors="ignore"
        ).splitlines():
            match = _PYATB_POLARIZATION_RE.search(line)
            if match:
                axis = match.group(1).lower()
                values[axis] = float(match.group(2))
                quanta[axis] = float(match.group(3))
        if set(values) == {"a", "b", "c"}:
            return (
                np.asarray([values[axis] for axis in "abc"], dtype=float),
                np.asarray([quanta[axis] for axis in "abc"], dtype=float),
                candidate.resolve(),
            )
    raise FileNotFoundError(
        f"No complete PYATB Polarization/polarization.dat found under {root}"
    )


def _manifest_stage_path(root: Path, entry: dict, sign: str) -> Path:
    recorded = Path(entry[sign])
    if recorded.is_dir():
        return recorded.resolve()
    fallback = root / f"mode-{int(entry['mode']):04d}" / sign
    if fallback.is_dir():
        return fallback.resolve()
    raise FileNotFoundError(
        f"Molecular mode {entry['mode']} {sign} directory was not found"
    )


def _stage_polarization(
    stage: Path,
    preferred_subdir: Optional[str],
) -> tuple[np.ndarray, np.ndarray, Path]:
    names = [preferred_subdir] if preferred_subdir else []
    names.extend(["pyatb", "pyatb-polar"])
    errors: list[str] = []
    for name in dict.fromkeys(item for item in names if item):
        try:
            return read_pyatb_polarization(stage / name)
        except FileNotFoundError as exc:
            errors.append(str(exc))
    try:
        return read_pyatb_polarization(stage)
    except FileNotFoundError as exc:
        errors.append(str(exc))
    raise FileNotFoundError("; ".join(errors))


def collect_molecular_dipole_derivatives(
    displacement_dir: str | Path,
    *,
    cell_volume_angstrom3: float,
    cell_lattice_angstrom: Optional[np.ndarray] = None,
    polarization_subdir: Optional[str] = None,
) -> tuple[np.ndarray, np.ndarray, str]:
    """Collect central-difference molecular dipole derivatives from PYATB P."""

    root = Path(displacement_dir).resolve()
    manifest_path = root / "raman_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    amplitude = float(manifest["amplitude_A_sqrt_amu"])
    if amplitude <= 0.0:
        raise ValueError("Molecular displacement amplitude must be positive")
    volume_m3 = float(cell_volume_angstrom3) * 1.0e-30
    if volume_m3 <= 0.0:
        raise ValueError("Molecular supercell volume must be positive")
    if cell_lattice_angstrom is None:
        basis_to_cartesian = np.eye(3)
    else:
        lattice = np.asarray(cell_lattice_angstrom, dtype=float)
        if lattice.shape != (3, 3):
            raise ValueError("Molecular cell lattice must have shape (3, 3)")
        lengths = np.linalg.norm(lattice, axis=1)
        if np.any(lengths <= np.finfo(float).tiny):
            raise ValueError("Molecular cell lattice vectors must be nonzero")
        basis_to_cartesian = lattice / lengths[:, None]

    mode_numbers: list[int] = []
    derivatives: list[np.ndarray] = []
    for entry in manifest["modes"]:
        plus_stage = _manifest_stage_path(root, entry, "plus")
        minus_stage = _manifest_stage_path(root, entry, "minus")
        plus, plus_quanta, plus_source = _stage_polarization(
            plus_stage, polarization_subdir
        )
        minus, minus_quanta, minus_source = _stage_polarization(
            minus_stage, polarization_subdir
        )
        delta = plus - minus
        quanta = 0.5 * (np.abs(plus_quanta) + np.abs(minus_quanta))
        valid = quanta > np.finfo(float).tiny
        delta[valid] -= np.rint(delta[valid] / quanta[valid]) * quanta[valid]
        delta_cartesian = delta @ basis_to_cartesian
        derivative = (
            delta_cartesian * volume_m3 / (2.0 * amplitude * DEBYE_C_M)
        )
        mode_numbers.append(int(entry["mode"]))
        derivatives.append(derivative)
        entry["plus_polarization"] = str(plus_source)
        entry["minus_polarization"] = str(minus_source)

    derivative_array = np.asarray(derivatives, dtype=float)
    tensor_kind = "molecular dipole derivative (Debye per Angstrom sqrt(amu))"
    output = {
        "mode_numbers": mode_numbers,
        "dipole_derivatives": derivative_array.tolist(),
        "tensor_kind": tensor_kind,
        "cell_volume_A3": float(cell_volume_angstrom3),
        "amplitude_A_sqrt_amu": amplitude,
        "conversion": "dmu/dQ = V * dP/dQ",
        "polarization_basis_to_cartesian": basis_to_cartesian.tolist(),
    }
    (root / "molecular_ir_derivatives.json").write_text(
        json.dumps(output, indent=2), encoding="utf-8"
    )
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return np.asarray(mode_numbers, dtype=int), derivative_array, tensor_kind


def calculate_molecular_ir_spectrum(
    modes: GammaModes,
    mode_numbers: Sequence[int],
    dipole_derivatives: np.ndarray,
    *,
    broadening_cm1: float = 10.0,
    max_frequency_cm1: Optional[float] = None,
    points: int = 2001,
) -> MolecularIRSpectrumResult:
    """Calculate a normalized molecular IR spectrum from d(mu)/dQ."""

    numbers = np.asarray(mode_numbers, dtype=int)
    derivatives = np.asarray(dipole_derivatives, dtype=float)
    if derivatives.shape != (len(numbers), 3):
        raise ValueError(
            f"Molecular dipole derivatives have invalid shape: {derivatives.shape}"
        )
    indices = numbers - 1
    if np.any(indices < 0) or np.any(indices >= len(modes.frequencies_cm1)):
        raise IndexError("Molecular IR mode number is outside qpoints.yaml")
    frequencies = modes.frequencies_cm1[indices]
    activities = np.einsum("mi,mi->m", derivatives, derivatives)
    maximum = float(np.max(activities)) if len(activities) else 0.0
    normalized = activities / maximum if maximum > 0.0 else activities.copy()
    upper = max_frequency_cm1
    if upper is None:
        upper = max(100.0, float(np.max(frequencies)) + 5.0 * broadening_cm1)
    grid = np.linspace(0.0, float(upper), int(points))
    spectrum = _lorentzian(grid, frequencies, broadening_cm1) @ activities
    spectrum_maximum = float(np.max(spectrum)) if len(spectrum) else 0.0
    if spectrum_maximum > 0.0:
        spectrum /= spectrum_maximum
    return MolecularIRSpectrumResult(
        mode_numbers=numbers,
        frequencies_cm1=frequencies,
        dipole_derivatives=derivatives,
        activities=activities,
        normalized_activities=normalized,
        frequency_grid_cm1=grid,
        spectrum=spectrum,
    )


def _write_tensor_frequency_file(
    path: Path, frequency: np.ndarray, values: np.ndarray
) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# frequency_cm-1 xx xy xz yx yy yz zx zy zz\n")
        for omega, tensor in zip(frequency, values):
            flat = " ".join(f"{value:.10e}" for value in tensor.reshape(9))
            handle.write(f"{omega:.8f} {flat}\n")


def _save_figure_bundle(fig, output: Path, stem: str) -> dict[str, str]:
    base = output / stem
    paths = {
        "plot": base.with_suffix(".png"),
        "plot_pdf": base.with_suffix(".pdf"),
        "plot_svg": base.with_suffix(".svg"),
    }
    fig.savefig(paths["plot"], dpi=300, bbox_inches="tight")
    fig.savefig(paths["plot_pdf"], bbox_inches="tight")
    fig.savefig(paths["plot_svg"], bbox_inches="tight")
    return {key: str(path.resolve()) for key, path in paths.items()}


def write_ir_outputs(
    outdir: str | Path,
    result: IRSpectrumResult,
    *,
    plot: bool = True,
) -> dict:
    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)

    modes_path = output / "ir_modes.csv"
    with modes_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "mode",
                "frequency_cm-1",
                "Zmode_x",
                "Zmode_y",
                "Zmode_z",
                "intensity_x",
                "intensity_y",
                "intensity_z",
                "intensity_total",
            ]
        )
        for number, frequency, charge, intensity in zip(
            result.mode_numbers,
            result.frequencies_cm1,
            result.effective_charges,
            result.intensities,
        ):
            writer.writerow(
                [
                    int(number),
                    float(frequency),
                    *charge.tolist(),
                    *intensity.tolist(),
                    float(np.sum(intensity)),
                ]
            )

    spectrum_path = output / "ir_spectrum.dat"
    spectrum_values = np.column_stack(
        [result.frequency_grid_cm1, result.spectrum, np.sum(result.spectrum, axis=1)]
    )
    np.savetxt(
        spectrum_path,
        spectrum_values,
        header="frequency_cm-1 intensity_x intensity_y intensity_z intensity_total",
    )
    real_path = output / "ir_response_real.dat"
    imag_path = output / "ir_response_imag.dat"
    _write_tensor_frequency_file(
        real_path, result.frequency_grid_cm1, result.response_real
    )
    _write_tensor_frequency_file(
        imag_path, result.frequency_grid_cm1, result.response_imag
    )

    plot_files: dict[str, str] = {}
    if plot:
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        with mpl.rc_context(
            {
                "font.family": "sans-serif",
                "font.sans-serif": ["Arial", "DejaVu Sans"],
                "font.size": 9,
                "axes.spines.right": True,
                "axes.spines.top": True,
                "axes.linewidth": 0.8,
                "xtick.direction": "in",
                "ytick.direction": "in",
                "xtick.top": True,
                "ytick.right": True,
                "legend.frameon": False,
                "pdf.fonttype": 42,
                "svg.fonttype": "none",
            }
        ):
            fig, ax = plt.subplots(figsize=(7.2, 4.5))
            total = np.sum(result.spectrum, axis=1)
            scale = max(float(np.max(total)), np.finfo(float).tiny)
            colors = ("#2f6b9a", "#5f9c76", "#c06b32")
            for index, (label, color) in enumerate(zip(("x", "y", "z"), colors)):
                ax.plot(
                    result.frequency_grid_cm1,
                    result.spectrum[:, index] / scale,
                    color=color,
                    linewidth=1.1,
                    label=label,
                )
            ax.plot(
                result.frequency_grid_cm1,
                total / scale,
                color="#202124",
                linewidth=1.7,
                label="total",
            )
            mode_strength = np.sum(result.intensities, axis=1)
            if float(np.max(mode_strength)) > 0.0:
                mode_strength = mode_strength / float(np.max(mode_strength))
                ax.vlines(
                    result.frequencies_cm1,
                    -0.075,
                    -0.075 + 0.06 * mode_strength,
                    color="#4b5563",
                    linewidth=0.8,
                )
            ax.set_xlabel(r"Wavenumber (cm$^{-1}$)")
            ax.set_ylabel("Normalized IR intensity")
            ax.set_xlim(
                result.frequency_grid_cm1[0], result.frequency_grid_cm1[-1]
            )
            ax.set_ylim(-0.085, 1.05)
            ax.grid(axis="y", alpha=0.18, linewidth=0.6)
            ax.legend(ncol=4, loc="upper right")
            fig.tight_layout()
            plot_files = _save_figure_bundle(fig, output, "ir_spectrum")
            plt.close(fig)

    summary = {
        "dimensionality": result.dimensionality,
        "response_kind": result.response_kind,
        "modes": len(result.mode_numbers),
        "files": {
            "modes": str(modes_path.resolve()),
            "spectrum": str(spectrum_path.resolve()),
            "response_real": str(real_path.resolve()),
            "response_imag": str(imag_path.resolve()),
            "plot": plot_files.get("plot"),
            "plot_pdf": plot_files.get("plot_pdf"),
            "plot_svg": plot_files.get("plot_svg"),
        },
    }
    (output / "ir_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    return summary


def write_molecular_ir_outputs(
    outdir: str | Path,
    result: MolecularIRSpectrumResult,
    *,
    plot: bool = True,
) -> dict:
    """Write molecular dipole derivatives, normalized spectrum, and plots."""

    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    modes_path = output / "ir_modes.csv"
    with modes_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "mode",
                "frequency_cm-1",
                "dmu_x",
                "dmu_y",
                "dmu_z",
                "activity_Debye2_per_A2_amu",
                "activity_normalized",
            ]
        )
        for number, frequency, derivative, activity, normalized in zip(
            result.mode_numbers,
            result.frequencies_cm1,
            result.dipole_derivatives,
            result.activities,
            result.normalized_activities,
        ):
            writer.writerow(
                [
                    int(number),
                    float(frequency),
                    *derivative.tolist(),
                    float(activity),
                    float(normalized),
                ]
            )

    spectrum_path = output / "ir_spectrum.dat"
    np.savetxt(
        spectrum_path,
        np.column_stack([result.frequency_grid_cm1, result.spectrum]),
        header="frequency_cm-1 intensity_normalized",
    )
    plot_files: dict[str, str] = {}
    if plot:
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        with mpl.rc_context(
            {
                "font.family": "sans-serif",
                "font.sans-serif": ["Arial", "DejaVu Sans"],
                "font.size": 9,
                "axes.spines.right": True,
                "axes.spines.top": True,
                "axes.linewidth": 0.8,
                "xtick.direction": "in",
                "ytick.direction": "in",
                "xtick.top": True,
                "ytick.right": True,
                "pdf.fonttype": 42,
                "svg.fonttype": "none",
            }
        ):
            fig, ax = plt.subplots(figsize=(7.2, 4.5))
            ax.plot(
                result.frequency_grid_cm1,
                result.spectrum,
                color="#b2472f",
                linewidth=1.7,
            )
            ax.vlines(
                result.frequencies_cm1,
                -0.075,
                -0.075 + 0.06 * result.normalized_activities,
                color="#4b5563",
                linewidth=0.8,
            )
            ax.set(
                xlabel=r"Wavenumber (cm$^{-1}$)",
                ylabel="Normalized molecular IR intensity",
                xlim=(
                    result.frequency_grid_cm1[0],
                    result.frequency_grid_cm1[-1],
                ),
                ylim=(-0.085, 1.05),
            )
            ax.grid(axis="y", alpha=0.18, linewidth=0.6)
            fig.tight_layout()
            plot_files = _save_figure_bundle(fig, output, "ir_spectrum")
            plt.close(fig)

    summary = {
        "dimensionality": 0,
        "response_kind": "molecular dipole derivative",
        "modes": len(result.mode_numbers),
        "dipole_derivative_unit": "Debye per Angstrom sqrt(amu)",
        "activity_unit": "Debye^2 per Angstrom^2 amu",
        "normalization": "maximum activity equals one",
        "files": {
            "modes": str(modes_path.resolve()),
            "spectrum": str(spectrum_path.resolve()),
            "plot": plot_files.get("plot"),
            "plot_pdf": plot_files.get("plot_pdf"),
            "plot_svg": plot_files.get("plot_svg"),
        },
    }
    (output / "ir_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    return summary


def _flatten_structure_positions(
    element_symbols: Sequence[str],
    element_coordinates: dict[str, list[list[float]]],
) -> np.ndarray:
    return np.asarray(
        [
            coordinate
            for symbol in element_symbols
            for coordinate in element_coordinates[symbol]
        ],
        dtype=float,
    )


def _write_displaced_stru(
    source: Path,
    destination: Path,
    displacement_angstrom: np.ndarray,
) -> None:
    (
        lattice_constant,
        lattice_vectors,
        element_symbols,
        element_atomnumber,
        coordinate_type,
        element_coordinates,
        element_movements,
        element_magnetisms,
        _element_mass,
        _element_pp,
        _element_orb,
    ) = stru_analyzer(str(source))
    cell = (
        np.asarray(lattice_vectors, dtype=float)
        * float(lattice_constant)
        * BOHR_ANGSTROM
    )
    positions = _flatten_structure_positions(element_symbols, element_coordinates)
    if positions.shape != displacement_angstrom.shape:
        raise ValueError(
            f"STRU has {len(positions)} atoms but mode has "
            f"{len(displacement_angstrom)} atoms"
        )
    if coordinate_type == "Direct":
        fractional = positions
    else:
        cartesian = positions * float(lattice_constant) * BOHR_ANGSTROM
        fractional = cartesian @ np.linalg.inv(cell)
    displaced = (fractional + displacement_angstrom @ np.linalg.inv(cell)) % 1.0

    text = source.read_text(encoding="utf-8")
    match = re.search(r"(?m)^ATOMIC_POSITIONS\s*$", text)
    if not match:
        raise ValueError(f"ATOMIC_POSITIONS not found in {source}")
    prefix = text[: match.end()]

    lines = [prefix, "\nDirect\n"]
    atom_index = 0
    for symbol in element_symbols:
        lines.extend([f"\n{symbol}\n", "0.0\n", f"{element_atomnumber[symbol]}\n"])
        for local_index in range(element_atomnumber[symbol]):
            position = displaced[atom_index]
            movement = element_movements[symbol][local_index]
            magnetism = element_magnetisms[symbol][local_index]
            line = "  ".join(f"{value:.12f}" for value in position)
            line += " m " + " ".join(str(int(value)) for value in movement)
            if isinstance(magnetism, (int, float)) and abs(float(magnetism)) > 0:
                line += f" mag {float(magnetism):g}"
            elif isinstance(magnetism, list) and any(
                abs(float(value)) > 0 for value in magnetism
            ):
                line += " mag " + " ".join(f"{float(value):g}" for value in magnetism)
            lines.append(line + "\n")
            atom_index += 1
    destination.write_text("".join(lines), encoding="utf-8")


def prepare_raman_displacements(
    stru_path: str | Path,
    modes: GammaModes,
    outdir: str | Path,
    *,
    amplitude: float = 0.02,
    mode_numbers: Optional[Sequence[int]] = None,
    acoustic_cutoff_cm1: float = 5.0,
    copy_files: Optional[Sequence[str | Path]] = None,
) -> dict:
    """Write central-difference structures displaced along Gamma modes.

    ``amplitude`` is the normal-coordinate step in Angstrom sqrt(amu), so each
    atomic displacement is ``amplitude * eigenvector / sqrt(mass)``.
    """

    source = Path(stru_path).resolve()
    output = Path(outdir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    indices = _selected_mode_indices(
        modes.frequencies_cm1, mode_numbers, acoustic_cutoff_cm1
    )
    eigenvectors = _mode_phase_real(modes.eigenvectors)
    mass_weighted = eigenvectors / np.sqrt(modes.masses_amu)[None, :, None]
    copied_files = [Path(path).resolve() for path in (copy_files or [])]

    entries = []
    for index in indices:
        mode_dir = output / f"mode-{index + 1:04d}"
        entry = {
            "mode": int(index + 1),
            "frequency_cm-1": float(modes.frequencies_cm1[index]),
            "amplitude_A_sqrt_amu": float(amplitude),
            "plus": str((mode_dir / "plus").resolve()),
            "minus": str((mode_dir / "minus").resolve()),
        }
        for sign_name, sign in (("plus", 1.0), ("minus", -1.0)):
            target = mode_dir / sign_name
            target.mkdir(parents=True, exist_ok=True)
            displacement = sign * float(amplitude) * mass_weighted[index]
            _write_displaced_stru(source, target / "STRU", displacement)
            for file_path in copied_files:
                if file_path.is_file():
                    shutil.copy2(file_path, target / file_path.name)
                elif file_path.is_dir():
                    shutil.copytree(
                        file_path,
                        target / file_path.name,
                        dirs_exist_ok=True,
                    )
                else:
                    raise FileNotFoundError(
                        f"Raman input to copy was not found: {file_path}"
                    )
        entries.append(entry)

    manifest = {
        "schema": 1,
        "source_stru": str(source),
        "qpoints_atoms": len(modes.masses_amu),
        "amplitude_A_sqrt_amu": float(amplitude),
        "modes": entries,
    }
    manifest_path = output / "raman_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def collect_raman_tensors(
    raman_dir: str | Path,
    *,
    dimensionality: int = 3,
    cell_height_angstrom: Optional[float] = None,
    cell_volume_angstrom3: Optional[float] = None,
) -> tuple[np.ndarray, np.ndarray, str]:
    """Central-difference PYATB dielectric tensors in ``mode-*/plus|minus``."""

    root = Path(raman_dir)
    manifest_path = root / "raman_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    amplitude = float(manifest["amplitude_A_sqrt_amu"])
    mode_numbers: list[int] = []
    tensors: list[np.ndarray] = []
    for entry in manifest["modes"]:
        plus_root = Path(entry["plus"])
        minus_root = Path(entry["minus"])
        plus, plus_source = read_static_dielectric(
            plus_root / "pyatb" if (plus_root / "pyatb").is_dir() else plus_root
        )
        minus, minus_source = read_static_dielectric(
            minus_root / "pyatb" if (minus_root / "pyatb").is_dir() else minus_root
        )
        derivative = (plus - minus) / (2.0 * amplitude)
        if dimensionality == 2:
            if cell_height_angstrom is None:
                raise ValueError(
                    "cell_height_angstrom is required for 2D Raman collection"
                )
            derivative = derivative * float(cell_height_angstrom)
        elif dimensionality == 0:
            if cell_volume_angstrom3 is None:
                raise ValueError(
                    "cell_volume_angstrom3 is required for molecular Raman collection"
                )
            derivative = (
                derivative * float(cell_volume_angstrom3) / (4.0 * math.pi)
            )
        elif dimensionality != 3:
            raise ValueError("dimensionality must be 0, 2, or 3")
        mode_numbers.append(int(entry["mode"]))
        tensors.append(0.5 * (derivative + derivative.T))
        entry["plus_dielectric"] = str(plus_source)
        entry["minus_dielectric"] = str(minus_source)
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    kind = {
        0: "molecular polarizability derivative (Angstrom^3 per Angstrom sqrt(amu))",
        2: "2D sheet susceptibility derivative",
        3: "dielectric tensor derivative",
    }[dimensionality]
    tensor_array = np.asarray(tensors, dtype=float)
    np.save(root / "raman_tensors.npy", tensor_array)
    (root / "raman_tensors.json").write_text(
        json.dumps(
            {
                "mode_numbers": mode_numbers,
                "tensors": tensor_array.tolist(),
                "tensor_kind": kind,
                "cell_volume_A3": (
                    float(cell_volume_angstrom3)
                    if dimensionality == 0
                    else None
                ),
                "conversion": (
                    "dalpha/dQ = V/(4*pi) * d(epsilon_r)/dQ"
                    if dimensionality == 0
                    else None
                ),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return np.asarray(mode_numbers, dtype=int), tensor_array, kind


def load_raman_tensors(path: str | Path) -> tuple[np.ndarray, Optional[np.ndarray], str]:
    tensor_path = Path(path)
    if tensor_path.suffix.lower() == ".npy":
        tensors = np.asarray(np.load(tensor_path), dtype=float)
        mode_numbers = np.arange(1, len(tensors) + 1, dtype=int)
        return mode_numbers, tensors, "user-supplied Raman tensor"
    data = yaml.safe_load(tensor_path.read_text(encoding="utf-8"))
    if isinstance(data, dict):
        raw = data.get("tensors") or data.get("raman_tensors")
        numbers = data.get("mode_numbers")
        kind = str(data.get("tensor_kind", "user-supplied Raman tensor"))
    else:
        raw = data
        numbers = None
        kind = "user-supplied Raman tensor"
    tensors = np.asarray(raw, dtype=float)
    mode_numbers = (
        np.asarray(numbers, dtype=int)
        if numbers is not None
        else np.arange(1, len(tensors) + 1, dtype=int)
    )
    return mode_numbers, tensors, kind


def calculate_raman_spectrum(
    modes: GammaModes,
    mode_numbers: Sequence[int],
    tensors: np.ndarray,
    *,
    tensor_kind: str = "Raman tensor",
    temperature_K: float = 300.0,
    laser_nm: float = 532.0,
    broadening_cm1: float = 8.0,
    max_frequency_cm1: Optional[float] = None,
    points: int = 2001,
) -> RamanSpectrumResult:
    """Calculate powder-averaged, non-resonant Placzek Raman intensities."""

    numbers = np.asarray(mode_numbers, dtype=int)
    raman = np.asarray(tensors, dtype=float)
    if raman.shape != (len(numbers), 3, 3):
        raise ValueError(f"Raman tensors have invalid shape: {raman.shape}")
    indices = numbers - 1
    if np.any(indices < 0) or np.any(indices >= len(modes.frequencies_cm1)):
        raise IndexError("Raman mode number is outside qpoints.yaml")
    frequencies = modes.frequencies_cm1[indices]
    symmetric = 0.5 * (raman + np.swapaxes(raman, 1, 2))
    alpha = np.trace(symmetric, axis1=1, axis2=2) / 3.0
    gamma2 = 0.5 * (
        (symmetric[:, 0, 0] - symmetric[:, 1, 1]) ** 2
        + (symmetric[:, 1, 1] - symmetric[:, 2, 2]) ** 2
        + (symmetric[:, 2, 2] - symmetric[:, 0, 0]) ** 2
        + 6.0
        * (
            symmetric[:, 0, 1] ** 2
            + symmetric[:, 1, 2] ** 2
            + symmetric[:, 0, 2] ** 2
        )
    )
    placzek = 45.0 * alpha * alpha + 7.0 * gamma2
    depolarization = 3.0 * gamma2 / np.maximum(
        45.0 * alpha * alpha + 4.0 * gamma2, np.finfo(float).tiny
    )

    positive = frequencies > 1.0e-12
    laser_cm1 = 1.0e7 / float(laser_nm)
    bose = np.zeros_like(frequencies)
    exponent = (
        PLANCK
        * SPEED_OF_LIGHT
        * frequencies[positive]
        * 100.0
        / (BOLTZMANN * float(temperature_K))
    )
    bose[positive] = 1.0 / np.expm1(exponent)
    activities = np.zeros_like(frequencies)
    activities[positive] = (
        np.maximum(laser_cm1 - frequencies[positive], 0.0) ** 4
        * (bose[positive] + 1.0)
        * placzek[positive]
        / frequencies[positive]
    )
    maximum = float(np.max(activities)) if len(activities) else 0.0
    if maximum > 0.0:
        activities /= maximum

    upper = max_frequency_cm1
    if upper is None:
        upper = max(100.0, float(np.max(frequencies)) + 5.0 * broadening_cm1)
    grid = np.linspace(0.0, float(upper), int(points))
    spectrum = _lorentzian(grid, frequencies, broadening_cm1) @ activities
    if float(np.max(spectrum)) > 0.0:
        spectrum /= float(np.max(spectrum))

    return RamanSpectrumResult(
        mode_numbers=numbers,
        frequencies_cm1=frequencies,
        tensors=symmetric,
        activities=activities,
        depolarization_ratios=depolarization,
        frequency_grid_cm1=grid,
        spectrum=spectrum,
        tensor_kind=tensor_kind,
    )


def write_raman_outputs(
    outdir: str | Path,
    result: RamanSpectrumResult,
    *,
    plot: bool = True,
) -> dict:
    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    modes_path = output / "raman_modes.csv"
    with modes_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "mode",
                "frequency_cm-1",
                "activity_normalized",
                "depolarization_ratio",
                "Rxx",
                "Rxy",
                "Rxz",
                "Ryx",
                "Ryy",
                "Ryz",
                "Rzx",
                "Rzy",
                "Rzz",
            ]
        )
        for number, frequency, activity, ratio, tensor in zip(
            result.mode_numbers,
            result.frequencies_cm1,
            result.activities,
            result.depolarization_ratios,
            result.tensors,
        ):
            writer.writerow(
                [
                    int(number),
                    float(frequency),
                    float(activity),
                    float(ratio),
                    *tensor.reshape(9).tolist(),
                ]
            )
    tensor_path = output / "raman_tensors.npy"
    np.save(tensor_path, result.tensors)
    spectrum_path = output / "raman_spectrum.dat"
    np.savetxt(
        spectrum_path,
        np.column_stack([result.frequency_grid_cm1, result.spectrum]),
        header="frequency_cm-1 intensity_normalized",
    )

    plot_files: dict[str, str] = {}
    if plot:
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        with mpl.rc_context(
            {
                "font.family": "sans-serif",
                "font.sans-serif": ["Arial", "DejaVu Sans"],
                "font.size": 9,
                "axes.spines.right": True,
                "axes.spines.top": True,
                "axes.linewidth": 0.8,
                "xtick.direction": "in",
                "ytick.direction": "in",
                "xtick.top": True,
                "ytick.right": True,
                "legend.frameon": False,
                "pdf.fonttype": 42,
                "svg.fonttype": "none",
            }
        ):
            fig, ax = plt.subplots(figsize=(7.2, 4.5))
            ax.plot(
                result.frequency_grid_cm1,
                result.spectrum,
                color="#2f6b9a",
                linewidth=1.7,
            )
            ax.vlines(
                result.frequencies_cm1,
                -0.075,
                -0.075 + 0.06 * result.activities,
                color="#b2472f",
                linewidth=0.9,
            )
            ax.set_xlabel(r"Raman shift (cm$^{-1}$)")
            ax.set_ylabel("Normalized Raman intensity")
            ax.set_xlim(
                result.frequency_grid_cm1[0], result.frequency_grid_cm1[-1]
            )
            ax.set_ylim(-0.085, 1.05)
            ax.grid(axis="y", alpha=0.18, linewidth=0.6)
            fig.tight_layout()
            plot_files = _save_figure_bundle(fig, output, "raman_spectrum")
            plt.close(fig)

    summary = {
        "tensor_kind": result.tensor_kind,
        "modes": len(result.mode_numbers),
        "files": {
            "modes": str(modes_path.resolve()),
            "tensors": str(tensor_path.resolve()),
            "spectrum": str(spectrum_path.resolve()),
            "plot": plot_files.get("plot"),
            "plot_pdf": plot_files.get("plot_pdf"),
            "plot_svg": plot_files.get("plot_svg"),
        },
    }
    (output / "raman_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    return summary
