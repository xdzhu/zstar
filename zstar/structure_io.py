"""Small structure readers used by calculator-independent ZStar paths.

The module intentionally covers only the data needed by workflow orchestration:
cell vectors, fractional coordinates, element labels, and cell volume.  It is
not intended to replace a full crystallographic parser.  In particular, rich
VASP XML, CHGCAR, and POTCAR parsing remains in the optional VASP adapter.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np

try:
    import spglib
except Exception:  # pragma: no cover - declared package dependency
    spglib = None


_FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


@dataclass(frozen=True)
class StructureData:
    """Minimal periodic structure representation used by ZStar."""

    lattice_angstrom: np.ndarray
    positions_fractional: np.ndarray
    symbols: tuple[str, ...]

    @property
    def volume(self) -> float:
        return abs(float(np.linalg.det(self.lattice_angstrom)))

    @property
    def cart_coords(self) -> np.ndarray:
        return np.asarray(self.positions_fractional) @ np.asarray(self.lattice_angstrom)


def _dataset_value(dataset, name: str):
    return getattr(dataset, name) if hasattr(dataset, name) else dataset[name]


def _float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def read_poscar(path: str | Path) -> StructureData:
    """Read the common POSCAR/CONTCAR format without a pymatgen dependency."""

    source = Path(path)
    lines = source.read_text(encoding="utf-8", errors="ignore").splitlines()
    if len(lines) < 8:
        raise ValueError(f"POSCAR is too short: {source}")
    scale = _float(lines[1].split()[0])
    raw_lattice = np.asarray(
        [[_float(token) for token in lines[index].split()[:3]] for index in range(2, 5)],
        dtype=float,
    )
    if scale < 0:
        raw_volume = abs(float(np.linalg.det(raw_lattice)))
        if raw_volume <= 0:
            raise ValueError(f"POSCAR has a zero lattice volume: {source}")
        scale = (-scale / raw_volume) ** (1.0 / 3.0)
    lattice = raw_lattice * scale

    cursor = 5
    species_tokens = lines[cursor].split()
    try:
        counts = [int(token) for token in species_tokens]
    except ValueError:
        symbols = species_tokens
        cursor += 1
        counts = [int(token) for token in lines[cursor].split()]
    if not counts or any(count < 0 for count in counts):
        raise ValueError(f"POSCAR has invalid atom counts: {source}")
    if len(symbols) != len(counts):
        raise ValueError(f"POSCAR species/count mismatch: {source}")
    cursor += 1
    if cursor < len(lines) and lines[cursor].strip().lower().startswith("s"):
        cursor += 1
    if cursor >= len(lines):
        raise ValueError(f"POSCAR has no coordinate mode: {source}")
    mode = lines[cursor].strip().lower()
    cursor += 1
    n_atoms = sum(counts)
    coordinates: list[list[float]] = []
    for line in lines[cursor:]:
        fields = line.split()
        if len(fields) < 3:
            continue
        try:
            coordinates.append([_float(token) for token in fields[:3]])
        except ValueError:
            continue
        if len(coordinates) == n_atoms:
            break
    if len(coordinates) != n_atoms:
        raise ValueError(f"POSCAR contains {len(coordinates)} coordinates, expected {n_atoms}: {source}")
    positions = np.asarray(coordinates, dtype=float)
    if mode.startswith("c") or mode.startswith("k"):
        positions = positions @ np.linalg.inv(lattice)
    labels = tuple(symbol for symbol, count in zip(symbols, counts) for _ in range(count))
    return StructureData(lattice, positions, labels)


def read_abacus_stru(path: str | Path) -> StructureData:
    """Read the lattice and atomic positions from an ABACUS ``STRU`` file."""

    source = Path(path)
    lines = source.read_text(encoding="utf-8", errors="ignore").splitlines()
    try:
        constant_index = next(index for index, line in enumerate(lines) if line.strip().upper() == "LATTICE_CONSTANT")
        lattice_constant = _float(lines[constant_index + 1].split()[0])
        vector_index = next(index for index, line in enumerate(lines) if line.strip().upper() == "LATTICE_VECTORS")
        lattice = np.asarray(
            [[_float(token) for token in lines[vector_index + 1 + offset].split()[:3]] for offset in range(3)],
            dtype=float,
        ) * lattice_constant / 1.889726125
        species_index = next(index for index, line in enumerate(lines) if line.strip().upper() == "ATOMIC_SPECIES")
        position_index = next(index for index, line in enumerate(lines) if line.strip().upper() == "ATOMIC_POSITIONS")
    except (StopIteration, IndexError, ValueError) as exc:
        raise ValueError(f"Cannot parse ABACUS STRU header: {source}") from exc

    species: list[str] = []
    for line in lines[species_index + 1 : position_index]:
        fields = line.split()
        if len(fields) >= 2 and not line.lstrip().startswith("#"):
            try:
                _float(fields[1])
            except ValueError:
                continue
            species.append(fields[0])
    if not species:
        raise ValueError(f"STRU contains no atomic species: {source}")
    mode = lines[position_index + 1].strip().lower()
    cursor = position_index + 2
    positions: list[list[float]] = []
    labels: list[str] = []
    while cursor < len(lines):
        label = lines[cursor].split("#", 1)[0].strip()
        if not label:
            cursor += 1
            continue
        if label not in species:
            break
        if cursor + 2 >= len(lines):
            raise ValueError(f"Incomplete atomic block for {label} in {source}")
        try:
            count = int(lines[cursor + 2].split("#", 1)[0].strip().split()[0])
        except (IndexError, ValueError) as exc:
            raise ValueError(f"Invalid atom count for {label} in {source}") from exc
        for raw in lines[cursor + 3 : cursor + 3 + count]:
            fields = raw.split()
            if len(fields) < 3:
                raise ValueError(f"Invalid coordinate line for {label} in {source}")
            positions.append([_float(token) for token in fields[:3]])
            labels.append(label)
        cursor += 3 + count
    coordinates = np.asarray(positions, dtype=float)
    if mode.startswith("cart"):
        coordinates = coordinates @ np.linalg.inv(lattice)
    if len(labels) == 0:
        raise ValueError(f"STRU contains no atomic coordinates: {source}")
    return StructureData(lattice, coordinates, tuple(labels))


def read_structure(path: str | Path) -> StructureData:
    """Read POSCAR-like or ABACUS structure files based on their content/name."""

    source = Path(path)
    if source.name.upper() == "STRU" or source.suffix.lower() in {".stru", ".abacus"}:
        return read_abacus_stru(source)
    return read_poscar(source)


def wyckoff_summary(path: str | Path, *, symprec: float = 1.0e-3) -> dict:
    """Return space-group and per-site Wyckoff data using only ``spglib``.

    Chemical symbols are mapped to stable integer type labels because spglib
    needs species identity, not atomic numbers, for symmetry detection.
    """

    if spglib is None:
        raise RuntimeError("spglib is required for Wyckoff analysis")
    if not np.isfinite(float(symprec)) or float(symprec) <= 0.0:
        raise ValueError("symprec must be finite and positive")
    structure = read_structure(path)
    type_ids: dict[str, int] = {}
    numbers: list[int] = []
    for symbol in structure.symbols:
        type_ids.setdefault(symbol, len(type_ids) + 1)
        numbers.append(type_ids[symbol])
    dataset = spglib.get_symmetry_dataset(
        (
            np.asarray(structure.lattice_angstrom, dtype=float),
            np.mod(np.asarray(structure.positions_fractional, dtype=float), 1.0),
            numbers,
        ),
        symprec=float(symprec),
    )
    if dataset is None:
        raise RuntimeError(
            f"spglib could not identify a symmetry dataset for {Path(path).resolve()}"
        )
    wyckoffs = tuple(str(value) for value in _dataset_value(dataset, "wyckoffs"))
    equivalent = tuple(
        int(value) for value in _dataset_value(dataset, "equivalent_atoms")
    )
    site_symmetry = tuple(
        str(value) for value in _dataset_value(dataset, "site_symmetry_symbols")
    )
    sites = []
    for index, (symbol, letter, representative, symmetry, position) in enumerate(
        zip(
            structure.symbols,
            wyckoffs,
            equivalent,
            site_symmetry,
            structure.positions_fractional,
        ),
        start=1,
    ):
        sites.append(
            {
                "index": index,
                "symbol": symbol,
                "wyckoff": letter,
                "site_symmetry": symmetry,
                "representative": representative + 1,
                "fractional": tuple(float(value % 1.0) for value in position),
            }
        )
    return {
        "space_group": str(_dataset_value(dataset, "international")),
        "number": int(_dataset_value(dataset, "number")),
        "hall_number": int(_dataset_value(dataset, "hall_number")),
        "volume_angstrom3": structure.volume,
        "sites": sites,
    }


def format_wyckoff_summary(summary: dict) -> str:
    """Format :func:`wyckoff_summary` for the command-line interface."""

    lines = [
        f"Space group: {summary['space_group']} (No. {summary['number']})",
        f"Cell volume: {summary['volume_angstrom3']:.8f} Angstrom^3",
        "Index  Element  Wyckoff  Site symmetry  Representative  Fractional coordinates",
    ]
    for site in summary["sites"]:
        coordinates = " ".join(f"{value:.8f}" for value in site["fractional"])
        lines.append(
            f"{site['index']:5d}  {site['symbol']:<7s}  {site['wyckoff']:<7s}  "
            f"{site['site_symmetry']:<13s}  {site['representative']:14d}  {coordinates}"
        )
    return "\n".join(lines)


def write_poscar(path: str | Path, structure: StructureData, comment: str = "Generated by ZStar") -> Path:
    """Write a minimal, portable POSCAR from :class:`StructureData`."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    ordered_symbols: list[str] = []
    for symbol in structure.symbols:
        if symbol not in ordered_symbols:
            ordered_symbols.append(symbol)
    lines = [comment, "1.0"]
    lines.extend("  ".join(f"{value:.16f}" for value in row) for row in structure.lattice_angstrom)
    lines.append("  ".join(ordered_symbols))
    lines.append("  ".join(str(structure.symbols.count(symbol)) for symbol in ordered_symbols))
    lines.append("Direct")
    for symbol in ordered_symbols:
        for position, site_symbol in zip(structure.positions_fractional, structure.symbols):
            if site_symbol == symbol:
                lines.append("  ".join(f"{value:.16f}" for value in position))
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return destination
