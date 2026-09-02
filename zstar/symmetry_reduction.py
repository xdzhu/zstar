"""Structure-level symmetry reduction for finite-displacement workflows."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .verify_born_symmetry import read_abacus_stru

try:
    import spglib
except Exception:  # pragma: no cover - dependency is declared by the package
    spglib = None


@dataclass(frozen=True)
class SymmetryReduction:
    """Symmetry information needed to choose finite-displacement atoms."""

    representatives: tuple[int, ...]
    equivalent_atoms: tuple[int, ...]
    symbols: tuple[str, ...]
    space_group: str | None
    hall_number: int | None
    dimensionality: int
    engine: str

    @property
    def atom_count(self) -> int:
        return len(self.symbols)


def _dataset_value(dataset, name: str):
    return getattr(dataset, name) if hasattr(dataset, name) else dataset[name]


def reduce_abacus_atoms(
    structure: str | Path,
    *,
    symprec: float = 1.0e-3,
    dimensionality: int = 3,
) -> SymmetryReduction:
    """Return one-based independent atom indices for an ABACUS ``STRU``.

    ``spglib`` operates on a periodic cell. For a molecular workflow
    (``dimensionality=0``), no periodic symmetry reduction is applied here;
    molecular point-group expansion remains the responsibility of the APT
    post-processor. This prevents a molecule's vacuum box from being
    mistaken for a crystal symmetry.
    """

    dimension = int(dimensionality)
    if dimension not in {0, 1, 2, 3}:
        raise ValueError("dimensionality must be 0, 1, 2, or 3")
    if not np.isfinite(float(symprec)) or float(symprec) <= 0.0:
        raise ValueError("symprec must be finite and positive")
    lattice, fractional, symbols = read_abacus_stru(str(structure))
    labels = tuple(str(symbol) for symbol in symbols)
    all_atoms = tuple(range(1, len(labels) + 1))
    if dimension == 0:
        return SymmetryReduction(
            representatives=all_atoms,
            equivalent_atoms=tuple(index - 1 for index in all_atoms),
            symbols=labels,
            space_group=None,
            hall_number=None,
            dimensionality=dimension,
            engine="none-molecular",
        )
    if spglib is None:
        raise RuntimeError("spglib is required for periodic symmetry reduction")

    # Atomic numbers are unnecessary for symmetry; stable per-element labels
    # preserve chemical species without introducing an ASE dependency.
    type_ids: dict[str, int] = {}
    numbers = []
    for symbol in labels:
        type_ids.setdefault(symbol, len(type_ids) + 1)
        numbers.append(type_ids[symbol])
    dataset = spglib.get_symmetry_dataset(
        (
            np.asarray(lattice, dtype=float),
            np.asarray(fractional, dtype=float),
            numbers,
        ),
        symprec=float(symprec),
    )
    if dataset is None:
        raise RuntimeError(
            "spglib could not identify a symmetry dataset for "
            f"{Path(structure).resolve()}. Rerun with --all to generate "
            "all-atom displacements, or adjust --symmprec."
        )
    equivalent = tuple(
        int(value) for value in _dataset_value(dataset, "equivalent_atoms")
    )
    if len(equivalent) != len(labels):
        raise RuntimeError("spglib returned an invalid equivalent-atom mapping")
    representatives = tuple(
        index + 1
        for index, representative in enumerate(equivalent)
        if representative == index
    )
    if not representatives:
        raise RuntimeError("spglib returned no symmetry representatives")
    return SymmetryReduction(
        representatives=representatives,
        equivalent_atoms=equivalent,
        symbols=labels,
        space_group=str(_dataset_value(dataset, "international")),
        hall_number=int(_dataset_value(dataset, "hall_number")),
        dimensionality=dimension,
        engine="spglib",
    )


def write_reduction_report(path: str | Path, result: SymmetryReduction) -> Path:
    """Write a compact human-readable record of the selected representatives."""

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# ZStar symmetry reduction",
        f"engine = {result.engine}",
        f"dimensionality = {result.dimensionality}",
        f"space_group = {result.space_group or 'none'}",
        f"hall_number = {result.hall_number if result.hall_number is not None else 'none'}",
        "# one-based atom index, element, equivalent-atom representative (zero-based)",
    ]
    lines.extend(
        f"{index:4d} {symbol:<4s} {representative:4d}"
        for index, (symbol, representative) in enumerate(
            zip(result.symbols, result.equivalent_atoms), start=1
        )
    )
    lines.append("representatives = " + " ".join(map(str, result.representatives)))
    target.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return target
