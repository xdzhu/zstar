"""Write ABACUS structures in the strict form expected by Phonopy."""

from __future__ import annotations

from pathlib import Path

from .stru_analyzer import stru_analyzer


def write_phonopy_compatible_stru(
    source: str | Path,
    destination: str | Path,
    *,
    include_magnetism: bool = True,
) -> None:
    """Normalize optional ABACUS coordinate fields for Phonopy."""

    (
        lattice_constant,
        lattice_vectors,
        element_symbols,
        element_atomnumber,
        coordinate_type,
        element_coordinates,
        element_movements,
        element_magnetisms,
        element_mass,
        element_pp,
        element_orb,
    ) = stru_analyzer(str(source))

    lines = ["ATOMIC_SPECIES\n"]
    for symbol, mass, pp in zip(element_symbols, element_mass, element_pp):
        lines.append(f"{symbol} {mass} {pp}\n")
    if any(element_orb):
        lines.append("\nNUMERICAL_ORBITAL\n")
        lines.extend(f"{orb}\n" for orb in element_orb)
    lines.extend(
        [
            "\nLATTICE_CONSTANT\n",
            f"{lattice_constant:.12f}\n",
            "\nLATTICE_VECTORS\n",
        ]
    )
    lines.extend(
        " ".join(f"{value:.12f}" for value in vector) + "\n"
        for vector in lattice_vectors
    )
    lines.extend(["\nATOMIC_POSITIONS\n", f"{coordinate_type}\n"])

    magnetic_values = [
        value
        for symbol in element_symbols
        for value in element_magnetisms[symbol]
    ]
    magnetic_width = 3 if any(isinstance(value, list) for value in magnetic_values) else 1
    has_magnetism = any(
        any(abs(float(component)) > 0 for component in value)
        if isinstance(value, list)
        else abs(float(value)) > 0
        for value in magnetic_values
    )

    for symbol in element_symbols:
        lines.extend(
            [f"\n{symbol}\n", "0.0\n", f"{element_atomnumber[symbol]}\n"]
        )
        for coords, movement, magnetism in zip(
            element_coordinates[symbol],
            element_movements[symbol],
            element_magnetisms[symbol],
        ):
            line = " ".join(f"{value:.12f}" for value in coords)
            line += " m " + " ".join(str(int(value)) for value in movement)
            # Phonopy converts the coordinate rows to a rectangular NumPy
            # array.  If only magnetic atoms carry a ``mag`` field, mixed
            # magnetic/non-magnetic structures produce ragged rows and fail
            # before symmetry analysis.  Once any moment is present, emit the
            # same-width magnetic field for every atom, including explicit 0s.
            if include_magnetism and has_magnetism:
                values = magnetism if isinstance(magnetism, list) else [magnetism]
                if magnetic_width == 3:
                    values = list(values[:3]) + [0.0] * (3 - len(values[:3]))
                else:
                    values = values[:1]
                line += " mag " + " ".join(f"{float(value):g}" for value in values)
            lines.append(line + "\n")

    Path(destination).write_text("".join(lines), encoding="utf-8")
