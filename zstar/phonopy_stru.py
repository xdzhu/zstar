"""Write ABACUS structures in the strict form expected by Phonopy."""

from __future__ import annotations

from pathlib import Path

from .stru_analyzer import stru_analyzer


def write_phonopy_compatible_stru(
    source: str | Path,
    destination: str | Path,
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
            if isinstance(magnetism, (int, float)) and abs(float(magnetism)) > 0:
                line += f" mag {float(magnetism):g}"
            elif isinstance(magnetism, list) and any(
                abs(float(value)) > 0 for value in magnetism
            ):
                line += " mag " + " ".join(
                    f"{float(value):g}" for value in magnetism
                )
            lines.append(line + "\n")

    Path(destination).write_text("".join(lines), encoding="utf-8")
