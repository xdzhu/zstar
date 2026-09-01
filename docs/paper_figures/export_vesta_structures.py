"""Export manuscript ABACUS structures as VASP5 files for VESTA rendering."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np


BOHR_TO_ANGSTROM = 0.529177210903


def read_abacus_stru(path: Path) -> tuple[np.ndarray, list[str], np.ndarray]:
    lines = [line.split("#", 1)[0].strip() for line in path.read_text().splitlines()]
    scale_index = lines.index("LATTICE_CONSTANT")
    scale = float(lines[scale_index + 1]) * BOHR_TO_ANGSTROM
    vectors_index = lines.index("LATTICE_VECTORS")
    lattice = scale * np.array(
        [[float(value) for value in lines[vectors_index + offset].split()[:3]] for offset in range(1, 4)]
    )

    positions_index = lines.index("ATOMIC_POSITIONS")
    if lines[positions_index + 1].lower() != "direct":
        raise ValueError(f"Only Direct coordinates are supported: {path}")
    species: list[str] = []
    fractional: list[list[float]] = []
    cursor = positions_index + 2
    while cursor < len(lines):
        if not lines[cursor]:
            cursor += 1
            continue
        element = lines[cursor].split()[0]
        cursor += 1
        while cursor < len(lines) and not lines[cursor]:
            cursor += 1
        cursor += 1  # species magnetism line
        while cursor < len(lines) and not lines[cursor]:
            cursor += 1
        count = int(lines[cursor].split()[0])
        cursor += 1
        for _ in range(count):
            while cursor < len(lines) and not lines[cursor]:
                cursor += 1
            fractional.append([float(value) for value in lines[cursor].split()[:3]])
            species.append(element)
            cursor += 1
    return lattice, species, np.asarray(fractional, dtype=float)


def ordered_counts(species: list[str]) -> tuple[list[str], list[int]]:
    order: list[str] = []
    counts: list[int] = []
    for element in species:
        if not order or element != order[-1]:
            if element in order:
                raise ValueError(f"Species {element} is not contiguous")
            order.append(element)
            counts.append(1)
        else:
            counts[-1] += 1
    return order, counts


def write_vasp(path: Path, title: str, lattice: np.ndarray, species: list[str], fractional: np.ndarray) -> None:
    elements, counts = ordered_counts(species)
    lines = [title, "1.0"]
    lines.extend("  " + "  ".join(f"{value:16.10f}" for value in vector) for vector in lattice)
    lines.append("  " + "  ".join(elements))
    lines.append("  " + "  ".join(str(count) for count in counts))
    lines.append("Direct")
    lines.extend("  " + "  ".join(f"{value:16.10f}" for value in position) for position in fractional)
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    root = Path(__file__).resolve().parent
    data = root / "source_data"
    output = root / "vesta_structures"
    output.mkdir(parents=True, exist_ok=True)
    structures = {
        "CH4_molecule.vasp": (data / "molecular/ch4/STRU", "CH4 molecule; Figure 3a and molecular APT panel b"),
        "GaAs_nanowire.vasp": (data / "gaas_nanowire/STRU", "H-passivated wurtzite GaAs nanowire; Figure 3d"),
        "MoS2_monolayer.vasp": (data / "mos2/STRU", "Monolayer 2H-MoS2; Figure 3g"),
        "PbTiO3_tetragonal.vasp": (data / "pto/STRU", "Tetragonal P4mm PbTiO3; archived spectroscopy validation"),
        "BaTiO3_cubic.vasp": (data / "bec_structures/bto_cubic/STRU", "Cubic Pm-3m BaTiO3; bulk BEC panel a"),
        "HfO2_tetragonal.vasp": (data / "bec_structures/hfo2_tetragonal/STRU", "Tetragonal P42/nmc HfO2; Figure 6j and bulk BEC panel b"),
        "hBN_monolayer.vasp": (data / "hbn/STRU", "Monolayer hBN; 2D BEC panel a"),
        "alpha-In2Se3_monolayer.vasp": (data / "bec_structures/in2se3_monolayer/STRU", "Ferroelectric alpha-In2Se3 monolayer; 2D BEC panel b"),
        "H2O_molecule.vasp": (data / "bec_structures/h2o_molecule/STRU", "H2O molecule; molecular APT panel a"),
    }

    records = []
    for name, (source, title) in structures.items():
        lattice, species, fractional = read_abacus_stru(source)
        destination = output / name
        write_vasp(destination, title, lattice, species, fractional)
        records.append(
            {
                "output": destination.name,
                "source": str(source.relative_to(root)).replace("\\", "/"),
                "atoms": len(species),
                "elements": ordered_counts(species)[0],
                "cell_volume_A3": float(abs(np.linalg.det(lattice))),
                "source_sha256": sha256(source),
                "output_sha256": sha256(destination),
            }
        )

    (output / "manifest.json").write_text(json.dumps({"schema": 1, "structures": records}, indent=2) + "\n")
    print(json.dumps({"output": str(output), "count": len(records)}, indent=2))


if __name__ == "__main__":
    main()
