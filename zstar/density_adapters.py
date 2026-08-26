"""Calculator-specific exporters into ZStar's calculator-neutral cube contract."""

from __future__ import annotations

import json
from pathlib import Path
import re
from typing import Sequence


_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


def write_cube_sidecar(
    cube_path: str | Path,
    ionic_valence_charges: Sequence[float],
    *,
    backend: str,
    density_sign: str = "positive_electron_density",
) -> Path:
    cube = Path(cube_path).resolve()
    charges = [float(value) for value in ionic_valence_charges]
    if not charges or any(value <= 0.0 for value in charges):
        raise ValueError("ionic_valence_charges must be a non-empty positive sequence")
    target = cube.with_suffix(cube.suffix + ".zstar.json")
    target.write_text(
        json.dumps(
            {
                "schema": "zstar-cube",
                "schema_version": "1.0",
                "backend": str(backend),
                "cube": cube.name,
                "density_sign": density_sign,
                "ionic_valence_charges": charges,
                "length_unit": "bohr",
                "density_unit": "electron/bohr^3",
            },
            indent=2,
        ),
        encoding="utf-8",
        newline="\n",
    )
    return target


def vasp_chgcar_to_cube(
    chgcar_path: str | Path,
    output_path: str | Path,
    *,
    potcar_path: str | Path | None = None,
) -> Path:
    """Convert VASP CHGCAR to cube and record valence ionic charges."""

    from pymatgen.io.vasp.inputs import Potcar
    from pymatgen.io.vasp.outputs import Chgcar

    source = Path(chgcar_path).resolve()
    target = Path(output_path).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    charge = Chgcar.from_file(source)
    charge.to_cube(target, comment="ZStar electron density converted from VASP CHGCAR")
    if potcar_path is not None:
        potcar = Potcar.from_file(Path(potcar_path).resolve())
        valence_by_symbol = {
            str(item.element): float(item.nelectrons) for item in potcar
        }
        charges = [valence_by_symbol[str(site.specie.symbol)] for site in charge.structure]
        write_cube_sidecar(target, charges, backend="vasp")
    return target


def qe_pp_cube_input(
    *,
    prefix: str,
    outdir: str = "./tmp",
    output_cube: str = "charge-density.cube",
    filplot: str = "zstar_charge_density",
) -> str:
    """Return pp.x input for the total valence electron density in cube format."""

    return (
        "&inputpp\n"
        f"  prefix = {prefix!r},\n"
        f"  outdir = {outdir!r},\n"
        f"  filplot = {filplot!r},\n"
        "  plot_num = 0,\n"
        "/\n"
        "&plot\n"
        "  iflag = 3,\n"
        "  output_format = 6,\n"
        f"  fileout = {output_cube!r},\n"
        "/\n"
    )


def _qe_species_and_positions(text: str) -> tuple[dict[str, str], list[str]]:
    lines = text.splitlines()
    species: dict[str, str] = {}
    labels: list[str] = []
    for index, line in enumerate(lines):
        if re.match(r"^\s*ATOMIC_SPECIES\b", line, re.I):
            for raw in lines[index + 1 :]:
                fields = raw.split()
                if len(fields) < 3 or not re.fullmatch(_NUMBER, fields[1]):
                    break
                species[fields[0]] = fields[2]
        if re.match(r"^\s*ATOMIC_POSITIONS\b", line, re.I):
            for raw in lines[index + 1 :]:
                fields = raw.split()
                if len(fields) < 4 or not all(re.fullmatch(_NUMBER, item) for item in fields[1:4]):
                    break
                labels.append(fields[0])
    if not species or not labels:
        raise ValueError("QE input must contain ATOMIC_SPECIES and ATOMIC_POSITIONS")
    return species, labels


def _upf_z_valence(path: Path) -> float:
    text = path.read_text(encoding="utf-8", errors="ignore")[:100000]
    patterns = (
        rf"z_valence\s*=\s*[\"']({_NUMBER})",
        rf"Z\s*valence\s*=\s*({_NUMBER})",
    )
    for pattern in patterns:
        match = re.search(pattern, text, re.I)
        if match:
            return float(match.group(1).replace("D", "E").replace("d", "e"))
    raise ValueError(f"Cannot determine z_valence from UPF: {path}")


def write_qe_cube_sidecar(
    cube_path: str | Path,
    pw_input: str | Path,
    *,
    pseudo_dir: str | Path,
) -> Path:
    source = Path(pw_input).resolve()
    species, labels = _qe_species_and_positions(source.read_text(encoding="utf-8"))
    root = Path(pseudo_dir).resolve()
    values = {label: _upf_z_valence(root / filename) for label, filename in species.items()}
    return write_cube_sidecar(
        cube_path, [values[label] for label in labels], backend="qe"
    )


def cp2k_density_cube_block(*, stride: tuple[int, int, int] = (1, 1, 1)) -> str:
    if len(stride) != 3 or any(int(value) < 1 for value in stride):
        raise ValueError("stride must contain three positive integers")
    return (
        "&E_DENSITY_CUBE\n"
        f"  STRIDE {int(stride[0])} {int(stride[1])} {int(stride[2])}\n"
        "&END E_DENSITY_CUBE\n"
    )
