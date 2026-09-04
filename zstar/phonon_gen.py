"""Generate finite-displacement phonon directories with Phonopy."""

from __future__ import annotations

from pathlib import Path
import re
import shutil
import subprocess
from typing import Optional

import numpy as np


_ABACUS_SECTION_HEADERS = {
    "ATOMIC_SPECIES",
    "NUMERICAL_ORBITAL",
    "NUMERICAL_DESCRIPTOR",
    "LATTICE_CONSTANT",
    "LATTICE_VECTORS",
    "ATOMIC_POSITIONS",
}


def create_symlink(source: str | Path, link_name: str | Path) -> None:
    source_path = Path(source).resolve()
    link_path = Path(link_name)
    if link_path.is_symlink() or link_path.is_file():
        link_path.unlink()
    elif link_path.exists():
        raise IsADirectoryError(f"Cannot replace directory with symlink: {link_path}")
    link_path.symlink_to(source_path)


def _copy_optional_script(script: Optional[str | Path], target: Path) -> None:
    if not script:
        return
    source = Path(script)
    if source.is_file():
        shutil.copy2(source, target / source.name)


def _referenced_abacus_assets(structure: Path, workdir: Path) -> list[Path]:
    """Find distributable PP/orbital files referenced by an ABACUS STRU."""

    section = None
    references: list[str] = []
    for raw in structure.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        upper = line.upper()
        if upper in _ABACUS_SECTION_HEADERS:
            section = upper
            continue
        tokens = line.split()
        if section == "ATOMIC_SPECIES" and len(tokens) >= 3:
            references.append(tokens[2])
        elif section == "NUMERICAL_ORBITAL" and tokens:
            references.append(tokens[0])

    assets: list[Path] = []
    for reference in references:
        candidate = Path(reference).expanduser()
        choices = [candidate] if candidate.is_absolute() else [
            workdir / candidate,
            workdir / "assets" / candidate.name,
        ]
        resolved = next((path.resolve() for path in choices if path.is_file()), None)
        if resolved is not None and resolved not in assets:
            assets.append(resolved)
    return assets


def _copy_abacus_assets(assets: list[Path], target: Path) -> None:
    for source in assets:
        destination = target / source.name
        if source != destination.resolve():
            shutil.copy2(source, destination)


def _phonopy_displacements(
    structure: Path,
    *,
    dim: str,
    symm_tol: float,
    calculator: str,
) -> None:
    normalized = None
    source = structure
    if calculator == "abacus":
        from .phonopy_stru import write_phonopy_compatible_stru

        normalized = structure.with_name(f".{structure.name}.zstar-phonopy")
        write_phonopy_compatible_stru(structure, normalized)
        source = normalized
    command = [
        "phonopy",
        f"--dim={dim}",
        "-v",
        "-d",
        f"--{calculator}",
        "-c",
        str(source),
        f"--tolerance={float(symm_tol):g}",
    ]
    try:
        result = subprocess.run(command, text=True, capture_output=True)
    finally:
        if normalized is not None:
            normalized.unlink(missing_ok=True)
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"Phonopy displacement generation failed with exit code "
            f"{result.returncode}"
        )


def run_phonopy_and_process_files(
    f_stru: str = "STRU",
    symm_tol: float = 1e-3,
    dim: str = "1 1 1",
    abacus_sub: Optional[str] = "abacus_x.sh",
    vasp_sub: Optional[str] = "vasp_scf.sh",
    node: str = "s1",
) -> list[str]:
    """Generate displacement folders without requiring a copied job script.

    ``node`` is retained for CLI compatibility. Execution and scheduler
    selection are intentionally separate from directory generation.
    """

    del node
    structure = Path(f_stru).resolve()
    if not structure.is_file():
        raise FileNotFoundError(f"Structure file not found: {structure}")
    from .shared_abacus import MANIFEST, load_manifest
    if Path(MANIFEST).is_file():
        from .shared_response import read_structure
        metadata = load_manifest('.')
        if [int(n) for n in str(dim).split()] != [1, 1, 1]:
            raise ValueError('Prepare finite-q/supercell phonons in a separate directory; this directory contains a shared Gamma ensemble')
        given = read_structure(structure)
        reference = read_structure('0.no-move/STRU')
        if given.symbols != reference.symbols or not np.allclose(given.cell, reference.cell, atol=1e-8, rtol=0) or not np.allclose(given.positions, reference.positions, atol=1e-8, rtol=0):
            raise ValueError('The supplied STRU differs from the shared Gamma reference')
        print('[SHARED] Reusing existing displacements. Run zstar bec run, then zstar phonon post.')
        return [s['name'] for s in metadata['stages']]
    structure_text = structure.read_text(encoding="utf-8", errors="ignore")
    is_abacus = "ATOMIC_SPECIES" in structure_text
    if is_abacus:
        for required in ("INPUT", "KPT"):
            if not (Path.cwd() / required).is_file():
                raise FileNotFoundError(
                    f"{required} is required beside STRU to prepare ABACUS "
                    "phonon displacement folders"
                )
        input_text = (Path.cwd() / "INPUT").read_text(
            encoding="utf-8",
            errors="ignore",
        )
        force_enabled = re.search(
            r"(?im)^[ \t]*cal_force[ \t]+(?:1|true|t|yes)[ \t]*(?:#.*)?$",
            input_text,
        )
        if force_enabled is None:
            raise ValueError(
                "ABACUS phonon INPUT must enable `cal_force 1` before "
                "running `zstar ph`"
            )
    _phonopy_displacements(
        structure,
        dim=dim,
        symm_tol=symm_tol,
        calculator="abacus" if is_abacus else "vasp",
    )

    workdir = Path.cwd().resolve()
    stru_files = sorted(
        path
        for path in workdir.glob("STRU-*")
        if path.is_file() and "unitcell" not in path.name
    )
    poscar_files = sorted(
        path for path in workdir.glob("POSCAR-*") if path.is_file()
    )
    generated: list[str] = []

    if stru_files:
        assets = _referenced_abacus_assets(structure, workdir)
        for displaced in stru_files:
            number = displaced.name.split("-", 1)[1]
            target = workdir / f"disp-{number}"
            target.mkdir(parents=True, exist_ok=True)
            create_symlink(workdir / "INPUT", target / "INPUT")
            create_symlink(workdir / "KPT", target / "KPT")
            create_symlink(displaced, target / "STRU")
            _copy_abacus_assets(assets, target)
            _copy_optional_script(abacus_sub, target)
            generated.append(str(target))
    elif poscar_files:
        for required in ("INCAR", "KPOINTS", "POTCAR"):
            if not (workdir / required).is_file():
                raise FileNotFoundError(
                    f"{required} is required to prepare VASP phonon folders"
                )
        for displaced in poscar_files:
            number = displaced.name.split("-", 1)[1]
            target = workdir / f"disp-{number}"
            target.mkdir(parents=True, exist_ok=True)
            create_symlink(workdir / "INCAR", target / "INCAR")
            create_symlink(workdir / "KPOINTS", target / "KPOINTS")
            create_symlink(workdir / "POTCAR", target / "POTCAR")
            create_symlink(displaced, target / "POSCAR")
            _copy_optional_script(vasp_sub, target)
            generated.append(str(target))
    else:
        raise RuntimeError("Phonopy did not generate STRU-* or POSCAR-* files")

    print(f"Generated {len(generated)} phonon displacement folders.")
    if abacus_sub and not Path(abacus_sub).is_file() and stru_files:
        print(
            f"[INFO] Optional job script {abacus_sub!r} was not found; "
            "directories were generated without per-folder scripts."
        )
    return generated


if __name__ == "__main__":
    run_phonopy_and_process_files(f_stru="STRU", symm_tol=0.01, dim="4 4 1")
