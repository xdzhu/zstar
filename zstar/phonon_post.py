"""Collect finite-displacement forces and generate Gamma-point phonon data."""

from __future__ import annotations

from pathlib import Path
import subprocess
from typing import Optional

import yaml

from .phonopy_stru import write_phonopy_compatible_stru


def _dim_text(value) -> str:
    if isinstance(value, (list, tuple)):
        return " ".join(str(int(item)) for item in value)
    return str(value).strip()


def get_phonopy_params(
    yaml_file: str | Path = "phonopy_disp.yaml",
) -> tuple[str, Optional[float], Optional[str]]:
    """Read supercell dimensions and symmetry metadata from Phonopy YAML."""

    path = Path(yaml_file)
    if not path.is_file():
        raise FileNotFoundError(f"Phonopy displacement metadata not found: {path}")
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    configuration = (data.get("phonopy") or {}).get("configuration") or {}
    dim = configuration.get("dim")
    if dim is None:
        matrix = (data.get("supercell_matrix") or [])
        if len(matrix) == 3 and all(len(row) == 3 for row in matrix):
            dim = [matrix[index][index] for index in range(3)]
    if dim is None:
        raise ValueError(f"No supercell dimension was found in {path}")
    tolerance = configuration.get("symmetry_tolerance")
    space_group = (data.get("space_group") or {}).get("type")
    result = (
        _dim_text(dim),
        float(tolerance) if tolerance is not None else None,
        str(space_group) if space_group is not None else None,
    )
    print(
        "Detected Phonopy settings: "
        f"DIM={result[0]}, tolerance={result[1]}, space_group={result[2]}"
    )
    return result


def _run_phonopy(command: list[str], *, stage: str) -> None:
    print("$ " + " ".join(command))
    result = subprocess.run(command, text=True, capture_output=True)
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"Phonopy {stage} failed with exit code {result.returncode}"
        )


def run_eigen_irrep(
    f_stru: str = "STRU",
    symm_tol: float = 1e-3,
    nac: bool = False,
    dim: Optional[str] = None,
) -> dict:
    """Collect forces and write Gamma eigenvectors plus irreducible modes."""

    dim_auto, tolerance_auto, space_group = get_phonopy_params()
    dimension = _dim_text(dim) if dim is not None else dim_auto
    tolerance = (
        float(symm_tol)
        if symm_tol is not None
        else (tolerance_auto if tolerance_auto is not None else 1.0e-3)
    )
    structure = Path(f_stru)
    if not structure.is_file():
        raise FileNotFoundError(f"Structure file not found: {structure}")

    force_logs = sorted(Path(".").glob("disp-*/OUT*/running*.log"))
    if not force_logs:
        raise FileNotFoundError(
            "No disp-*/OUT*/running*.log files were found for force collection"
        )
    _run_phonopy(
        ["phonopy", "-f", *(str(path) for path in force_logs)],
        stage="force collection",
    )

    normalized = structure.with_name(f".{structure.name}.zstar-phonopy")
    write_phonopy_compatible_stru(structure, normalized)
    try:
        common = [
            f"--dim={dimension}",
            "-v",
            "-c",
            str(normalized),
            f"--tolerance={tolerance:g}",
            "--abacus",
            "--qpoints=0 0 0",
        ]
        qpoint_command = ["phonopy", *common, "--eigenvectors"]
        if nac:
            qpoint_command.extend(["--nac", "--nac_method=GONZE"])
        _run_phonopy(qpoint_command, stage="Gamma eigenvector calculation")

        irrep_command = [
            "phonopy",
            f"--dim={dimension}",
            "-c",
            str(normalized),
            f"--tolerance={tolerance:g}",
            "--abacus",
            "--pa=auto",
            "--qpoints=0 0 0",
            "--irreps=0 0 0 1e-3",
        ]
        if nac:
            irrep_command.extend(["--nac", "--nac_method=GONZE"])
        _run_phonopy(
            irrep_command,
            stage="irreducible-representation calculation",
        )
    finally:
        normalized.unlink(missing_ok=True)

    outputs = {
        "dimension": dimension,
        "symmetry_tolerance": tolerance,
        "space_group": space_group,
        "force_logs": len(force_logs),
        "qpoints": str(Path("qpoints.yaml").resolve()),
        "irreps": str(Path("irreps.yaml").resolve()),
        "nac": bool(nac),
    }
    print(f"Collected {len(force_logs)} force calculations.")
    return outputs


if __name__ == "__main__":
    run_eigen_irrep()
