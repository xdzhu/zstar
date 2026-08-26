"""Installation and machine-readable preflight support for the ZStar Agent Skill."""

from __future__ import annotations

from collections import Counter
import json
import os
from pathlib import Path
import shutil
import sys
from typing import Any

from . import __version__


SKILL_NAME = "run-zstar-workflows"
LANES = ("bec", "phonon", "ir", "raman", "dielectric", "md", "cp2k", "database")
DIMENSIONS = ("molecule", "1d", "2d", "bulk")


def packaged_skill_path() -> Path:
    path = Path(__file__).resolve().parent / "agent_skills" / SKILL_NAME
    if not (path / "SKILL.md").is_file():
        raise FileNotFoundError(f"Packaged Agent Skill is incomplete: {path}")
    return path


def default_codex_skill_root() -> Path:
    codex_home = os.environ.get("CODEX_HOME")
    return (Path(codex_home).expanduser() if codex_home else Path.home() / ".codex") / "skills"


def install_agent_skill(
    destination_root: str | Path | None = None,
    *,
    force: bool = False,
) -> Path:
    """Install the packaged skill beneath a Codex-compatible skills directory."""

    root = Path(destination_root).expanduser() if destination_root else default_codex_skill_root()
    destination = root / SKILL_NAME
    source = packaged_skill_path()
    root.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if not force:
            raise FileExistsError(
                f"Skill already exists: {destination}. Re-run with --force to replace it."
            )
    temporary = root / f".{SKILL_NAME}.tmp"
    if temporary.exists():
        shutil.rmtree(temporary)
    shutil.copytree(
        source,
        temporary,
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )
    if destination.exists():
        shutil.rmtree(destination)
    try:
        temporary.replace(destination)
    except PermissionError:
        # Directory replacement can be denied transiently on Windows even
        # after a successful removal. A fresh copy preserves force semantics.
        shutil.copytree(temporary, destination)
        shutil.rmtree(temporary)
    return destination.resolve()


def _command_record(name: str) -> dict[str, Any]:
    resolved = shutil.which(name)
    return {"available": resolved is not None, "path": resolved}


def _state_summary(root: Path) -> dict[str, Any]:
    stage_dir = root / ".zstar" / "stages"
    counts: Counter[str] = Counter()
    errors: list[dict[str, str]] = []
    if stage_dir.is_dir():
        for path in sorted(stage_dir.glob("*.json")):
            try:
                data = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                errors.append({"file": path.name, "error": str(exc)})
                continue
            status = str(data.get("status", "unknown"))
            counts[status] += 1
            if data.get("error"):
                errors.append({"stage": str(data.get("name", path.stem)), "error": str(data["error"])})
    return {
        "present": stage_dir.is_dir(),
        "counts": dict(sorted(counts.items())),
        "errors": errors,
        "event_log": str(root / ".zstar" / "workflow.jsonl"),
    }


def _artifact_record(root: Path, relative: str) -> dict[str, Any]:
    path = root / relative
    return {
        "path": relative,
        "exists": path.exists(),
        "bytes": path.stat().st_size if path.is_file() else None,
    }


def preflight_report(
    root: str | Path = ".",
    *,
    lane: str = "bec",
    dimensionality: str = "bulk",
) -> dict[str, Any]:
    """Return a non-mutating readiness and provenance report for an agent."""

    if lane not in LANES:
        raise ValueError(f"Unsupported lane {lane!r}; choose from {', '.join(LANES)}")
    if dimensionality not in DIMENSIONS:
        raise ValueError(
            f"Unsupported dimensionality {dimensionality!r}; choose from {', '.join(DIMENSIONS)}"
        )
    root_path = Path(root).expanduser().resolve()
    blockers: list[str] = []
    warnings: list[str] = []
    if not root_path.is_dir():
        blockers.append(f"Workspace does not exist: {root_path}")

    command_names = {
        "bec": ("abacus", "pyatb_input", "pyatb"),
        "phonon": ("abacus", "phonopy"),
        "ir": ("phonopy",),
        "raman": ("abacus", "pyatb_input", "pyatb"),
        "dielectric": ("phonopy",),
        "md": (),
        "cp2k": ("cp2k.psmp",),
        "database": (),
    }[lane]
    commands = {name: _command_record(name) for name in command_names}
    missing_commands = [name for name, record in commands.items() if not record["available"]]
    if missing_commands:
        warnings.append(
            "Executables were not found on PATH: "
            + ", ".join(missing_commands)
            + ". Absolute commands or environment activation may still satisfy the workflow."
        )

    artifacts = {
        relative: _artifact_record(root_path, relative)
        for relative in (
            "STRU", "0.no-move", "BORN", "Z-BORN-symm.out",
            "phonopy_disp.yaml", "FORCE_SETS", "phonopy.yaml",
            "qpoints.yaml", "irreps.yaml",
            "ir_spectrum/ir_summary.json",
            "raman_spectrum/raman_summary.json",
            "dielectric_response/ir_summary.json",
            "md_dielectric/md_dielectric_summary.json",
        )
    }

    if lane in {"bec", "phonon"} and not artifacts["STRU"]["exists"]:
        blockers.append("STRU is required in the selected workspace.")
    if lane in {"ir", "raman", "dielectric"} and not artifacts["qpoints.yaml"]["exists"]:
        blockers.append("qpoints.yaml is required for the selected mode-resolved workflow.")
    if lane in {"ir", "dielectric"} and not (
        artifacts["BORN"]["exists"] or artifacts["Z-BORN-symm.out"]["exists"]
    ):
        blockers.append("A BORN or Z-BORN-symm.out tensor file is required.")
    if lane == "raman" and not artifacts["0.no-move"]["exists"]:
        blockers.append("A completed or prepared 0.no-move reference directory is required.")
    if lane == "cp2k" and root_path.is_dir() and not any(root_path.glob("*.inp")):
        warnings.append("No CP2K *.inp template was detected in the workspace root.")
    if lane == "database" and not (root_path / "candidates.csv").is_file():
        warnings.append("No candidates.csv manifest was detected; create one with zstar db init.")
    if lane == "md":
        warnings.append(
            "MD readiness also requires an explicit trajectory and either fixed or frame-resolved BEC tensors."
        )
    if dimensionality == "2d":
        warnings.append(
            "Confirm that the slab normal is Cartesian z; out-of-plane BEC requires charge-density cubes."
        )
    if dimensionality == "1d":
        warnings.append(
            "Use the current z-periodic orthogonal-wire convention: transverse BEC "
            "from neutral charge-density cubes and longitudinal BEC from the PYATB "
            "Berry phase. Do not apply bulk NAC; finite-q polar dispersion still "
            "requires a genuine 1D Coulomb cutoff."
        )
    if dimensionality == "molecule" and lane not in {"ir", "raman"}:
        warnings.append("Molecule mode is currently intended for the IR/Raman workflow.")

    return {
        "schema_version": 1,
        "skill": SKILL_NAME,
        "zstar_version": __version__,
        "python": sys.version.split()[0],
        "workspace": str(root_path),
        "lane": lane,
        "dimensionality": dimensionality,
        "ready": not blockers,
        "blockers": blockers,
        "warnings": warnings,
        "commands": commands,
        "state": _state_summary(root_path),
        "artifacts": artifacts,
    }


def write_preflight_json(
    root: str | Path = ".",
    *,
    lane: str = "bec",
    dimensionality: str = "bulk",
) -> str:
    return json.dumps(
        preflight_report(root, lane=lane, dimensionality=dimensionality),
        indent=2,
        sort_keys=True,
    )
