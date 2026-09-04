"""Layered configuration for calculator executables and launchers."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shlex
import shutil
from typing import Any, Mapping

try:  # Python 3.11+
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - exercised on Python 3.9/3.10
    import tomli as tomllib


DEFAULT_EXECUTABLES = {
    "abacus": "abacus",
    "pyatb": "pyatb",
    "vasp": "vasp_std",
    "cp2k": "cp2k.psmp",
    "qe_pw": "pw.x",
    "qe_ph": "ph.x",
    "qe_dynmat": "dynmat.x",
    "phonopy": "phonopy",
}
EXECUTION_SYSTEM_ALIASES = {
    "shell": "shell",
    "local": "shell",
    "bash": "shell",
    "slurm": "slurm",
    "torque": "torque",
    "pbs": "torque",
    "openpbs": "torque",
}
SUPPORTED_EXECUTION_SYSTEMS = ("shell", "slurm", "torque")
ENV_EXECUTABLES = {
    name: f"ZSTAR_{name.upper()}_EXECUTABLE" for name in DEFAULT_EXECUTABLES
}


def normalize_execution_system(system: str) -> str:
    """Return ZStar's canonical execution-system name.

    ``shell`` and ``torque`` are the implementation names used in manifests.
    ``local``/``bash`` and ``pbs``/``openpbs`` are accepted user-facing
    aliases so the CLI matches common cluster terminology without duplicating
    scheduler implementations.
    """

    key = str(system).strip().lower()
    try:
        return EXECUTION_SYSTEM_ALIASES[key]
    except KeyError as exc:
        choices = ", ".join(sorted(EXECUTION_SYSTEM_ALIASES))
        raise ValueError(
            f"Unknown execution system {system!r}; choose one of: {choices}"
        ) from exc


def user_config_path() -> Path:
    override = os.environ.get("ZSTAR_CONFIG")
    if override:
        return Path(override).expanduser().resolve()
    if os.name == "nt" and os.environ.get("APPDATA"):
        return Path(os.environ["APPDATA"]) / "zstar" / "config.toml"
    root = Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config"))
    return root / "zstar" / "config.toml"


def project_config_path(root: str | Path = ".") -> Path:
    return Path(root).resolve() / ".zstar" / "config.toml"


def _deep_merge(base: dict[str, Any], update: Mapping[str, Any]) -> dict[str, Any]:
    for key, value in update.items():
        if isinstance(value, Mapping) and isinstance(base.get(key), dict):
            _deep_merge(base[key], value)
        else:
            base[str(key)] = value
    return base


def _read_toml(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    with path.open("rb") as handle:
        return dict(tomllib.load(handle))


def load_config(root: str | Path = ".") -> dict[str, Any]:
    """Load defaults, user config, project config, then environment overrides."""

    data: dict[str, Any] = {
        "executables": dict(DEFAULT_EXECUTABLES),
        "execution": {"system": "shell", "tasks": 1, "cpus_per_task": 1},
        "abacus": {"pseudo_dir": "", "orbital_dir": ""},
    }
    _deep_merge(data, _read_toml(user_config_path()))
    _deep_merge(data, _read_toml(project_config_path(root)))
    executables = data.setdefault("executables", {})
    for name, variable in ENV_EXECUTABLES.items():
        if os.environ.get(variable):
            executables[name] = os.environ[variable]
    return data


def resolve_executable(name: str, *, root: str | Path = ".", override: str | None = None) -> str:
    if override:
        return str(override)
    key = str(name).lower()
    if key not in DEFAULT_EXECUTABLES:
        raise KeyError(f"Unknown ZStar executable key: {name}")
    return str(load_config(root).get("executables", {}).get(key, DEFAULT_EXECUTABLES[key]))


def launcher_command(
    name: str,
    *,
    root: str | Path = ".",
    system: str = "shell",
    tasks: int | None = None,
    override: str | None = None,
) -> str:
    """Build an MPI launcher around a configured executable."""

    executable = shlex.quote(resolve_executable(name, root=root, override=override))
    tasks, _ = resolve_parallelism(root, tasks=tasks)
    if int(tasks) < 1:
        raise ValueError("tasks must be positive")
    system_key = normalize_execution_system(system)
    if system_key == "slurm":
        return f"srun --ntasks={int(tasks)} {executable}"
    if system_key in {"shell", "torque"}:
        return f"mpirun -np {int(tasks)} {executable}"
    raise AssertionError("normalize_execution_system returned an invalid system")


def resolve_parallelism(root='.', *, tasks=None, cpus_per_task=None):
    execution = load_config(root).get('execution', {})
    mpi = execution.get('mpi', execution.get('tasks', 1)) if tasks is None else tasks
    omp = execution.get('omp', execution.get('cpus_per_task', 1)) if cpus_per_task is None else cpus_per_task
    if isinstance(mpi, bool) or isinstance(omp, bool) or int(mpi) != float(mpi) or int(omp) != float(omp):
        raise ValueError('MPI ranks and OMP threads must be positive integers')
    if min(int(mpi), int(omp)) < 1:
        raise ValueError('MPI ranks and OMP threads must be positive integers')
    return int(mpi), int(omp)


def executable_available(command: str) -> tuple[bool, str | None]:
    """Return whether the first executable in a configured command is available."""

    try:
        token = shlex.split(str(command), posix=os.name != "nt")[0]
    except (ValueError, IndexError):
        return False, None
    expanded = Path(token).expanduser()
    if expanded.is_file():
        return True, str(expanded.resolve())
    found = shutil.which(token.strip('"'))
    return found is not None, found


def config_report(root: str | Path = ".") -> dict[str, Any]:
    from .job_headers import header_locations
    data = load_config(root)
    checks = {}
    for name, command in data["executables"].items():
        available, resolved = executable_available(str(command))
        checks[name] = {
            "command": str(command),
            "available": available,
            "resolved": resolved,
        }
    asset_checks = {}
    for name in ("pseudo_dir", "orbital_dir"):
        configured = str(data.get("abacus", {}).get(name, "") or "")
        if configured:
            path = Path(configured).expanduser().resolve()
            asset_checks[name] = {"path": str(path), "available": path.is_dir()}
        else:
            asset_checks[name] = {"path": "", "available": False}
    return {
        "user_config": str(user_config_path()),
        "project_config": str(project_config_path(root)),
        "configuration": data,
        "executables": checks,
        "abacus_assets": asset_checks,
        "job_headers": header_locations(root),
    }


def _toml_scalar(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    return json.dumps(str(value), ensure_ascii=True)


def write_config(path: str | Path, data: Mapping[str, Any]) -> Path:
    target = Path(path).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for section, values in data.items():
        if not isinstance(values, Mapping):
            continue
        if lines:
            lines.append("")
        lines.append(f"[{section}]")
        for key, value in values.items():
            if isinstance(value, Mapping):
                raise ValueError("Nested TOML sections deeper than one level are unsupported")
            lines.append(f"{key} = {_toml_scalar(value)}")
    target.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return target


def initialize_config(*, root: str | Path = ".", user: bool = False, force: bool = False) -> Path:
    target = user_config_path() if user else project_config_path(root)
    if target.exists() and not force:
        raise FileExistsError(target)
    return write_config(
        target,
        {
            "executables": DEFAULT_EXECUTABLES,
            "execution": {"system": "shell", "mpi": 1, "omp": 1},
            "abacus": {"pseudo_dir": "", "orbital_dir": ""},
        },
    )


def set_config_value(
    key: str,
    value: str,
    *,
    root: str | Path = ".",
    user: bool = False,
) -> Path:
    parts = str(key).split(".", 1)
    if len(parts) != 2 or not all(parts):
        raise ValueError("Configuration keys must use section.name notation")
    target = user_config_path() if user else project_config_path(root)
    data = _read_toml(target)
    section = data.setdefault(parts[0], {})
    if not isinstance(section, dict):
        raise ValueError(f"Configuration section is not a table: {parts[0]}")
    parsed: Any = value
    if str(value).lower() in {"true", "false"}:
        parsed = str(value).lower() == "true"
    else:
        try:
            parsed = int(value)
        except ValueError:
            try:
                parsed = float(value)
            except ValueError:
                pass
    section[parts[1]] = parsed
    return write_config(target, data)
