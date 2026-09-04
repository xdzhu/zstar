"""Serial, resumable execution for finite-displacement phonon folders."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
from typing import Iterable

from .configuration import launcher_command, normalize_execution_system
from .workflow import scf_is_complete


@dataclass
class PhononStageState:
    name: str
    path: str
    calculator: str
    status: str
    started_at: str | None = None
    finished_at: str | None = None
    error: str | None = None


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def discover_phonon_stages(root: str | Path = ".") -> list[tuple[Path, str]]:
    root_path = Path(root).resolve()
    stages = []
    for path in root_path.glob("disp-*"):
        if not path.is_dir():
            continue
        calculator = "abacus" if (path / "STRU").is_file() else (
            "vasp" if (path / "POSCAR").is_file() else ""
        )
        if calculator:
            stages.append((path, calculator))
    stages.sort(
        key=lambda item: (
            int(re.search(r"(\d+)$", item[0].name).group(1))
            if re.search(r"(\d+)$", item[0].name)
            else 10**12,
            item[0].name,
        )
    )
    if not stages:
        raise FileNotFoundError(f"No disp-* phonon folders found under {root_path}")
    return stages


def _vasp_complete(path: Path) -> bool:
    outcar = path / "OUTCAR"
    if not outcar.is_file() or outcar.stat().st_size == 0:
        return False
    tail = outcar.read_text(encoding="utf-8", errors="ignore")[-100000:]
    return (
        "General timing and accounting informations" in tail
        or "Voluntary context switches" in tail
    )


def phonon_stage_complete(path: Path, calculator: str) -> bool:
    return scf_is_complete(path) if calculator == "abacus" else _vasp_complete(path)


def _state_path(root: Path) -> Path:
    return root / ".zstar" / "phonon_state.json"


def _write_states(root: Path, states: Iterable[PhononStageState]) -> None:
    target = _state_path(root)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(
            {"updated_at": _utc_now(), "stages": [asdict(item) for item in states]},
            indent=2,
        ),
        encoding="utf-8",
        newline="\n",
    )


def phonon_workflow_status(root: str | Path = ".") -> list[PhononStageState]:
    root_path = Path(root).resolve()
    saved = {}
    state_path = _state_path(root_path)
    if state_path.is_file():
        saved = {
            item["name"]: item
            for item in json.loads(state_path.read_text(encoding="utf-8")).get("stages", [])
        }
    states = []
    for path, calculator in discover_phonon_stages(root_path):
        previous = saved.get(path.name, {})
        complete = phonon_stage_complete(path, calculator)
        status = "completed" if complete else previous.get("status", "pending")
        if status == "completed" and not complete:
            status = "pending"
        states.append(
            PhononStageState(
                name=path.name,
                path=str(path),
                calculator=calculator,
                status=status,
                started_at=previous.get("started_at"),
                finished_at=previous.get("finished_at"),
                error=previous.get("error"),
            )
        )
    return states


def run_phonon_workflow(
    root: str | Path = ".",
    *,
    command: str | None = None,
    omp_threads: int = 1,
    dry_run: bool = False,
    stop_after: int | None = None,
) -> list[PhononStageState]:
    root_path = Path(root).resolve()
    from .shared_abacus import MANIFEST
    if (root_path / MANIFEST).is_file():
        raise ValueError('This is a shared BEC/Gamma ensemble. Use zstar bec run to preserve reference-first charge reuse and the insulating gate, then zstar phonon post.')
    states = phonon_workflow_status(root_path)
    executed = 0
    for state in states:
        path = Path(state.path)
        if phonon_stage_complete(path, state.calculator):
            state.status = "completed"
            continue
        if stop_after is not None and executed >= int(stop_after):
            state.status = "pending"
            continue
        stage_command = command or (
            launcher_command("abacus", root=root_path)
            if state.calculator == "abacus"
            else launcher_command("vasp", root=root_path)
        )
        state.started_at = _utc_now()
        state.error = None
        if dry_run:
            state.status = "dry-run"
        else:
            state.status = "running"
            _write_states(root_path, states)
            environment = os.environ.copy()
            environment["OMP_NUM_THREADS"] = str(max(1, int(omp_threads)))
            result = subprocess.run(stage_command, cwd=path, env=environment, shell=True)
            if result.returncode != 0 or not phonon_stage_complete(path, state.calculator):
                state.status = "failed"
                state.error = f"{state.calculator} exited with code {result.returncode} or produced no complete output"
                _write_states(root_path, states)
                raise RuntimeError(f"{state.name}: {state.error}")
            state.status = "completed"
            state.finished_at = _utc_now()
        executed += 1
        _write_states(root_path, states)
    _write_states(root_path, states)
    return states


def format_phonon_status(states: Iterable[PhononStageState]) -> str:
    rows = [(item.name, item.calculator, item.status, item.error or "") for item in states]
    headers = ("stage", "calculator", "status", "error")
    widths = [max(len(headers[i]), *(len(row[i]) for row in rows)) for i in range(4)]
    lines = ["  ".join(headers[i].ljust(widths[i]) for i in range(4))]
    lines.append("  ".join("-" * width for width in widths))
    lines.extend("  ".join(row[i].ljust(widths[i]) for i in range(4)) for row in rows)
    return "\n".join(lines)


def generate_phonon_script(
    root: str | Path = ".",
    *,
    backend: str = "shell",
    output: str | Path | None = None,
    job_name: str = "zstar-phonon",
    nodes: int = 1,
    tasks: int | None = None,
    cpus_per_task: int | None = None,
    walltime: str = "24:00:00",
    queue: str | None = None,
    account: str | None = None,
    env_script: str | Path | None = None,
    header_file: str | Path | None = None,
    command: str | None = None,
    dry_run: bool = False,
) -> Path:
    from .configuration import resolve_parallelism
    from .job_headers import compose_job_script, torque_ppn
    tasks, cpus_per_task = resolve_parallelism(root, tasks=tasks, cpus_per_task=cpus_per_task)
    backend_key = normalize_execution_system(backend)
    if min(int(nodes), int(tasks), int(cpus_per_task)) < 1:
        raise ValueError("nodes, tasks, and cpus_per_task must be positive")
    root_path = Path(root).resolve()
    stages = discover_phonon_stages(root_path)
    calculators = {calculator for _, calculator in stages}
    if len(calculators) != 1:
        raise ValueError("A phonon workflow cannot mix calculator types")
    calculator = next(iter(calculators))
    if command is None:
        command = launcher_command(
            calculator, root=root_path, system=backend_key, tasks=int(tasks)
        )
    suffix = {"shell": "sh", "slurm": "slurm", "torque": "pbs"}[backend_key]
    target = Path(output).resolve() if output else root_path / f"run_zstar_phonon.{suffix}"
    header = ["#!/usr/bin/env bash", f"# ZStar phonon execution system: {backend_key}"]
    if backend_key == "slurm":
        header.extend([
            f"#SBATCH --job-name={job_name}", f"#SBATCH --nodes={int(nodes)}",
            f"#SBATCH --ntasks={int(tasks)}", f"#SBATCH --cpus-per-task={int(cpus_per_task)}",
            f"#SBATCH --time={walltime}", f"#SBATCH --output={root_path}/.zstar/slurm-%j.out",
        ])
        if queue:
            header.append(f"#SBATCH --partition={queue}")
        if account:
            header.append(f"#SBATCH --account={account}")
    elif backend_key == "torque":
        header.extend([
            f"#PBS -N {job_name}", f"#PBS -l nodes={int(nodes)}:ppn={torque_ppn(nodes, tasks, cpus_per_task)}",
            f"#PBS -l walltime={walltime}", f"#PBS -o {root_path}/.zstar/torque.out",
            f"#PBS -e {root_path}/.zstar/torque.err",
        ])
        if queue:
            header.append(f"#PBS -q {queue}")
        if account:
            header.append(f"#PBS -A {account}")
    body = [
        "set -euo pipefail", f"ROOT={shlex.quote(str(root_path))}", 'mkdir -p "$ROOT/.zstar"', 'cd "$ROOT"'
    ]
    if env_script:
        body.append(f"source {shlex.quote(str(Path(env_script).expanduser().resolve()))}")
    body.extend([
        f"export OMP_NUM_THREADS={int(cpus_per_task)}",
        f"zstar phonon run --root \"$ROOT\" --command {shlex.quote(command)} --omp-threads {int(cpus_per_task)}"
        f"{' --dry-run' if dry_run else ''} 2>&1 | tee -a \"$ROOT/.zstar/phonon.log\"",
    ])
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(compose_job_script(root_path, backend_key, header, body, specified=header_file),
                      encoding="utf-8", newline="\n")
    if os.name != "nt":
        target.chmod(target.stat().st_mode | 0o111)
    return target
