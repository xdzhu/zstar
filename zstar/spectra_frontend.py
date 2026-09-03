"""Calculator-aware routing for the canonical spectroscopy lifecycle."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shlex
from typing import Callable, Sequence

from .configuration import (
    launcher_command,
    normalize_execution_system,
    resolve_executable,
)
from .project_manifest import manifest_path, read_manifest, write_manifest


LegacyRunner = Callable[[Sequence[str]], None]
ALIASES = {"prepare": "pre", "status": "stat", "collect": "post", "script": "job"}


def _option(arguments: Sequence[str], *names: str, default=None):
    for index, token in enumerate(arguments):
        for name in names:
            if token == name and index + 1 < len(arguments):
                return arguments[index + 1]
            if token.startswith(name + "="):
                return token.split("=", 1)[1]
    return default


def _drop(arguments: Sequence[str], *names: str) -> list[str]:
    result = []
    skip = False
    for token in arguments:
        if skip:
            skip = False
            continue
        if token in names:
            skip = True
            continue
        if any(token.startswith(name + "=") for name in names):
            continue
        result.append(token)
    return result


def _has(arguments: Sequence[str], name: str) -> bool:
    return any(token == name or token.startswith(name + "=") for token in arguments)


def _saved(root: str, calculator: str | None, kind: str | None, dim: int | None):
    path = manifest_path(root, "spectra")
    options = {}
    if path.is_file():
        data = read_manifest("spectra", root)
        calculator = calculator or str(data["calculator"])
        dim = int(data["dimensionality"]) if dim is None else dim
        options = dict(data.get("options", {}))
        kind = kind or str(options.get("kind", "all"))
    elif (Path(root) / "spectra_manifest.json").is_file():
        data = json.loads((Path(root) / "spectra_manifest.json").read_text(encoding="utf-8"))
        calculator = calculator or str(data["calculator"])
        dim = int(data.get("dimensionality", 3)) if dim is None else dim
        kind = kind or "all"
    elif calculator is None:
        raise FileNotFoundError(
            f"No spectroscopy manifest found under {Path(root).resolve()}; "
            "run `zstar spectra pre` first"
        )
    return calculator or "abacus", kind or "all", 3 if dim is None else dim, options


def _write_driver(
    root: str,
    *,
    system: str,
    output: str | None,
    job_name: str,
    nodes: int,
    tasks: int,
    cpus_per_task: int,
    walltime: str,
    queue: str | None,
    account: str | None,
    env_script: str | None,
    dry_run: bool,
) -> Path:
    root_path = Path(root).resolve()
    key = normalize_execution_system(system)
    suffix = {"shell": "sh", "slurm": "slurm", "torque": "pbs"}[key]
    target = Path(output).resolve() if output else root_path / f"run_zstar_spectra.{suffix}"
    header = ["#!/usr/bin/env bash", f"# ZStar spectroscopy execution system: {key}"]
    if key == "slurm":
        header.extend([
            f"#SBATCH --job-name={job_name}", f"#SBATCH --nodes={nodes}",
            f"#SBATCH --ntasks={tasks}", f"#SBATCH --cpus-per-task={cpus_per_task}",
            f"#SBATCH --time={walltime}", f"#SBATCH --output={root_path}/.zstar/slurm-%j.out",
        ])
        if queue:
            header.append(f"#SBATCH --partition={queue}")
        if account:
            header.append(f"#SBATCH --account={account}")
    elif key == "torque":
        header.extend([
            f"#PBS -N {job_name}", f"#PBS -l nodes={nodes}:ppn={tasks * cpus_per_task}",
            f"#PBS -l walltime={walltime}", f"#PBS -o {root_path}/.zstar/torque.out",
            f"#PBS -e {root_path}/.zstar/torque.err",
        ])
        if queue:
            header.append(f"#PBS -q {queue}")
        if account:
            header.append(f"#PBS -A {account}")
    body = ["set -euo pipefail", f"ROOT={shlex.quote(str(root_path))}", 'cd "$ROOT"']
    if env_script:
        body.append(f"source {shlex.quote(str(Path(env_script).expanduser().resolve()))}")
    body.extend([
        f"export OMP_NUM_THREADS={cpus_per_task}",
        f"zstar spectra run --root \"$ROOT\" --omp-threads {cpus_per_task}"
        f"{' --dry-run' if dry_run else ''} 2>&1 | tee -a \"$ROOT/.zstar/spectra.log\"",
        'zstar spectra post --root "$ROOT" 2>&1 | tee -a "$ROOT/.zstar/spectra.log"',
    ])
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(header + [""] + body) + "\n", encoding="utf-8", newline="\n")
    if os.name != "nt":
        target.chmod(target.stat().st_mode | 0o111)
    return target


def run_spectra_cli(arguments: Sequence[str], legacy: LegacyRunner) -> None:
    if not arguments or arguments[0] in {"-h", "--help"}:
        print("usage: zstar spectra <action> [options]")
        print("actions: pre, job, run, stat, post")
        return
    action = ALIASES.get(arguments[0], arguments[0])
    if action not in {"pre", "run", "stat", "post", "job"}:
        raise SystemExit(f"Unknown zstar spectra action: {arguments[0]}")
    rest = list(arguments[1:])
    if any(token in {"-h", "--help"} for token in rest):
        print(f"usage: zstar spectra {action} [options]")
        print("workflow options: --calculator, --kind, --root, --dim")
        if action == "job":
            print("job options: --system, --tasks, --cpus-per-task, --walltime")
        return
    legacy_root = str(_option(rest, "--root", default="calculator_spectra"))
    if (
        arguments[0] in {"status", "collect", "script"}
        and not manifest_path(legacy_root, "spectra").is_file()
    ):
        legacy(["spectra", arguments[0], *rest])
        return
    if (
        arguments[0] == "run"
        and _has(rest, "--command")
        and not manifest_path(legacy_root, "spectra").is_file()
    ):
        legacy(["spectra", "run", *rest])
        return
    calculator = _option(rest, "--calculator", "--calc")
    kind = _option(rest, "--kind")
    root_given = _option(rest, "--root")
    dim_given = _option(rest, "--dim")
    dim = None if dim_given is None else int(dim_given)
    clean = _drop(rest, "--calculator", "--calc", "--kind")

    if action == "pre":
        calculator = str(calculator or "abacus").lower()
        kind = str(kind or "all").lower()
        if kind not in {"ir", "raman", "all"}:
            raise SystemExit("--kind must be ir, raman, or all")
        root = str(root_given or ("raman" if calculator == "abacus" else "calculator_spectra"))
        dim = 3 if dim is None else dim
        if calculator == "abacus":
            clean = _drop(clean, "--root", "--born", "--dielectric", "--dim")
            if kind in {"raman", "all"}:
                if not _has(clean, "--outdir"):
                    clean.extend(["--outdir", root])
                legacy(["raman", "prepare", *clean])
        elif calculator in {"vasp", "cp2k"}:
            if not _has(clean, "--root"):
                clean.extend(["--root", root])
            legacy(["spectra", "prepare", "--calculator", calculator, *clean])
        elif calculator == "qe":
            if not _has(clean, "--root"):
                clean.extend(["--root", root])
            legacy(["qe-bec", "prepare", *clean])
        else:
            raise SystemExit(f"Unsupported spectroscopy calculator: {calculator}")
        write_manifest(
            "spectra",
            root=root,
            calculator=calculator,
            dimensionality=dim,
            options={
                "kind": kind,
                "qpoints": str(Path(_option(rest, "--qpoints", default="qpoints.yaml")).resolve()),
                "born": str(Path(_option(rest, "--born", default="Z-BORN-symm.out")).resolve()),
                "dielectric": (
                    None
                    if _option(rest, "--dielectric") is None
                    else str(Path(_option(rest, "--dielectric")).resolve())
                ),
            },
        )
        print(f"[MANIFEST] {manifest_path(root, 'spectra')}")
        return

    root = str(root_given or "calculator_spectra")
    calculator, kind, dim, options = _saved(root, calculator, kind, dim)
    if action == "job":
        system = str(_option(rest, "--system", "--backend", default="shell"))
        output = _option(rest, "--output")
        target = _write_driver(
            root,
            system=system,
            output=output,
            job_name=str(_option(rest, "--job-name", default="zstar-spectra")),
            nodes=int(_option(rest, "--nodes", default=1)),
            tasks=int(_option(rest, "--tasks", default=1)),
            cpus_per_task=int(_option(rest, "--cpus-per-task", default=1)),
            walltime=str(_option(rest, "--walltime", default="24:00:00")),
            queue=_option(rest, "--queue"),
            account=_option(rest, "--account"),
            env_script=_option(rest, "--env-script"),
            dry_run="--dry-run" in rest,
        )
        print(f"[OUT] {target}")
        return

    if calculator in {"vasp", "cp2k"}:
        mapping = {"run": "run", "stat": "status", "post": "collect"}
        clean = _drop(clean, "--calculator", "--dim")
        if not _has(clean, "--root"):
            clean.extend(["--root", root])
        if action == "run" and not _has(clean, "--command"):
            executable = resolve_executable(calculator, root=root)
            command = executable if calculator == "vasp" else f"{executable} -i input.inp -o output.log"
            clean.extend(["--command", command])
        legacy(["spectra", mapping[action], *clean])
        return
    if calculator == "qe":
        mapping = {"run": "run", "stat": "status", "post": "collect"}
        clean = _drop(clean, "--calculator", "--kind", "--dim")
        if not _has(clean, "--root"):
            clean.extend(["--root", root])
        if action == "run":
            for flag, key in (("--pw-command", "qe_pw"), ("--ph-command", "qe_ph"), ("--dynmat-command", "qe_dynmat")):
                if not _has(clean, flag):
                    clean.extend([flag, resolve_executable(key, root=root)])
        legacy(["qe-bec", mapping[action], *clean])
        return

    # ABACUS/PYATB uses the established Raman displacement and IR postprocessors.
    if action == "run":
        if kind == "ir":
            print("IR post-processing requires no additional displaced SCF stages.")
            return
        clean = _drop(clean, "--root", "--calculator", "--kind")
        if not _has(clean, "--raman-dir"):
            clean.extend(["--raman-dir", root])
        if not _has(clean, "--dim"):
            clean.extend(["--dim", str(dim)])
        if not _has(clean, "--abacus-command"):
            clean.extend(["--abacus-command", launcher_command("abacus", root=root)])
        if not _has(clean, "--pyatb-command"):
            clean.extend(["--pyatb-command", launcher_command("pyatb", root=root)])
        legacy(["raman", "run", *clean])
    elif action == "stat":
        if kind == "ir":
            print("IR inputs are ready for post-processing.")
        else:
            legacy(["raman", "status", "--raman-dir", root])
    elif action == "post":
        qpoints = str(options.get("qpoints", "qpoints.yaml"))
        if kind in {"ir", "all"}:
            ir_args = ["ir", "--qpoints", qpoints, "--born", str(options.get("born", "Z-BORN-symm.out")), "--dim", str(dim)]
            if options.get("dielectric"):
                ir_args.extend(["--dielectric", str(options["dielectric"])])
            legacy(ir_args)
        if kind in {"raman", "all"}:
            legacy(["raman", "spectrum", "--raman-dir", root, "--qpoints", qpoints, "--dim", str(dim)])
