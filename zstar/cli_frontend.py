"""Canonical command-family routing layered over the compatibility CLI.

The established command implementations remain the source of behavior.  This
module gives them a coherent public vocabulary while old commands continue to
work during the deprecation period.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Callable, Sequence

from .configuration import (
    config_report,
    initialize_config,
    launcher_command,
    load_config,
    resolve_executable,
    set_config_value,
)
from .project_manifest import manifest_path, read_manifest, write_manifest


LegacyRunner = Callable[[Sequence[str]], None]

ACTION_ALIASES = {
    "prepare": "pre",
    "status": "stat",
    "collect": "post",
    "deal": "post",
    "script": "job",
}

FAMILY_HELP = {
    "bec": "pre, run, job, stat, post",
    "phonon": "pre, run, job, stat, post, irrep",
    "spectra": "pre, run, job, stat, post",
    "dielectric": "static, freq, optics",
    "stru": "convert, wyckoff",
    "data": "qnep, db",
    "skill": "install, path, preflight",
    "config": "init, show, set, check",
}


def _print_family_help(family: str) -> None:
    print(f"usage: zstar {family} <action> [options]")
    print(f"actions: {FAMILY_HELP[family]}")


def _option(arguments: Sequence[str], *names: str, default=None):
    for index, token in enumerate(arguments):
        for name in names:
            if token == name and index + 1 < len(arguments):
                return arguments[index + 1]
            if token.startswith(name + "="):
                return token.split("=", 1)[1]
    return default


def _drop_options(arguments: Sequence[str], *names: str) -> list[str]:
    output: list[str] = []
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
        output.append(token)
    return output


def _replace_option(arguments: Sequence[str], old: str, new: str) -> list[str]:
    output = []
    for token in arguments:
        if token == old:
            output.append(new)
        elif token.startswith(old + "="):
            output.append(new + "=" + token.split("=", 1)[1])
        else:
            output.append(token)
    return output


def _has_option(arguments: Sequence[str], *names: str) -> bool:
    return any(
        token in names or any(token.startswith(name + "=") for name in names)
        for token in arguments
    )


def _manifest_defaults(root: str, calculator: str | None, dimensionality: int | None):
    path = manifest_path(root, "bec")
    if path.is_file():
        saved = read_manifest("bec", root)
        calculator = calculator or str(saved["calculator"])
        dimensionality = (
            int(saved["dimensionality"])
            if dimensionality is None
            else dimensionality
        )
        options = dict(saved.get("options", {}))
    else:
        options = {}
    return calculator or "abacus", 3 if dimensionality is None else dimensionality, options


def _run_bec(arguments: Sequence[str], legacy: LegacyRunner) -> None:
    if not arguments or arguments[0] in {"-h", "--help"}:
        _print_family_help("bec")
        return
    action = ACTION_ALIASES.get(arguments[0], arguments[0])
    if action not in {"pre", "run", "stat", "post", "job"}:
        raise SystemExit(f"Unknown zstar bec action: {arguments[0]}")
    rest = list(arguments[1:])
    calculator = _option(rest, "--calculator", "--calc")
    root = str(_option(rest, "--root", default="."))
    dim_text = _option(rest, "--dim", "--dimensionality")
    dimensionality = None if dim_text is None else int(dim_text)
    method = _option(rest, "--method")
    clean = _drop_options(rest, "--calculator", "--calc")

    if action == "pre":
        calculator = str(calculator or "abacus").lower()
        dimensionality = 3 if dimensionality is None else dimensionality
        if calculator == "abacus":
            clean = _drop_options(clean, "--root")
            manifest_root = str(Path(root).resolve())
            Path(manifest_root).mkdir(parents=True, exist_ok=True)
            previous = Path.cwd()
            try:
                os.chdir(manifest_root)
                legacy(["gen", *clean])
            finally:
                os.chdir(previous)
        elif calculator == "cp2k":
            legacy(["cp2k-bec", "prepare", *clean])
            manifest_root = root if root != "." else "cp2k_bec"
        elif calculator == "vasp":
            clean = _drop_options(clean, "--dim", "--dimensionality")
            legacy(["vasp-bec", "prepare", *clean])
            manifest_root = root if root != "." else "vasp_bec"
            dimensionality = 3
        elif calculator == "qe":
            legacy(["qe-bec", "prepare", *clean])
            manifest_root = root if root != "." else "qe_response"
        else:
            raise SystemExit(f"Unsupported BEC calculator: {calculator}")
        write_manifest(
            "bec",
            root=manifest_root,
            calculator=calculator,
            dimensionality=dimensionality,
            options={"method": method or ("central" if calculator == "cp2k" else "forward")},
        )
        print(f"[MANIFEST] {manifest_path(manifest_root, 'bec')}")
        return

    calculator, dimensionality, saved_options = _manifest_defaults(
        root, None if calculator is None else str(calculator).lower(), dimensionality
    )
    if not _has_option(clean, "--root"):
        clean.extend(["--root", root])

    action_map = {
        "abacus": {
            "run": ["workflow", "run"],
            "stat": ["workflow", "status"],
            "post": ["deal"],
            "job": ["workflow", "script"],
        },
        "cp2k": {
            "run": ["cp2k-bec", "run"],
            "stat": ["cp2k-bec", "status"],
            "post": ["cp2k-bec", "collect"],
            "job": ["cp2k-bec", "script"],
        },
        "vasp": {
            "run": ["vasp-bec", "run"],
            "stat": ["vasp-bec", "status"],
            "post": ["vasp-bec", "collect"],
            "job": ["vasp-bec", "script"],
        },
        "qe": {
            "run": ["qe-bec", "run"],
            "stat": ["qe-bec", "status"],
            "post": ["qe-bec", "collect"],
            "job": ["qe-bec", "script"],
        },
    }
    try:
        target = action_map[calculator][action]
    except KeyError as exc:
        raise SystemExit(f"Unsupported BEC route: {calculator} {action}") from exc
    if action == "job":
        clean = _replace_option(clean, "--system", "--backend")
    system = str(_option(clean, "--backend", default="shell"))
    tasks = int(_option(clean, "--tasks", default=1))
    if action == "run":
        if calculator == "abacus":
            if not _has_option(clean, "--abacus-command"):
                clean.extend(["--abacus-command", launcher_command("abacus", root=root)])
            if not _has_option(clean, "--pyatb-command"):
                clean.extend(["--pyatb-command", launcher_command("pyatb", root=root)])
        elif calculator == "cp2k" and not _has_option(clean, "--cp2k-command"):
            clean.extend(["--cp2k-command", resolve_executable("cp2k", root=root)])
        elif calculator == "vasp" and not _has_option(clean, "--vasp-command"):
            clean.extend(["--vasp-command", resolve_executable("vasp", root=root)])
        elif calculator == "qe":
            for flag, key in (
                ("--pw-command", "qe_pw"),
                ("--ph-command", "qe_ph"),
                ("--dynmat-command", "qe_dynmat"),
            ):
                if not _has_option(clean, flag):
                    clean.extend([flag, resolve_executable(key, root=root)])
    elif action == "job":
        if calculator == "abacus":
            if not _has_option(clean, "--abacus-command"):
                clean.extend([
                    "--abacus-command",
                    launcher_command("abacus", root=root, system=system, tasks=tasks),
                ])
            if not _has_option(clean, "--pyatb-command"):
                clean.extend([
                    "--pyatb-command",
                    launcher_command("pyatb", root=root, system=system, tasks=tasks),
                ])
        elif calculator == "cp2k" and not _has_option(clean, "--cp2k-command"):
            clean.extend([
                "--cp2k-command",
                launcher_command("cp2k", root=root, system=system, tasks=tasks),
            ])
        elif calculator == "vasp" and not _has_option(clean, "--vasp-command"):
            clean.extend([
                "--vasp-command",
                launcher_command("vasp", root=root, system=system, tasks=tasks),
            ])
        elif calculator == "qe":
            for flag, key in (
                ("--pw-command", "qe_pw"),
                ("--ph-command", "qe_ph"),
                ("--dynmat-command", "qe_dynmat"),
            ):
                if not _has_option(clean, flag):
                    clean.extend([
                        flag,
                        launcher_command(key, root=root, system=system, tasks=tasks),
                    ])
    if calculator == "abacus":
        if action in {"run", "job"} and not _has_option(clean, "--dim", "--dimensionality"):
            clean.extend(["--dimensionality", str(dimensionality)])
        if action == "post":
            clean = _drop_options(clean, "--root")
            if not _has_option(clean, "--dim", "--dimensionality"):
                clean.extend(["--dim", str(dimensionality)])
            if not _has_option(clean, "--method"):
                clean.extend(["--method", str(saved_options.get("method", "forward"))])
            previous = Path.cwd()
            try:
                os.chdir(Path(root).resolve())
                legacy([*target, *clean])
            finally:
                os.chdir(previous)
            return
    legacy([*target, *clean])


def _run_config(arguments: Sequence[str]) -> None:
    parser = argparse.ArgumentParser(prog="zstar config")
    actions = parser.add_subparsers(dest="action", required=True)
    init = actions.add_parser("init")
    init.add_argument("--root", default=".")
    init.add_argument("--user", action="store_true")
    init.add_argument("--force", action="store_true")
    show = actions.add_parser("show")
    show.add_argument("--root", default=".")
    show.add_argument("--json", action="store_true")
    set_value = actions.add_parser("set")
    set_value.add_argument("key")
    set_value.add_argument("value")
    set_value.add_argument("--root", default=".")
    set_value.add_argument("--user", action="store_true")
    check = actions.add_parser("check")
    check.add_argument("--root", default=".")
    check.add_argument("--json", action="store_true")
    args = parser.parse_args(list(arguments))
    if args.action == "init":
        print(initialize_config(root=args.root, user=args.user, force=args.force))
    elif args.action == "set":
        print(set_config_value(args.key, args.value, root=args.root, user=args.user))
    elif args.action == "show":
        report = load_config(args.root)
        print(json.dumps(report, indent=2) if args.json else _format_config(report))
    elif args.action == "check":
        report = config_report(args.root)
        if args.json:
            print(json.dumps(report, indent=2))
        else:
            for name, item in report["executables"].items():
                state = "available" if item["available"] else "missing"
                print(f"{name:<10} {state:<9} {item['command']}")


def _format_config(data: dict) -> str:
    lines = []
    for section, values in data.items():
        lines.append(f"[{section}]")
        lines.extend(f"{key} = {value}" for key, value in values.items())
        lines.append("")
    return "\n".join(lines).rstrip()


def _run_phonon(arguments: Sequence[str], legacy: LegacyRunner) -> None:
    if not arguments or arguments[0] in {"-h", "--help"}:
        _print_family_help("phonon")
        return
    action = ACTION_ALIASES.get(arguments[0], arguments[0])
    if action not in {"pre", "run", "stat", "post", "irrep", "job"}:
        raise SystemExit(f"Unknown zstar phonon action: {arguments[0]}")
    rest = list(arguments[1:])
    root = str(_option(rest, "--root", default="."))
    root_path = Path(root).resolve()
    if action in {"run", "stat", "job"}:
        from .phonon_workflow import (
            format_phonon_status,
            generate_phonon_script,
            phonon_workflow_status,
            run_phonon_workflow,
        )

        parser = argparse.ArgumentParser(prog=f"zstar phonon {action}")
        parser.add_argument('--root', default='.')
        if action == 'run':
            parser.add_argument('--command', default=None)
            parser.add_argument('--omp-threads', type=int, default=1)
            parser.add_argument('--dry-run', action='store_true')
            parser.add_argument('--stop-after', type=int, default=None)
            args = parser.parse_args(rest)
            states = run_phonon_workflow(
                args.root,
                command=args.command,
                omp_threads=args.omp_threads,
                dry_run=args.dry_run,
                stop_after=args.stop_after,
            )
            print(format_phonon_status(states))
        elif action == 'stat':
            args = parser.parse_args(rest)
            print(format_phonon_status(phonon_workflow_status(args.root)))
        else:
            parser.add_argument('--system', '--backend', choices=['shell', 'local', 'slurm', 'torque', 'pbs'], default='shell')
            parser.add_argument('--output', default=None)
            parser.add_argument('--job-name', default='zstar-phonon')
            parser.add_argument('--nodes', type=int, default=1)
            parser.add_argument('--tasks', type=int, default=1)
            parser.add_argument('--cpus-per-task', type=int, default=1)
            parser.add_argument('--walltime', default='24:00:00')
            parser.add_argument('--queue', default=None)
            parser.add_argument('--account', default=None)
            parser.add_argument('--env-script', default=None)
            parser.add_argument('--command', default=None)
            parser.add_argument('--dry-run', action='store_true')
            args = parser.parse_args(rest)
            output = generate_phonon_script(
                args.root,
                backend=args.system,
                output=args.output,
                job_name=args.job_name,
                nodes=args.nodes,
                tasks=args.tasks,
                cpus_per_task=args.cpus_per_task,
                walltime=args.walltime,
                queue=args.queue,
                account=args.account,
                env_script=args.env_script,
                command=args.command,
                dry_run=args.dry_run,
            )
            print(f"[OUT] {output}")
        return

    calculator = _option(rest, "--calculator", "--calc")
    physical_dim_text = _option(rest, "--physical-dim", default=None)
    clean = _drop_options(rest, "--root", "--calculator", "--calc")
    if action == "pre":
        physical_dim = 3 if physical_dim_text is None else int(physical_dim_text)
        clean = _drop_options(clean, "--physical-dim")
        structure = Path(_option(clean, "--stru", default="STRU"))
        structure_path = structure if structure.is_absolute() else root_path / structure
        if calculator is None and structure_path.is_file():
            text = structure_path.read_text(encoding="utf-8", errors="ignore")
            calculator = "abacus" if "ATOMIC_SPECIES" in text else "vasp"
        calculator = str(calculator or "abacus").lower()
        root_path.mkdir(parents=True, exist_ok=True)
        previous = Path.cwd()
        try:
            os.chdir(root_path)
            legacy(["ph", *clean])
        finally:
            os.chdir(previous)
        write_manifest(
            "phonon",
            root=root_path,
            calculator=calculator,
            dimensionality=physical_dim,
            options={
                "supercell": str(_option(clean, "--dim", default="1 1 1")),
                "structure": str(_option(clean, "--stru", default="STRU")),
            },
        )
        print(f"[MANIFEST] {manifest_path(root_path, 'phonon')}")
        return

    if manifest_path(root_path, "phonon").is_file():
        saved = read_manifest("phonon", root_path)
        physical_dim = int(saved["dimensionality"])
        saved_options = dict(saved.get("options", {}))
    else:
        physical_dim = 3 if physical_dim_text is None else int(physical_dim_text)
        saved_options = {}
    previous = Path.cwd()
    try:
        os.chdir(root_path)
        if action == "post":
            if not _has_option(clean, "--physical-dim"):
                clean.extend(["--physical-dim", str(physical_dim)])
            if not _has_option(clean, "--stru") and saved_options.get("structure"):
                clean.extend(["--stru", str(saved_options["structure"])])
            legacy(["postph", *clean])
        else:
            clean = _drop_options(clean, "--physical-dim")
            legacy(["irrep", *clean])
    finally:
        os.chdir(previous)


def handle_canonical_cli(arguments: Sequence[str], legacy: LegacyRunner) -> bool:
    """Handle a canonical family and return whether it consumed the command."""

    if not arguments:
        return False
    family = arguments[0]
    rest = list(arguments[1:])
    if family == "bec":
        _run_bec(rest, legacy)
        return True
    phonon_actions = set(ACTION_ALIASES) | {"pre", "run", "stat", "post", "irrep", "job"}
    if family == "phonon" or (
        family == "ph" and rest and rest[0] in phonon_actions | {"-h", "--help"}
    ):
        _run_phonon(rest, legacy)
        return True
    if family == "spectra":
        if not rest or rest[0] in {"-h", "--help"}:
            _print_family_help("spectra")
            return True
        from .spectra_frontend import run_spectra_cli

        run_spectra_cli(rest, legacy)
        return True
    if family in {"dielectric", "diel"}:
        if not rest or rest[0] in {"-h", "--help"}:
            _print_family_help("dielectric")
            return True
        action = "static" if rest[0] == "zero" else rest[0]
        mapping = {"static": "calc", "freq": "freq", "optics": "optics"}
        if action not in mapping:
            raise SystemExit(f"Unknown zstar dielectric action: {rest[0]}")
        legacy([mapping[action], *rest[1:]])
        return True
    if family == "stru":
        if not rest or rest[0] in {"-h", "--help"}:
            _print_family_help("stru")
            return True
        mapping = {"convert": "vasp", "wyckoff": "wyckoff"}
        if rest[0] not in mapping:
            raise SystemExit(f"Unknown zstar stru action: {rest[0]}")
        tail = list(rest[1:])
        if rest[0] == "convert":
            target = _option(tail, "--to", default="vasp")
            if str(target).lower() not in {"vasp", "poscar"}:
                raise SystemExit("zstar stru convert currently supports only --to vasp")
            tail = _drop_options(tail, "--to")
        legacy([mapping[rest[0]], *tail])
        return True
    if family == "data":
        if not rest or rest[0] in {"-h", "--help"}:
            _print_family_help("data")
            return True
        if rest[0] not in {"qnep", "db"}:
            raise SystemExit(f"Unknown zstar data action: {rest[0]}")
        legacy(rest)
        return True
    if family == "skill":
        if not rest or rest[0] in {"-h", "--help"}:
            _print_family_help("skill")
            return True
        legacy(["agent-skill", *rest])
        return True
    if family == "config":
        _run_config(rest)
        return True
    return False
