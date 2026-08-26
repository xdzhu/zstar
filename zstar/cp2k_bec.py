"""CP2K finite-displacement and native APT workflows for Born charges."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
from typing import Iterable, Optional, Sequence

import numpy as np


DEBYE_PER_E_ANGSTROM = 4.80320471257
_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
_SECTION_RE = re.compile(r"^\s*&\s*([A-Za-z][A-Za-z0-9_]*)\b", re.I)
_END_RE = re.compile(r"^\s*&\s*END(?:\s+([A-Za-z][A-Za-z0-9_]*))?\b", re.I)
_MOMENT_RE = re.compile(
    rf"X\s*=\s*({_NUMBER})\s+Y\s*=\s*({_NUMBER})\s+Z\s*=\s*({_NUMBER})",
    re.I,
)
_QUANTUM_RE = re.compile(
    rf"\[([XYZ])\]\s*=?\s*\[\s*({_NUMBER})\s+({_NUMBER})\s+({_NUMBER})\s*\]",
    re.I,
)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def _without_comment(line: str) -> str:
    return line.split("#", 1)[0].split("!", 1)[0]


@dataclass(frozen=True)
class Coordinate:
    label: str
    xyz: tuple[float, float, float]
    line_index: int
    indent: str = "    "
    suffix: str = ""


@dataclass(frozen=True)
class Cp2kMoment:
    dipole_debye: tuple[float, float, float]
    quantum_debye: tuple[tuple[float, float, float], ...]


@dataclass
class Cp2kStageState:
    name: str
    path: str
    status: str = "pending"
    started_at: Optional[str] = None
    finished_at: Optional[str] = None
    output: Optional[str] = None
    error: Optional[str] = None


def _section_spans(lines: Sequence[str]) -> list[tuple[tuple[str, ...], int, int]]:
    stack: list[tuple[str, int, tuple[str, ...]]] = []
    spans: list[tuple[tuple[str, ...], int, int]] = []
    for index, raw in enumerate(lines):
        line = _without_comment(raw)
        if _END_RE.match(line):
            if not stack:
                raise ValueError(f"Unmatched CP2K &END at line {index + 1}")
            _name, start, path = stack.pop()
            spans.append((path, start, index))
            continue
        match = _SECTION_RE.match(line)
        if not match:
            continue
        name = match.group(1).upper()
        path = (*((entry[0] for entry in stack)), name)
        stack.append((name, index, path))
    if stack:
        raise ValueError(f"Unclosed CP2K section: {'/'.join(stack[-1][2])}")
    return spans


def _find_section(lines: Sequence[str], path: Sequence[str]) -> tuple[int, int]:
    wanted = tuple(item.upper() for item in path)
    matches = [(start, end) for current, start, end in _section_spans(lines) if current == wanted]
    if not matches:
        raise ValueError(f"CP2K section not found: {'/'.join(wanted)}")
    if len(matches) > 1:
        raise ValueError(f"Repeated CP2K section is not supported: {'/'.join(wanted)}")
    return matches[0]


def _has_section(lines: Sequence[str], path: Sequence[str]) -> bool:
    wanted = tuple(item.upper() for item in path)
    return any(current == wanted for current, _start, _end in _section_spans(lines))


def validate_cp2k_bec_input(text: str) -> None:
    """Reject inputs outside the validated Berry-phase BEC domain."""

    upper = text.upper()
    if re.search(r"(?m)^\s*&\s*KPOINTS\b", upper):
        raise ValueError("CP2K BEC currently supports Gamma-point inputs only")
    if re.search(r"(?m)^\s*&\s*SMEAR\b", upper) or re.search(
        r"(?m)^\s*ADDED_MOS\s+[1-9]", upper
    ):
        raise ValueError("CP2K BEC requires an insulating, integer-occupation SCF")
    if not re.search(r"(?m)^\s*&\s*OT\b", upper):
        raise ValueError("CP2K BEC requires &SCF/&OT for periodic Berry-phase fields")
    if re.search(r"(?m)^\s*COORD_FILE_NAME\b", upper):
        raise ValueError("External CP2K coordinate files are not supported; use inline &COORD")
    if re.search(r"(?m)^\s*SCALED\s+(?:T|TRUE|\.TRUE\.)", upper):
        raise ValueError("Scaled CP2K coordinates are not supported")


def parse_cp2k_coordinates(text: str) -> tuple[list[str], list[Coordinate], float]:
    """Read inline CP2K coordinates and return input-units per Angstrom."""

    lines = text.splitlines()
    start, end = _find_section(lines, ("FORCE_EVAL", "SUBSYS", "COORD"))
    header = _without_comment(lines[start]).lower()
    units_per_angstrom = 1.0
    if "[bohr]" in header or "[a.u.]" in header:
        units_per_angstrom = 1.8897261254578281
    elif "[angstrom]" in header or "[ang]" in header:
        units_per_angstrom = 1.0

    coordinates: list[Coordinate] = []
    pattern = re.compile(
        rf"^(\s*)(\S+)\s+({_NUMBER})\s+({_NUMBER})\s+({_NUMBER})(.*)$"
    )
    for index in range(start + 1, end):
        clean = _without_comment(lines[index]).strip()
        if not clean:
            continue
        unit_match = re.match(r"(?i)^UNIT\s+(\S+)", clean)
        if unit_match:
            unit = unit_match.group(1).strip("[]").lower()
            if unit in {"bohr", "a.u.", "au"}:
                units_per_angstrom = 1.8897261254578281
            elif unit in {"angstrom", "ang"}:
                units_per_angstrom = 1.0
            else:
                raise ValueError(f"Unsupported CP2K coordinate unit: {unit}")
            continue
        match = pattern.match(lines[index])
        if not match:
            raise ValueError(f"Cannot parse CP2K coordinate at line {index + 1}: {lines[index]}")
        coordinates.append(
            Coordinate(
                label=match.group(2),
                xyz=tuple(_float(match.group(i)) for i in range(3, 6)),
                line_index=index,
                indent=match.group(1),
                suffix=match.group(6),
            )
        )
    if not coordinates:
        raise ValueError("No atoms found in CP2K &COORD")
    return lines, coordinates, units_per_angstrom


def _render_displacement(
    lines: Sequence[str], coordinate: Coordinate, direction: int, delta_input: float
) -> str:
    output = list(lines)
    xyz = list(coordinate.xyz)
    xyz[direction] += delta_input
    output[coordinate.line_index] = (
        f"{coordinate.indent}{coordinate.label:<6}"
        f" {xyz[0]: .12f} {xyz[1]: .12f} {xyz[2]: .12f}{coordinate.suffix}"
    )
    return "\n".join(output) + "\n"


def _insert_before_section_end(text: str, path: Sequence[str], block: Sequence[str]) -> str:
    lines = text.splitlines()
    _start, end = _find_section(lines, path)
    lines[end:end] = list(block)
    return "\n".join(lines) + "\n"


def _set_section_keyword(
    text: str, path: Sequence[str], keyword: str, value: str
) -> str:
    """Set a scalar keyword in a unique CP2K section."""

    lines = text.splitlines()
    start, end = _find_section(lines, path)
    pattern = re.compile(rf"(?i)^\s*{re.escape(keyword)}\b")
    for index in range(start + 1, end):
        if pattern.match(_without_comment(lines[index])):
            indent = re.match(r"^\s*", lines[index]).group(0)
            lines[index] = f"{indent}{keyword} {value}"
            return "\n".join(lines) + "\n"
    section_indent = re.match(r"^\s*", lines[start]).group(0)
    lines.insert(end, f"{section_indent}  {keyword} {value}")
    return "\n".join(lines) + "\n"


def ensure_periodic_moments(text: str) -> str:
    lines = text.splitlines()
    moments_path = ("FORCE_EVAL", "DFT", "PRINT", "MOMENTS")
    if _has_section(lines, moments_path):
        start, end = _find_section(lines, moments_path)
        has_periodic = False
        for index in range(start + 1, end):
            if re.match(r"(?i)^\s*PERIODIC\b", _without_comment(lines[index])):
                lines[index] = "        PERIODIC TRUE"
                has_periodic = True
        if not has_periodic:
            lines.insert(end, "        PERIODIC TRUE")
        return "\n".join(lines) + "\n"

    print_path = ("FORCE_EVAL", "DFT", "PRINT")
    block = ["      &MOMENTS", "        PERIODIC TRUE", "      &END MOMENTS"]
    if _has_section(lines, print_path):
        return _insert_before_section_end(text, print_path, block)
    return _insert_before_section_end(
        text,
        ("FORCE_EVAL", "DFT"),
        ["    &PRINT", *block, "    &END PRINT"],
    )


def _set_restart_guess(text: str) -> str:
    lines = text.splitlines()
    scf_path = ("FORCE_EVAL", "DFT", "SCF")
    start, end = _find_section(lines, scf_path)
    found_guess = False
    for index in range(start + 1, end):
        if re.match(r"(?i)^\s*SCF_GUESS\b", _without_comment(lines[index])):
            indent = re.match(r"^\s*", lines[index]).group(0)
            lines[index] = f"{indent}SCF_GUESS RESTART"
            found_guess = True
    if not found_guess:
        lines.insert(end, "      SCF_GUESS RESTART")
    text = "\n".join(lines) + "\n"

    lines = text.splitlines()
    dft_start, dft_end = _find_section(lines, ("FORCE_EVAL", "DFT"))
    for index in range(dft_start + 1, dft_end):
        if re.match(r"(?i)^\s*WFN_RESTART_FILE_NAME\b", _without_comment(lines[index])):
            indent = re.match(r"^\s*", lines[index]).group(0)
            lines[index] = f"{indent}WFN_RESTART_FILE_NAME reference-RESTART.wfn"
            return "\n".join(lines) + "\n"
    lines.insert(dft_end, "    WFN_RESTART_FILE_NAME reference-RESTART.wfn")
    return "\n".join(lines) + "\n"


def _parse_atom_selection(selection: str | Sequence[int] | None, natoms: int) -> list[int]:
    if selection is None or (isinstance(selection, str) and selection.strip().lower() == "all"):
        return list(range(1, natoms + 1))
    if isinstance(selection, str):
        values: list[int] = []
        for token in selection.replace(" ", "").split(","):
            if not token:
                continue
            if "-" in token:
                first, last = (int(value) for value in token.split("-", 1))
                values.extend(range(first, last + 1))
            else:
                values.append(int(token))
    else:
        values = [int(value) for value in selection]
    values = list(dict.fromkeys(values))
    if not values or min(values) < 1 or max(values) > natoms:
        raise ValueError(f"Atom selection must lie in 1..{natoms}")
    return values


def _referenced_assets(text: str, source_dir: Path) -> list[Path]:
    keywords = {
        "BASIS_SET_FILE_NAME",
        "POTENTIAL_FILE_NAME",
        "KERNEL_FILE_NAME",
        "PARAM_FILE_NAME",
    }
    assets: list[Path] = []
    for raw in text.splitlines():
        clean = _without_comment(raw).strip()
        if not clean or clean.startswith("&"):
            continue
        fields = shlex.split(clean, posix=True)
        if len(fields) < 2 or fields[0].upper() not in keywords:
            continue
        for value in fields[1:]:
            candidate = Path(os.path.expandvars(value)).expanduser()
            if candidate.is_absolute():
                continue
            resolved = source_dir / candidate
            if resolved.is_file() and resolved not in assets:
                assets.append(resolved)
    return assets


def prepare_cp2k_bec(
    input_path: str | Path,
    root: str | Path = "cp2k_bec",
    *,
    method: str = "central",
    displacement_angstrom: float = 0.01,
    atoms: str | Sequence[int] | None = "all",
    force: bool = False,
) -> Path:
    """Generate reference and finite-displacement CP2K stages."""

    source = Path(input_path).resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    method_key = method.lower()
    if method_key in {"center", "central"}:
        method_key = "central"
    elif method_key == "forward":
        method_key = "forward"
    else:
        raise ValueError("method must be central or forward")
    if displacement_angstrom <= 0:
        raise ValueError("displacement_angstrom must be positive")

    text = source.read_text(encoding="utf-8")
    validate_cp2k_bec_input(text)
    text = ensure_periodic_moments(text)
    lines, coordinates, units_per_angstrom = parse_cp2k_coordinates(text)
    selected = _parse_atom_selection(atoms, len(coordinates))
    assets = _referenced_assets(text, source.parent)

    target = Path(root).resolve()
    if target.exists() and force:
        shutil.rmtree(target)
    if target.exists() and any(target.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {target}")
    target.mkdir(parents=True, exist_ok=True)

    stages: list[dict] = []

    def write_stage(relative: Path, stage_text: str, **metadata) -> None:
        directory = target / relative
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "input.inp").write_text(stage_text, encoding="utf-8", newline="\n")
        for asset in assets:
            shutil.copy2(asset, directory / asset.name)
        stages.append({"name": relative.as_posix(), "path": relative.as_posix(), **metadata})

    write_stage(Path("reference"), text, reference=True)
    direction_names = ("x", "y", "z")
    signs = ((1, "plus"),) if method_key == "forward" else ((1, "plus"), (-1, "minus"))
    for atom_index in selected:
        coordinate = coordinates[atom_index - 1]
        safe_label = re.sub(r"[^A-Za-z0-9_-]+", "_", coordinate.label)
        for direction, direction_name in enumerate(direction_names):
            for sign, sign_name in signs:
                displaced = _render_displacement(
                    lines,
                    coordinate,
                    direction,
                    sign * displacement_angstrom * units_per_angstrom,
                )
                displaced = _set_restart_guess(displaced)
                relative = Path(f"atom-{atom_index:04d}-{safe_label}") / f"{direction_name}-{sign_name}"
                write_stage(
                    relative,
                    displaced,
                    reference=False,
                    atom_index=atom_index,
                    label=coordinate.label,
                    direction=direction_name,
                    sign=sign,
                )

    manifest = {
        "schema_version": 1,
        "backend": "cp2k",
        "created_at": _utc_now(),
        "source_input": str(source),
        "method": method_key,
        "displacement_angstrom": float(displacement_angstrom),
        "atoms": [
            {"index": index, "label": coordinates[index - 1].label}
            for index in selected
        ],
        "natoms_total": len(coordinates),
        "assets": [str(path) for path in assets],
        "stages": stages,
    }
    (target / "cp2k_bec_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8", newline="\n"
    )
    return target


def parse_cp2k_moment(path: str | Path) -> Cp2kMoment:
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    moments = _MOMENT_RE.findall(text)
    if not moments:
        raise ValueError(f"No CP2K periodic dipole moment found in {path}")
    quantum_rows: dict[str, tuple[float, float, float]] = {}
    for match in _QUANTUM_RE.findall(text):
        quantum_rows[match[0].upper()] = tuple(_float(value) for value in match[1:])
    if set(quantum_rows) != {"X", "Y", "Z"}:
        raise ValueError(f"No complete CP2K dipole quantum matrix found in {path}")
    dipole = tuple(_float(value) for value in moments[-1])
    quantum = tuple(quantum_rows[key] for key in ("X", "Y", "Z"))
    return Cp2kMoment(dipole_debye=dipole, quantum_debye=quantum)


def unwrap_dipole_delta(
    target_debye: Sequence[float],
    reference_debye: Sequence[float],
    quantum_debye: Sequence[Sequence[float]],
) -> tuple[np.ndarray, np.ndarray]:
    """Return the nearest Berry branch and its integer quantum shift."""

    delta = np.asarray(target_debye, dtype=float) - np.asarray(reference_debye, dtype=float)
    quantum = np.asarray(quantum_debye, dtype=float)
    if quantum.shape != (3, 3) or abs(np.linalg.det(quantum)) < 1e-12:
        raise ValueError("CP2K dipole quantum matrix must be nonsingular")
    coefficients = np.linalg.solve(quantum.T, delta)
    shifts = np.rint(coefficients).astype(int)
    return delta - shifts @ quantum, shifts


def _load_manifest(root: str | Path) -> tuple[Path, dict]:
    root_path = Path(root).resolve()
    path = root_path / "cp2k_bec_manifest.json"
    if not path.is_file():
        raise FileNotFoundError(path)
    return root_path, json.loads(path.read_text(encoding="utf-8"))


def cp2k_output_complete(path: str | Path) -> bool:
    output = Path(path)
    if not output.is_file() or output.stat().st_size == 0:
        return False
    tail = output.read_text(encoding="utf-8", errors="ignore")[-50000:]
    return "PROGRAM ENDED AT" in tail and "ABORT" not in tail


def _reference_wfn(reference: Path) -> Path:
    candidates = sorted(reference.glob("*-RESTART.wfn"), key=lambda path: path.stat().st_mtime)
    if not candidates:
        raise FileNotFoundError(f"No CP2K restart wavefunction found in {reference}")
    return candidates[-1]


def run_cp2k_bec(
    root: str | Path,
    *,
    cp2k_command: str = "cp2k.psmp",
    omp_threads: int = 1,
    dry_run: bool = False,
    stop_after: Optional[int] = None,
    extra_env: Optional[dict[str, str]] = None,
) -> list[Cp2kStageState]:
    """Run all CP2K BEC stages serially and resume completed stages."""

    root_path, manifest = _load_manifest(root)
    state_dir = root_path / ".zstar"
    state_dir.mkdir(exist_ok=True)
    state_path = state_dir / "cp2k_bec_state.json"
    old = {}
    if state_path.is_file():
        old = {item["name"]: item for item in json.loads(state_path.read_text())["stages"]}

    states: list[Cp2kStageState] = []
    executed = 0
    reference_dir = root_path / "reference"
    for stage in manifest["stages"]:
        directory = root_path / stage["path"]
        output = directory / "output.log"
        previous = old.get(stage["name"], {})
        state = Cp2kStageState(
            name=stage["name"],
            path=str(directory),
            status=previous.get("status", "pending"),
            started_at=previous.get("started_at"),
            finished_at=previous.get("finished_at"),
            output=str(output),
            error=previous.get("error"),
        )
        if cp2k_output_complete(output):
            parse_cp2k_moment(output)
            state.status = "completed"
            states.append(state)
            continue
        if stop_after is not None and executed >= stop_after:
            state.status = "pending"
            states.append(state)
            continue

        if not stage.get("reference"):
            source_wfn = _reference_wfn(reference_dir) if not dry_run else None
            if source_wfn is not None:
                shutil.copy2(source_wfn, directory / "reference-RESTART.wfn")

        command = (
            f"{cp2k_command} -i {shlex.quote('input.inp')} "
            f"-o {shlex.quote('output.log')}"
        )
        state.started_at = _utc_now()
        state.error = None
        if dry_run:
            state.status = "dry-run"
        else:
            state.status = "running"
            environment = os.environ.copy()
            environment["OMP_NUM_THREADS"] = str(max(1, int(omp_threads)))
            environment.setdefault("OMP_PROC_BIND", "close")
            environment.setdefault("OMP_PLACES", "cores")
            if extra_env:
                environment.update({str(key): str(value) for key, value in extra_env.items()})
            result = subprocess.run(command, cwd=directory, env=environment, shell=True)
            if result.returncode != 0 or not cp2k_output_complete(output):
                state.status = "failed"
                state.error = f"CP2K exited with code {result.returncode}"
                states.append(state)
                _write_states(state_path, states)
                raise RuntimeError(f"{state.name}: {state.error}")
            parse_cp2k_moment(output)
            state.status = "completed"
            state.finished_at = _utc_now()
        states.append(state)
        executed += 1
        _write_states(state_path, states)
    _write_states(state_path, states)
    return states


def _write_states(path: Path, states: Sequence[Cp2kStageState]) -> None:
    path.write_text(
        json.dumps({"updated_at": _utc_now(), "stages": [asdict(item) for item in states]}, indent=2),
        encoding="utf-8",
        newline="\n",
    )


def cp2k_bec_status(root: str | Path) -> list[Cp2kStageState]:
    root_path, manifest = _load_manifest(root)
    state_path = root_path / ".zstar" / "cp2k_bec_state.json"
    recorded = {}
    if state_path.is_file():
        saved = json.loads(state_path.read_text(encoding="utf-8"))
        recorded = {item["name"]: item for item in saved.get("stages", [])}
    states = []
    for stage in manifest["stages"]:
        directory = root_path / stage["path"]
        output = directory / "output.log"
        previous = recorded.get(stage["name"], {})
        complete = cp2k_output_complete(output)
        status = "completed" if complete else previous.get("status", "pending")
        if status == "completed" and not complete:
            status = "pending"
        states.append(
            Cp2kStageState(
                name=stage["name"],
                path=str(directory),
                status=status,
                started_at=previous.get("started_at"),
                finished_at=previous.get("finished_at"),
                output=str(output),
                error=previous.get("error"),
            )
        )
    return states


def format_cp2k_status(states: Iterable[Cp2kStageState]) -> str:
    rows = [(item.name, item.status, item.error or "") for item in states]
    widths = [max([len(title), *(len(row[index]) for row in rows)]) for index, title in enumerate(("stage", "status", "error"))]
    output = ["  ".join(title.ljust(widths[index]) for index, title in enumerate(("stage", "status", "error")))]
    output.append("  ".join("-" * width for width in widths))
    output.extend("  ".join(value.ljust(widths[index]) for index, value in enumerate(row)) for row in rows)
    return "\n".join(output)


def collect_cp2k_bec(
    root: str | Path,
    *,
    output: str | Path = "Z-BORN-all.out",
    json_output: str | Path = "cp2k_bec.json",
    response_output: str | Path | None = "zstar_response.json",
) -> dict:
    """Collect periodic dipole derivatives into ZStar-order BEC tensors."""

    root_path, manifest = _load_manifest(root)
    reference = parse_cp2k_moment(root_path / "reference" / "output.log")
    method = manifest["method"]
    displacement = float(manifest["displacement_angstrom"])
    tensors = []
    diagnostics = []
    for atom in manifest["atoms"]:
        tensor = np.zeros((3, 3), dtype=float)
        atom_diag = {"index": atom["index"], "label": atom["label"], "directions": {}}
        safe_label = re.sub(r"[^A-Za-z0-9_-]+", "_", atom["label"])
        atom_dir = root_path / f"atom-{atom['index']:04d}-{safe_label}"
        for direction_index, direction in enumerate(("x", "y", "z")):
            plus = parse_cp2k_moment(atom_dir / f"{direction}-plus" / "output.log")
            if method == "central":
                minus = parse_cp2k_moment(atom_dir / f"{direction}-minus" / "output.log")
                delta, shifts = unwrap_dipole_delta(
                    plus.dipole_debye, minus.dipole_debye, plus.quantum_debye
                )
                denominator = 2.0 * displacement * DEBYE_PER_E_ANGSTROM
            else:
                delta, shifts = unwrap_dipole_delta(
                    plus.dipole_debye, reference.dipole_debye, reference.quantum_debye
                )
                denominator = displacement * DEBYE_PER_E_ANGSTROM
            tensor[direction_index, :] = delta / denominator
            atom_diag["directions"][direction] = {
                "delta_dipole_debye": delta.tolist(),
                "branch_shifts": shifts.tolist(),
            }
        tensors.append(tensor)
        diagnostics.append(atom_diag)

    tensor_array = np.asarray(tensors)
    acoustic_sum = np.sum(tensor_array, axis=0)
    result = {
        "schema_version": 1,
        "backend": "cp2k",
        "method": method,
        "displacement_angstrom": displacement,
        "natoms_total": int(manifest["natoms_total"]),
        "natoms_selected": len(manifest["atoms"]),
        "sum_scope": (
            "all_atoms"
            if len(manifest["atoms"]) == int(manifest["natoms_total"])
            else "selected_atoms"
        ),
        "tensor_convention": "rows=atomic displacement/force; columns=polarization/electric field",
        "atoms": [
            {**atom, "tensor": tensor.tolist()}
            for atom, tensor in zip(manifest["atoms"], tensor_array)
        ],
        "acoustic_sum_tensor": acoustic_sum.tolist(),
        "diagnostics": diagnostics,
    }

    output_path = Path(output)
    if not output_path.is_absolute():
        output_path = root_path / output_path
    header = (
        f"{'No. Atom': <8} "
        + " ".join(f"{name:>11}" for name in ("xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"))
        + "\n"
    )
    rows = [header]
    for atom, tensor in zip(manifest["atoms"], tensor_array):
        values = " ".join(f"{value: 11.6f}" for value in tensor.reshape(9))
        rows.append(f" {atom['index']:>4} {atom['label']:<3} {values}\n")
    output_path.write_text("".join(rows), encoding="utf-8", newline="\n")

    json_path = Path(json_output)
    if not json_path.is_absolute():
        json_path = root_path / json_path
    json_path.write_text(json.dumps(result, indent=2), encoding="utf-8", newline="\n")
    response_path = None
    if response_output is not None:
        from .response_schema import response_record_from_bec_result

        response_path = Path(response_output)
        if not response_path.is_absolute():
            response_path = root_path / response_path
        response_record_from_bec_result(
            result,
            dimensionality=3,
            provenance={
                "collector": "zstar.cp2k_bec.collect_cp2k_bec",
                "source": str(root_path.resolve()),
                "legacy_result": str(json_path.resolve()),
            },
        ).write(response_path)
    result["output"] = str(output_path)
    result["json_output"] = str(json_path)
    result["response_output"] = None if response_path is None else str(response_path)
    return result


def prepare_native_apt(
    input_path: str | Path,
    output_dir: str | Path = "cp2k_native_apt",
    *,
    field_strength: float = 3.0e-4,
    force: bool = False,
) -> Path:
    """Prepare CP2K 2025.2+ native finite-field APT input."""

    source = Path(input_path).resolve()
    text = source.read_text(encoding="utf-8")
    validate_cp2k_bec_input(text)
    if field_strength <= 0:
        raise ValueError("field_strength must be positive")
    upper = text.upper()
    if "APT_FD" in upper:
        raise ValueError("Input already contains APT_FD")
    # CP2K's APT regression inputs use RUN_TYPE ENERGY. APT_FD performs the
    # six finite-field force evaluations internally.
    text = _set_section_keyword(text, ("GLOBAL",), "RUN_TYPE", "ENERGY")
    lines = text.splitlines()
    apt_print = ["        &PRINT", "          &APT ON", "          &END APT", "        &END PRINT"]
    apt_body = [
        "        APT_FD TRUE",
        f"        APT_FD_DE {field_strength:.12g}",
        *apt_print,
    ]
    dcdr_path = ("FORCE_EVAL", "PROPERTIES", "LINRES", "DCDR")
    print_path = (*dcdr_path, "PRINT")
    apt_path = (*print_path, "APT")
    if _has_section(lines, dcdr_path):
        text = _insert_before_section_end(
            text,
            dcdr_path,
            ["        APT_FD TRUE", f"        APT_FD_DE {field_strength:.12g}"],
        )
        lines = text.splitlines()
        if not _has_section(lines, apt_path):
            if _has_section(lines, print_path):
                text = _insert_before_section_end(
                    text, print_path, ["          &APT ON", "          &END APT"]
                )
            else:
                text = _insert_before_section_end(text, dcdr_path, apt_print)
    elif _has_section(lines, ("FORCE_EVAL", "PROPERTIES", "LINRES")):
        text = _insert_before_section_end(
            text,
            ("FORCE_EVAL", "PROPERTIES", "LINRES"),
            ["      &DCDR", *apt_body, "      &END DCDR"],
        )
    elif _has_section(lines, ("FORCE_EVAL", "PROPERTIES")):
        text = _insert_before_section_end(
            text,
            ("FORCE_EVAL", "PROPERTIES"),
            [
                "    &LINRES",
                "      &DCDR",
                *apt_body,
                "      &END DCDR",
                "    &END LINRES",
            ],
        )
    else:
        text = _insert_before_section_end(
            text,
            ("FORCE_EVAL",),
            [
                "  &PROPERTIES",
                "    &LINRES",
                "      &DCDR",
                *apt_body,
                "      &END DCDR",
                "    &END LINRES",
                "  &END PROPERTIES",
            ],
        )
    target = Path(output_dir).resolve()
    if target.exists() and force:
        shutil.rmtree(target)
    if target.exists() and any(target.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {target}")
    target.mkdir(parents=True, exist_ok=True)
    (target / "input.inp").write_text(text, encoding="utf-8", newline="\n")
    for asset in _referenced_assets(text, source.parent):
        shutil.copy2(asset, target / asset.name)
    return target


def run_native_apt(
    directory: str | Path,
    *,
    cp2k_command: str = "cp2k.ssmp",
    omp_threads: int = 1,
    extra_env: Optional[dict[str, str]] = None,
) -> Path:
    """Run a prepared CP2K native APT calculation and return its data file."""

    root = Path(directory).resolve()
    input_path = root / "input.inp"
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    output = root / "output.log"
    if cp2k_output_complete(output):
        existing = sorted(root.glob("*-apt*.data"))
        if existing:
            parse_native_apt(existing[-1])
            return existing[-1]

    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(max(1, int(omp_threads)))
    environment.setdefault("OMP_PROC_BIND", "close")
    environment.setdefault("OMP_PLACES", "cores")
    if extra_env:
        environment.update({str(key): str(value) for key, value in extra_env.items()})
    command = f"{cp2k_command} -i input.inp -o output.log"
    result = subprocess.run(command, cwd=root, env=environment, shell=True)
    if result.returncode != 0 or not cp2k_output_complete(output):
        raise RuntimeError(f"CP2K native APT failed with exit code {result.returncode}")
    candidates = sorted(root.glob("*-apt*.data"))
    if not candidates:
        raise FileNotFoundError(f"CP2K did not write an APT data file under {root}")
    parse_native_apt(candidates[-1])
    return candidates[-1]


def parse_native_apt(path: str | Path) -> list[dict]:
    """Parse CP2K APT_FD output into ZStar's row-force convention.

    CP2K stores rows as electric-field directions and columns as force
    directions. ZStar stores rows as displacement/force directions and columns
    as polarization/electric-field directions, so the raw matrix is transposed.
    """

    lines = Path(path).read_text(encoding="utf-8", errors="ignore").splitlines()
    atoms: list[dict] = []
    index = 0
    header = re.compile(rf"^\s*(\d+)\s+(\S+)\s+({_NUMBER})\s*$")
    while index < len(lines):
        match = header.match(lines[index])
        if not match:
            index += 1
            continue
        if index + 3 >= len(lines):
            raise ValueError(f"Truncated CP2K APT tensor in {path}")
        matrix = []
        for row in lines[index + 1 : index + 4]:
            values = re.findall(_NUMBER, row)
            if len(values) != 3:
                raise ValueError(f"Invalid CP2K APT tensor row: {row}")
            matrix.append([_float(value) for value in values])
        raw_matrix = np.asarray(matrix, dtype=float)
        atoms.append(
            {
                "index": int(match.group(1)),
                "label": match.group(2),
                "gapt": _float(match.group(3)),
                "tensor": raw_matrix.T.tolist(),
                "tensor_raw_cp2k": raw_matrix.tolist(),
            }
        )
        index += 4
    if not atoms:
        raise ValueError(f"No CP2K native APT tensors found in {path}")
    return atoms


def compare_cp2k_bec(zstar_json: str | Path, native_apt: str | Path) -> dict:
    zstar = json.loads(Path(zstar_json).read_text(encoding="utf-8"))
    native = parse_native_apt(native_apt)
    native_by_index = {int(atom["index"]): atom for atom in native}
    differences = []
    per_atom = []
    for left in zstar["atoms"]:
        atom_index = int(left["index"])
        if atom_index not in native_by_index:
            raise ValueError(f"Atom {atom_index} is missing from CP2K native APT")
        right = native_by_index[atom_index]
        delta = np.asarray(left["tensor"], dtype=float) - np.asarray(right["tensor"], dtype=float)
        differences.append(delta)
        per_atom.append(
            {
                "index": left["index"],
                "label": left["label"],
                "max_abs": float(np.max(np.abs(delta))),
                "rms": float(np.sqrt(np.mean(delta**2))),
                "difference": delta.tolist(),
            }
        )
    all_delta = np.asarray(differences)
    return {
        "max_abs": float(np.max(np.abs(all_delta))),
        "rms": float(np.sqrt(np.mean(all_delta**2))),
        "per_atom": per_atom,
    }
