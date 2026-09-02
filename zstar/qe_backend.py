"""Quantum ESPRESSO Gamma-point dielectric and spectroscopy workflow."""

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
from typing import Iterable, Mapping, Sequence

import numpy as np

from .dimensions import dimension_spec
from .configuration import normalize_execution_system
from .response_schema import ResponseQuantity, ResponseRecord
from .spectra import (
    calculate_native_line_spectrum,
    write_native_line_spectrum_outputs,
)


_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


def _float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass
class QeStageState:
    name: str
    path: str
    status: str = "pending"
    started_at: str | None = None
    finished_at: str | None = None
    output: str | None = None
    error: str | None = None
    gap_eV: float | None = None


def _namelist_bounds(lines: Sequence[str], section: str) -> tuple[int, int]:
    start_pattern = re.compile(rf"^\s*&\s*{re.escape(section)}\b", re.I)
    for start, line in enumerate(lines):
        if not start_pattern.search(line):
            continue
        for end in range(start + 1, len(lines)):
            if re.match(r"^\s*/\s*(?:!.*)?$", lines[end]):
                return start, end
        raise ValueError(f"Unterminated Quantum ESPRESSO &{section} namelist")
    raise ValueError(f"Quantum ESPRESSO input has no &{section} namelist")


def qe_namelist_values(text: str, section: str) -> dict[str, str]:
    lines = text.splitlines()
    start, end = _namelist_bounds(lines, section)
    values: dict[str, str] = {}
    pattern = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_]*)\s*=\s*(.*?)(?:,\s*)?(?:!.*)?$")
    for line in lines[start + 1 : end]:
        match = pattern.match(line)
        if match:
            values[match.group(1).lower()] = match.group(2).strip().rstrip(",")
    return values


def render_qe_namelist(text: str, section: str, updates: Mapping[str, str]) -> str:
    """Update scalar assignments in one QE namelist while retaining the input body."""

    had_newline = text.endswith("\n")
    lines = text.splitlines()
    start, end = _namelist_bounds(lines, section)
    remaining = {str(key).lower(): str(value) for key, value in updates.items()}
    key_pattern = re.compile(r"^(\s*)([A-Za-z][A-Za-z0-9_]*)\s*=", re.I)
    for index in range(start + 1, end):
        match = key_pattern.match(lines[index])
        if not match:
            continue
        key = match.group(2).lower()
        if key in remaining:
            lines[index] = f"{match.group(1)}{key} = {remaining.pop(key)},"
    indent = "  "
    additions = [f"{indent}{key} = {value}," for key, value in remaining.items()]
    lines[end:end] = additions
    output = "\n".join(lines)
    return output + "\n" if had_newline else output


def _unquote(value: str | None, default: str) -> str:
    if value is None:
        return default
    return value.strip().strip("'\"")


def _pseudo_files(text: str, source_dir: Path) -> list[Path]:
    control = qe_namelist_values(text, "control")
    pseudo_dir = Path(_unquote(control.get("pseudo_dir"), "."))
    if not pseudo_dir.is_absolute():
        pseudo_dir = source_dir / pseudo_dir
    lines = text.splitlines()
    try:
        start = next(
            index
            for index, line in enumerate(lines)
            if re.match(r"^\s*ATOMIC_SPECIES\b", line, re.I)
        )
    except StopIteration:
        raise ValueError("Quantum ESPRESSO input has no ATOMIC_SPECIES card") from None
    files: list[Path] = []
    for line in lines[start + 1 :]:
        fields = line.split()
        if len(fields) < 3 or not re.fullmatch(_NUMBER, fields[1]):
            break
        candidate = pseudo_dir / fields[2]
        if not candidate.is_file():
            candidate = source_dir / fields[2]
        if not candidate.is_file():
            raise FileNotFoundError(f"Referenced QE pseudopotential is missing: {fields[2]}")
        files.append(candidate.resolve())
    return files


def _disable_gamma_only_storage(text: str) -> str:
    """Use an explicit 1x1x1 grid so old ph.x can read pw.x restart data."""

    pattern = re.compile(r"^\s*K_POINTS\s+(?:\{|\()?gamma(?:\}|\))?\s*$", re.I | re.M)
    return pattern.sub("K_POINTS automatic\n1 1 1 0 0 0", text)


def prepare_qe_response(
    input_path: str | Path,
    root: str | Path = "qe_response",
    *,
    dimensionality: int = 3,
    periodic_axes: str | Iterable[str] | None = None,
    tr2_ph: float = 1.0e-14,
    raman: bool = True,
    force: bool = False,
) -> Path:
    """Prepare a serial pw.x -> ph.x -> dynmat.x response workflow."""

    dim = dimension_spec(dimensionality, periodic_axes)
    if dim.value not in {0, 3}:
        raise ValueError(
            "The initial QE response workflow supports dim=0 and dim=3; "
            "dim=1/2 intrinsic normalization is handled by the low-dimensional adapter"
        )
    if tr2_ph <= 0:
        raise ValueError("tr2_ph must be positive")
    source = Path(input_path).resolve()
    text = source.read_text(encoding="utf-8")
    control = qe_namelist_values(text, "control")
    prefix = _unquote(control.get("prefix"), "pwscf")
    text = render_qe_namelist(
        text,
        "control",
        {
            "calculation": "'scf'",
            "prefix": repr(prefix),
            "outdir": "'../scratch'",
            "pseudo_dir": "'../pseudo'",
            "restart_mode": "'from_scratch'",
        },
    )
    text = _disable_gamma_only_storage(text)
    target = Path(root).resolve()
    if target.exists() and force:
        shutil.rmtree(target)
    if target.exists() and any(target.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {target}")
    for directory in ("scf", "phonon", "dynmat", "scratch", "pseudo"):
        (target / directory).mkdir(parents=True, exist_ok=True)
    (target / "scf" / "pw.in").write_text(text, encoding="utf-8", newline="\n")
    for pseudo in _pseudo_files(source.read_text(encoding="utf-8"), source.parent):
        shutil.copy2(pseudo, target / "pseudo" / pseudo.name)

    ph_input = (
        "ZStar Gamma-point response\n"
        "&inputph\n"
        f"  tr2_ph = {tr2_ph:.10g},\n"
        f"  prefix = {prefix!r},\n"
        "  outdir = '../scratch',\n"
        "  fildyn = '../dynamical_matrix',\n"
        "  epsil = .true.,\n"
        "  trans = .true.,\n"
        "  recover = .true.,\n"
        f"  lraman = {'.true.' if raman else '.false.'},\n"
        "/\n"
        "0.0 0.0 0.0\n"
    )
    (target / "phonon" / "ph.in").write_text(ph_input, encoding="utf-8", newline="\n")
    asr = "zero-dim" if dim.value == 0 else "crystal"
    dynmat_input = (
        "&input\n"
        "  fildyn = '../dynamical_matrix',\n"
        f"  asr = '{asr}',\n"
        "/\n"
    )
    (target / "dynmat" / "dynmat.in").write_text(
        dynmat_input, encoding="utf-8", newline="\n"
    )
    manifest = {
        "schema_version": 1,
        "calculator": "qe",
        "created_at": _utc_now(),
        "source_input": str(source),
        "prefix": prefix,
        "dimensionality": dim.to_dict(),
        "raman_requested": bool(raman),
        "stages": [
            {"name": "scf", "path": "scf", "input": "pw.in", "output": "pw.out"},
            {"name": "phonon", "path": "phonon", "input": "ph.in", "output": "ph.out"},
            {"name": "dynmat", "path": "dynmat", "input": "dynmat.in", "output": "dynmat.out"},
        ],
    }
    (target / "qe_response_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8", newline="\n"
    )
    return target


def qe_output_complete(path: str | Path) -> bool:
    source = Path(path)
    if not source.is_file() or not source.stat().st_size:
        return False
    tail = source.read_text(encoding="utf-8", errors="ignore")[-100000:]
    return "JOB DONE." in tail and "Error in routine" not in tail


def parse_qe_gap(path: str | Path) -> float:
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    matches = re.findall(
        rf"highest occupied, lowest unoccupied level \(ev\):\s*({_NUMBER})\s+({_NUMBER})",
        text,
        re.I,
    )
    if not matches:
        raise ValueError(
            "QE did not print both occupied and unoccupied levels; set nbnd above the "
            "occupied count so ZStar can apply the insulating gap gate"
        )
    occupied, unoccupied = (_float(value) for value in matches[-1])
    return float(unoccupied - occupied)


def _load_manifest(root: str | Path) -> tuple[Path, dict]:
    root_path = Path(root).resolve()
    manifest_path = root_path / "qe_response_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)
    return root_path, json.loads(manifest_path.read_text(encoding="utf-8"))


def _state_path(root: Path) -> Path:
    return root / ".zstar" / "qe_response_state.json"


def _write_states(root: Path, states: Sequence[QeStageState]) -> None:
    path = _state_path(root)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {"updated_at": _utc_now(), "stages": [asdict(state) for state in states]},
            indent=2,
        ),
        encoding="utf-8",
        newline="\n",
    )


def run_qe_response(
    root: str | Path,
    *,
    pw_command: str = "pw.x",
    ph_command: str = "ph.x",
    dynmat_command: str = "dynmat.x",
    min_gap_eV: float = 0.01,
    omp_threads: int = 1,
    dry_run: bool = False,
) -> list[QeStageState]:
    """Run the prepared stages serially with gap gating and restart state."""

    root_path, manifest = _load_manifest(root)
    commands = {"scf": pw_command, "phonon": ph_command, "dynmat": dynmat_command}
    old: dict[str, dict] = {}
    if _state_path(root_path).is_file():
        saved = json.loads(_state_path(root_path).read_text(encoding="utf-8"))
        old = {item["name"]: item for item in saved.get("stages", [])}
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(int(omp_threads))
    states: list[QeStageState] = []
    for stage_data in manifest["stages"]:
        name = stage_data["name"]
        directory = root_path / stage_data["path"]
        output = directory / stage_data["output"]
        previous = old.get(name, {})
        state = QeStageState(
            name=name,
            path=str(directory),
            status=previous.get("status", "pending"),
            started_at=previous.get("started_at"),
            finished_at=previous.get("finished_at"),
            output=str(output),
            error=previous.get("error"),
            gap_eV=previous.get("gap_eV"),
        )
        if qe_output_complete(output):
            state.status = "completed"
            if name == "scf":
                try:
                    state.gap_eV = parse_qe_gap(output)
                except ValueError as exc:
                    state.status = "blocked_gap_unknown"
                    state.error = str(exc)
                if state.gap_eV is not None and state.gap_eV < min_gap_eV:
                    state.status = "rejected_metal"
                    state.error = (
                        f"Reference gap {state.gap_eV:.6g} eV is below "
                        f"the {min_gap_eV:.6g} eV threshold"
                    )
            states.append(state)
            _write_states(root_path, states)
            if state.status != "completed":
                break
            continue
        if name != "scf" and (not states or states[-1].status not in {"completed", "dry-run"}):
            state.status = "blocked"
            state.error = "Previous QE stage is not complete and accepted"
            states.append(state)
            _write_states(root_path, states)
            break
        state.started_at = _utc_now()
        state.error = None
        if dry_run:
            state.status = "dry-run"
        else:
            state.status = "running"
            _write_states(root_path, [*states, state])
            try:
                command = f"{commands[name]} -in {shlex.quote(stage_data['input'])}"
                with output.open("w", encoding="utf-8") as handle:
                    subprocess.run(
                        command,
                        cwd=directory,
                        shell=True,
                        check=True,
                        stdout=handle,
                        stderr=subprocess.STDOUT,
                        env=environment,
                    )
                if not qe_output_complete(output):
                    raise RuntimeError(f"QE {name} output has no JOB DONE footer")
                state.status = "completed"
                state.finished_at = _utc_now()
                if name == "scf":
                    state.gap_eV = parse_qe_gap(output)
                    if state.gap_eV < min_gap_eV:
                        state.status = "rejected_metal"
                        state.error = (
                            f"Reference gap {state.gap_eV:.6g} eV is below "
                            f"the {min_gap_eV:.6g} eV threshold"
                        )
            except Exception as exc:
                state.status = "failed"
                state.finished_at = _utc_now()
                state.error = str(exc)
        states.append(state)
        _write_states(root_path, states)
        if state.status not in {"completed", "dry-run"}:
            break
    return states


def qe_response_status(root: str | Path) -> list[QeStageState]:
    root_path, manifest = _load_manifest(root)
    saved: dict[str, dict] = {}
    if _state_path(root_path).is_file():
        data = json.loads(_state_path(root_path).read_text(encoding="utf-8"))
        saved = {item["name"]: item for item in data.get("stages", [])}
    states: list[QeStageState] = []
    for stage in manifest["stages"]:
        previous = saved.get(stage["name"], {})
        directory = root_path / stage["path"]
        output = directory / stage["output"]
        status = previous.get("status", "pending")
        if qe_output_complete(output) and status not in {"rejected_metal", "blocked_gap_unknown"}:
            status = "completed"
        states.append(
            QeStageState(
                name=stage["name"],
                path=str(directory),
                status=status,
                started_at=previous.get("started_at"),
                finished_at=previous.get("finished_at"),
                output=str(output),
                error=previous.get("error"),
                gap_eV=previous.get("gap_eV"),
            )
        )
    return states


def format_qe_status(states: Iterable[QeStageState]) -> str:
    rows = [
        (item.name, item.status, "" if item.gap_eV is None else f"{item.gap_eV:.6f}", item.error or "")
        for item in states
    ]
    headers = ("stage", "status", "gap_eV", "error")
    widths = [max(len(headers[i]), *(len(row[i]) for row in rows)) for i in range(4)]
    lines = ["  ".join(headers[i].ljust(widths[i]) for i in range(4))]
    lines.append("  ".join("-" * widths[i] for i in range(4)))
    lines.extend("  ".join(row[i].ljust(widths[i]) for i in range(4)) for row in rows)
    return "\n".join(lines)


def _triplet(line: str) -> list[float] | None:
    values = re.findall(_NUMBER, line)
    return [_float(value) for value in values[-3:]] if len(values) >= 3 else None


def parse_qe_ph_output(path: str | Path) -> dict:
    """Parse the final complete dielectric, BEC, polarizability, and frequency blocks."""

    source = Path(path)
    lines = source.read_text(encoding="utf-8", errors="ignore").splitlines()
    epsilon_candidates: list[np.ndarray] = []
    polar_candidates: list[np.ndarray] = []
    born_candidates: list[tuple[list[str], np.ndarray]] = []
    frequencies: dict[int, float] = {}
    for index, line in enumerate(lines):
        upper = line.upper()
        if "DIELECTRIC CONSTANT IN CARTESIAN AXIS" in upper:
            rows = [_triplet(raw) for raw in lines[index + 1 : index + 8]]
            rows = [row for row in rows if row is not None][:3]
            if len(rows) == 3:
                epsilon_candidates.append(np.asarray(rows))
        elif "POLARIZABILITY (A.U.)" in upper and "POLARIZABILITY (A^3)" in upper:
            rows: list[list[float]] = []
            for raw in lines[index + 1 : index + 6]:
                values = re.findall(_NUMBER, raw)
                if len(values) >= 6:
                    rows.append([_float(value) for value in values[-3:]])
                if len(rows) == 3:
                    break
            if len(rows) == 3:
                polar_candidates.append(np.asarray(rows))
        elif "EFFECTIVE CHARGES (D FORCE / DE)" in upper:
            labels: list[str] = []
            tensors: list[np.ndarray] = []
            cursor = index + 1
            while cursor < len(lines):
                atom = re.match(r"\s*atom\s+\d+\s+(\S+)", lines[cursor], re.I)
                if atom:
                    rows: list[list[float]] = []
                    for raw in lines[cursor + 1 : cursor + 7]:
                        if re.match(r"\s*E[xyz]\s*\(", raw, re.I):
                            row = _triplet(raw)
                            if row is not None:
                                rows.append(row)
                                if len(rows) == 3:
                                    break
                    if len(rows) == 3:
                        labels.append(atom.group(1))
                        # QE rows are electric field and columns are force.
                        tensors.append(np.asarray(rows).T)
                        cursor += 3
                elif tensors and ("DIELECTRIC" in lines[cursor].upper() or "DYNAMICAL" in lines[cursor].upper()):
                    break
                cursor += 1
            if tensors:
                born_candidates.append((labels, np.asarray(tensors)))
        match = re.search(
            rf"freq\s*\(\s*(\d+)\s*\)\s*=.*?=\s*({_NUMBER})\s*\[cm-1\]",
            line,
            re.I,
        )
        if match:
            frequencies[int(match.group(1))] = _float(match.group(2))
    if not epsilon_candidates:
        raise ValueError(f"No QE dielectric tensor found in {source}")
    if not born_candidates:
        raise ValueError(f"No QE Born effective charges found in {source}")
    labels, tensors = born_candidates[-1]
    return {
        "epsilon_infinity": epsilon_candidates[-1],
        "polarizability_A3": None if not polar_candidates else polar_candidates[-1],
        "labels": labels,
        "born_tensors": tensors,
        "frequencies_cm-1": np.asarray([frequencies[key] for key in sorted(frequencies)]),
    }


def parse_qe_dynmat_output(path: str | Path) -> dict:
    source = Path(path)
    lines = source.read_text(encoding="utf-8", errors="ignore").splitlines()
    header_index = None
    columns: list[str] = []
    for index, line in enumerate(lines):
        if re.match(r"\s*#\s*mode\s+\[cm-1\]", line, re.I):
            header_index = index
            columns = line.lstrip(" #").split()
    if header_index is None:
        raise ValueError(f"No QE dynmat mode table found in {source}")
    rows: list[list[float]] = []
    for line in lines[header_index + 1 :]:
        values = re.findall(_NUMBER, line)
        if len(values) < 4:
            if rows:
                break
            continue
        rows.append([_float(value) for value in values])
    if not rows:
        raise ValueError(f"QE dynmat mode table is empty in {source}")
    has_raman = any(column.lower() == "raman" for column in columns)
    return {
        "mode_numbers": [int(row[0]) for row in rows],
        "frequencies_cm-1": [row[1] for row in rows],
        "frequencies_thz": [row[2] for row in rows],
        "ir_activities": [row[3] for row in rows],
        "raman_activities": [row[4] for row in rows] if has_raman and all(len(row) >= 5 for row in rows) else None,
        "depolarization_ratios": [row[5] for row in rows] if has_raman and all(len(row) >= 6 for row in rows) else None,
    }


def collect_qe_response(
    root: str | Path,
    *,
    broadening_cm1: float = 8.0,
    points: int = 2001,
    plot: bool = True,
) -> dict:
    root_path, manifest = _load_manifest(root)
    ph = parse_qe_ph_output(root_path / "phonon" / "ph.out")
    modes = parse_qe_dynmat_output(root_path / "dynmat" / "dynmat.out")
    dim_data = manifest["dimensionality"]
    dim = dimension_spec(dim_data["value"], dim_data["periodic_axes"])
    quantities = [
        ResponseQuantity(
            name=(
                "atomic_polar_tensor"
                if dim.value == 0
                else "born_effective_charge"
            ),
            values=ph["born_tensors"],
            unit="e",
            normalization="per_atom",
            axes=("atom", "displacement", "polarization"),
            convention="rows=atomic displacement/force; columns=polarization/electric field",
            source="QE ph.x",
            metadata={"labels": ph["labels"], "source_transform": "transpose_each_tensor"},
        ),
        ResponseQuantity(
            name="electronic_dielectric" if dim.value == 3 else "supercell_electronic_dielectric",
            values=ph["epsilon_infinity"],
            unit="1",
            normalization="cell_volume",
            axes=("field", "polarization"),
            source="QE ph.x",
        ),
        ResponseQuantity(
            name="gamma_frequency",
            values=modes["frequencies_cm-1"],
            unit="cm^-1",
            normalization="mode",
            axes=("mode",),
            source="QE dynmat.x",
        ),
        ResponseQuantity(
            name="ir_activity",
            values=modes["ir_activities"],
            unit="(D/angstrom)^2/amu",
            normalization="mode",
            axes=("mode",),
            source="QE dynmat.x",
        ),
    ]
    if ph["polarizability_A3"] is not None:
        quantities.append(
            ResponseQuantity(
                name="molecular_polarizability",
                values=ph["polarizability_A3"],
                unit="angstrom^3",
                normalization="molecule",
                axes=("field", "dipole"),
                source="QE ph.x",
            )
        )
    if modes["raman_activities"] is not None:
        quantities.append(
            ResponseQuantity(
                name="raman_activity",
                values=modes["raman_activities"],
                unit="angstrom^4/amu",
                normalization="mode",
                axes=("mode",),
                source="QE dynmat.x",
            )
        )
    record = ResponseRecord(
        backend="qe",
        dimensionality=dim,
        quantities=tuple(quantities),
        provenance={
            "collector": "zstar.qe_backend.collect_qe_response",
            "ph_output": str((root_path / "phonon" / "ph.out").resolve()),
            "dynmat_output": str((root_path / "dynmat" / "dynmat.out").resolve()),
        },
        metadata={"prefix": manifest["prefix"], "raman_requested": manifest["raman_requested"]},
    )
    response_path = record.write(root_path / "zstar_response.json")
    common = {
        "mode_numbers": modes["mode_numbers"],
        "broadening_cm1": broadening_cm1,
        "points": points,
    }
    ir = calculate_native_line_spectrum(
        modes["frequencies_cm-1"],
        modes["ir_activities"],
        activity_kind="IR_activity",
        activity_unit="(D/Angstrom)^2/amu",
        **common,
    )
    ir_summary = write_native_line_spectrum_outputs(
        root_path / "ir_spectrum", ir, stem="ir", plot=plot
    )
    raman_summary = None
    if modes["raman_activities"] is not None:
        raman = calculate_native_line_spectrum(
            modes["frequencies_cm-1"],
            modes["raman_activities"],
            activity_kind="Raman_activity",
            activity_unit="Angstrom^4/amu",
            **common,
        )
        raman_summary = write_native_line_spectrum_outputs(
            root_path / "raman_spectrum", raman, stem="raman", plot=plot
        )
    result = {
        "schema_version": 1,
        "calculator": "qe",
        "dimensionality": dim.to_dict(),
        "mode_numbers": modes["mode_numbers"],
        "frequencies_cm-1": modes["frequencies_cm-1"],
        "epsilon_infinity": ph["epsilon_infinity"].tolist(),
        "polarizability_A3": None if ph["polarizability_A3"] is None else ph["polarizability_A3"].tolist(),
        "born_tensors": ph["born_tensors"].tolist(),
        "ir_activities": modes["ir_activities"],
        "raman_activities": modes["raman_activities"],
        "ir_summary": ir_summary,
        "raman_summary": raman_summary,
        "response_output": str(response_path),
    }
    (root_path / "qe_response.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8", newline="\n"
    )
    return result


def generate_qe_backend_script(
    root: str | Path,
    *,
    backend: str = "shell",
    output: str | Path | None = None,
    job_name: str = "zstar-qe-response",
    nodes: int = 1,
    tasks: int = 1,
    cpus_per_task: int = 1,
    walltime: str = "24:00:00",
    queue: str | None = None,
    account: str | None = None,
    env_script: str | Path | None = None,
    pw_command: str | None = None,
    ph_command: str | None = None,
    dynmat_command: str | None = None,
) -> Path:
    backend_key = normalize_execution_system(backend)
    if min(nodes, tasks, cpus_per_task) < 1:
        raise ValueError("nodes, tasks, and cpus_per_task must be positive")
    root_path, _manifest = _load_manifest(root)
    if output is None:
        suffix = {"shell": "sh", "slurm": "slurm", "torque": "pbs"}[backend_key]
        target = root_path / f"run_qe_response.{suffix}"
    else:
        target = Path(output).resolve()
    launcher = f"srun --ntasks={tasks}" if backend_key == "slurm" else f"mpirun -np {tasks}"
    pw_command = pw_command or f"{launcher} pw.x"
    ph_command = ph_command or f"{launcher} ph.x"
    dynmat_command = dynmat_command or f"{launcher} dynmat.x"
    header = ["#!/usr/bin/env bash", f"# ZStar QE response backend: {backend_key}"]
    if backend_key == "slurm":
        header.extend([
            f"#SBATCH --job-name={job_name}", f"#SBATCH --nodes={nodes}",
            f"#SBATCH --ntasks={tasks}", f"#SBATCH --cpus-per-task={cpus_per_task}",
            f"#SBATCH --time={walltime}", f"#SBATCH --output={root_path}/.zstar/slurm-%j.out",
        ])
        if queue:
            header.append(f"#SBATCH --partition={queue}")
        if account:
            header.append(f"#SBATCH --account={account}")
    elif backend_key == "torque":
        header.extend([
            f"#PBS -N {job_name}", f"#PBS -l nodes={nodes}:ppn={tasks * cpus_per_task}",
            f"#PBS -l walltime={walltime}", f"#PBS -o {root_path}/.zstar/torque.out",
            f"#PBS -e {root_path}/.zstar/torque.err",
        ])
        if queue:
            header.append(f"#PBS -q {queue}")
        if account:
            header.append(f"#PBS -A {account}")
    body = [
        "set -euo pipefail", f"ROOT={shlex.quote(str(root_path))}",
        'mkdir -p "$ROOT/.zstar"', 'cd "$ROOT"',
    ]
    if env_script:
        body.append(f"source {shlex.quote(str(Path(env_script).expanduser().resolve()))}")
    body.extend([
        f"export OMP_NUM_THREADS={cpus_per_task}",
        f"zstar qe-bec run --root \"$ROOT\" --pw-command {shlex.quote(pw_command)} "
        f"--ph-command {shlex.quote(ph_command)} --dynmat-command {shlex.quote(dynmat_command)} "
        f"--omp-threads {cpus_per_task} 2>&1 | tee -a \"$ROOT/.zstar/workflow.log\"",
        'zstar qe-bec collect --root "$ROOT" 2>&1 | tee -a "$ROOT/.zstar/workflow.log"',
    ])
    target.write_text("\n".join(header + [""] + body) + "\n", encoding="utf-8", newline="\n")
    target.chmod(target.stat().st_mode | 0o111)
    return target
