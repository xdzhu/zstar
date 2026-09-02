"""VASP linear-response workflows for Born effective charges."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess
from typing import Iterable, Optional

import numpy as np

from .configuration import normalize_execution_system


_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
_INCAR_KEY = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_]*)\s*=")
BEC_OUTPUT_PRECISION = 8
BEC_OUTPUT_WIDTH = 14


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


@dataclass
class VaspStageState:
    name: str
    path: str
    status: str = "pending"
    started_at: Optional[str] = None
    finished_at: Optional[str] = None
    output: Optional[str] = None
    error: Optional[str] = None
    gap_eV: Optional[float] = None


def _incar_keys(text: str) -> set[str]:
    keys: set[str] = set()
    for raw in text.splitlines():
        clean = raw.split("#", 1)[0].split("!", 1)[0]
        match = _INCAR_KEY.match(clean)
        if match:
            keys.add(match.group(1).upper())
    return keys


def _incar_value(text: str, key: str) -> str | None:
    pattern = re.compile(
        rf"^\s*{re.escape(key)}\s*=\s*([^#!;\s]+)", re.I | re.M
    )
    matches = pattern.findall(text)
    return matches[-1] if matches else None


def render_incar(
    text: str,
    *,
    updates: dict[str, str],
    remove: Iterable[str] = (),
    defaults: Optional[dict[str, str]] = None,
) -> str:
    """Update scalar INCAR tags while retaining unrelated user settings."""

    updates = {key.upper(): str(value) for key, value in updates.items()}
    remove_keys = {key.upper() for key in remove} | set(updates)
    output: list[str] = []
    for raw in text.splitlines():
        clean = raw.split("#", 1)[0].split("!", 1)[0]
        match = _INCAR_KEY.match(clean)
        if match and match.group(1).upper() in remove_keys:
            continue
        output.append(raw)
    if output and output[-1].strip():
        output.append("")
    output.append("# ZStar VASP Born-charge workflow")
    for key, value in updates.items():
        output.append(f"{key} = {value}")
    existing = _incar_keys("\n".join(output))
    for key, value in (defaults or {}).items():
        if key.upper() not in existing:
            output.append(f"{key.upper()} = {value}")
    return "\n".join(output) + "\n"


def _structure_labels(poscar: Path) -> list[str]:
    try:
        from pymatgen.core import Structure

        structure = Structure.from_file(poscar)
    except Exception as exc:
        raise ValueError(f"Cannot parse VASP structure {poscar}: {exc}") from exc
    return [site.specie.symbol for site in structure]


def prepare_vasp_bec(
    input_dir: str | Path,
    root: str | Path = "vasp_bec",
    *,
    method: str = "dfpt",
    field_strength: float = 0.001,
    force: bool = False,
) -> Path:
    """Prepare a reference-first VASP BEC workflow.

    ``dfpt`` uses LEPSILON. ``finite-field`` uses LCALCEPS and is useful for
    orbital-dependent functionals for which VASP DFPT is unavailable.
    """

    source = Path(input_dir).resolve()
    method_key = method.lower().replace("_", "-")
    if method_key not in {"dfpt", "finite-field"}:
        raise ValueError("method must be dfpt or finite-field")
    if field_strength <= 0:
        raise ValueError("field_strength must be positive")
    required = ("INCAR", "POSCAR", "KPOINTS", "POTCAR")
    missing = [name for name in required if not (source / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing VASP input files in {source}: {', '.join(missing)}")

    target = Path(root).resolve()
    if target.exists() and force:
        shutil.rmtree(target)
    if target.exists() and any(target.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {target}")
    target.mkdir(parents=True, exist_ok=True)

    original = (source / "INCAR").read_text(encoding="utf-8", errors="ignore")
    response_tags = ("LEPSILON", "LCALCEPS", "EFIELD_PEAD", "LOPTICS")
    reference_updates = {
        "NSW": "0",
        "IBRION": "-1",
        "LCHARG": ".TRUE.",
        "LWAVE": ".TRUE.",
    }
    convergence_override = None
    ediff_value = _incar_value(original, "EDIFF")
    try:
        ediff = abs(float(ediff_value.replace("D", "E").replace("d", "e"))) if ediff_value else None
    except ValueError:
        ediff = None
    if ediff is None or ediff > 1.0e-8:
        reference_updates["EDIFF"] = "1E-8"
        convergence_override = f"EDIFF={ediff_value or 'unset'} replaced by EDIFF=1E-8"
    occupation_override = None
    if method_key == "finite-field" and _incar_value(original, "ISMEAR") == "-5":
        # VASP warns that its PEAD minimization is not variational with the
        # tetrahedron method. Keep reference and response occupations aligned.
        reference_updates.update({"ISMEAR": "0", "SIGMA": "0.05"})
        occupation_override = "ISMEAR=-5 replaced by ISMEAR=0, SIGMA=0.05 for PEAD"
    reference_incar = render_incar(
        original,
        updates=reference_updates,
        remove=response_tags,
        defaults={"EDIFF": "1E-8"},
    )
    response_updates = {
        "NSW": "0",
        "IBRION": "-1",
        "ISTART": "1",
        "ICHARG": "1",
        "LREAL": ".FALSE.",
        "LRPA": ".FALSE.",
    }
    if convergence_override:
        response_updates["EDIFF"] = "1E-8"
    if method_key == "dfpt":
        response_updates["LEPSILON"] = ".TRUE."
    else:
        if occupation_override:
            response_updates.update({"ISMEAR": "0", "SIGMA": "0.05"})
        response_updates["LCALCEPS"] = ".TRUE."
        value = f"{field_strength:.10g}"
        response_updates["EFIELD_PEAD"] = f"{value} {value} {value}"
    response_incar = render_incar(
        original,
        updates=response_updates,
        remove=response_tags,
        defaults={"EDIFF": "1E-8"},
    )

    for stage, incar in (("reference", reference_incar), ("response", response_incar)):
        directory = target / stage
        directory.mkdir()
        (directory / "INCAR").write_text(incar, encoding="utf-8", newline="\n")
        for name in required[1:]:
            shutil.copy2(source / name, directory / name)

    labels = _structure_labels(source / "POSCAR")
    manifest = {
        "schema_version": 1,
        "backend": "vasp",
        "created_at": _utc_now(),
        "source_directory": str(source),
        "method": method_key,
        "field_strength_eV_per_angstrom": field_strength if method_key == "finite-field" else None,
        "occupation_override": occupation_override,
        "convergence_override": convergence_override,
        "natoms_total": len(labels),
        "labels": labels,
        "stages": [
            {"name": "reference", "path": "reference"},
            {"name": "response", "path": "response"},
        ],
    }
    (target / "vasp_bec_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8", newline="\n"
    )
    return target


def _load_manifest(root: str | Path) -> tuple[Path, dict]:
    root_path = Path(root).resolve()
    manifest_path = root_path / "vasp_bec_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)
    return root_path, json.loads(manifest_path.read_text(encoding="utf-8"))


def vasp_output_complete(path: str | Path) -> bool:
    output = Path(path)
    if not output.is_file() or output.stat().st_size == 0:
        return False
    tail = output.read_text(encoding="utf-8", errors="ignore")[-100000:]
    return (
        "General timing and accounting informations for this job" in tail
        or "Voluntary context switches" in tail
    )


def parse_vasp_gap(path: str | Path) -> float:
    """Return the fundamental gap from a completed vasprun.xml."""

    source = Path(path)
    try:
        from pymatgen.io.vasp.outputs import Vasprun

        run = Vasprun(
            source,
            parse_dos=False,
            parse_eigen=True,
            parse_projected_eigen=False,
            parse_potcar_file=False,
            exception_on_bad_xml=True,
        )
        gap, _cbm, _vbm, _direct = run.eigenvalue_band_properties
    except Exception as exc:
        raise ValueError(f"Cannot determine insulating gap from {source}: {exc}") from exc
    return float(gap)


def _write_states(path: Path, states: list[VaspStageState]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {"updated_at": _utc_now(), "stages": [asdict(stage) for stage in states]},
            indent=2,
        ),
        encoding="utf-8",
        newline="\n",
    )


def _copy_restart(reference: Path, response: Path) -> None:
    for name in ("WAVECAR", "CHGCAR"):
        source = reference / name
        if not source.is_file() or not source.stat().st_size:
            raise FileNotFoundError(f"Required VASP restart file is missing: {source}")
        shutil.copy2(source, response / name)


def run_vasp_bec(
    root: str | Path,
    *,
    vasp_command: str = "vasp_std",
    min_gap_eV: float = 0.01,
    dry_run: bool = False,
) -> list[VaspStageState]:
    """Run reference and response stages serially with an insulating gap gate."""

    root_path, manifest = _load_manifest(root)
    state_path = root_path / ".zstar" / "vasp_bec_state.json"
    old: dict[str, dict] = {}
    if state_path.is_file():
        data = json.loads(state_path.read_text(encoding="utf-8"))
        old = {item["name"]: item for item in data.get("stages", [])}
    states: list[VaspStageState] = []

    for stage_data in manifest["stages"]:
        name = stage_data["name"]
        directory = root_path / stage_data["path"]
        previous = old.get(name, {})
        state = VaspStageState(
            name=name,
            path=str(directory),
            status=previous.get("status", "pending"),
            started_at=previous.get("started_at"),
            finished_at=previous.get("finished_at"),
            output=str(directory / "vasp.log"),
            error=previous.get("error"),
            gap_eV=previous.get("gap_eV"),
        )
        if vasp_output_complete(directory / "OUTCAR"):
            state.status = "completed"
            if name == "reference":
                state.gap_eV = parse_vasp_gap(directory / "vasprun.xml")
                if state.gap_eV < min_gap_eV:
                    state.status = "rejected_metal"
                    state.error = (
                        f"Reference gap {state.gap_eV:.6g} eV is below "
                        f"the {min_gap_eV:.6g} eV threshold"
                    )
            states.append(state)
            _write_states(state_path, states)
            if state.status == "rejected_metal":
                break
            continue

        if name == "response":
            reference_state = next((item for item in states if item.name == "reference"), None)
            accepted_reference = {"completed", "dry-run"} if dry_run else {"completed"}
            if reference_state is None or reference_state.status not in accepted_reference:
                state.status = "blocked"
                state.error = "Reference SCF is not complete and insulating"
                states.append(state)
                _write_states(state_path, states)
                break
            if not dry_run:
                _copy_restart(root_path / "reference", directory)

        state.started_at = _utc_now()
        state.error = None
        if dry_run:
            state.status = "dry-run"
        else:
            state.status = "running"
            _write_states(state_path, [*states, state])
            try:
                with (directory / "vasp.log").open("w", encoding="utf-8") as log:
                    subprocess.run(
                        vasp_command,
                        cwd=directory,
                        shell=True,
                        check=True,
                        stdout=log,
                        stderr=subprocess.STDOUT,
                    )
                if not vasp_output_complete(directory / "OUTCAR"):
                    raise RuntimeError("VASP OUTCAR has no normal timing footer")
                state.status = "completed"
                state.finished_at = _utc_now()
                if name == "reference":
                    state.gap_eV = parse_vasp_gap(directory / "vasprun.xml")
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
        _write_states(state_path, states)
        if state.status not in {"completed", "dry-run"}:
            break
    return states


def vasp_bec_status(root: str | Path) -> list[VaspStageState]:
    root_path, manifest = _load_manifest(root)
    state_path = root_path / ".zstar" / "vasp_bec_state.json"
    saved: dict[str, dict] = {}
    if state_path.is_file():
        data = json.loads(state_path.read_text(encoding="utf-8"))
        saved = {item["name"]: item for item in data.get("stages", [])}
    states: list[VaspStageState] = []
    for stage in manifest["stages"]:
        previous = saved.get(stage["name"], {})
        directory = root_path / stage["path"]
        status = previous.get("status", "pending")
        if vasp_output_complete(directory / "OUTCAR") and status not in {"rejected_metal"}:
            status = "completed"
        states.append(
            VaspStageState(
                name=stage["name"],
                path=str(directory),
                status=status,
                started_at=previous.get("started_at"),
                finished_at=previous.get("finished_at"),
                output=str(directory / "vasp.log"),
                error=previous.get("error"),
                gap_eV=previous.get("gap_eV"),
            )
        )
    return states


def format_vasp_status(states: Iterable[VaspStageState]) -> str:
    rows = [
        (stage.name, stage.status, "" if stage.gap_eV is None else f"{stage.gap_eV:.6f}", stage.error or "")
        for stage in states
    ]
    headers = ("stage", "status", "gap_eV", "error")
    widths = [max(len(headers[i]), *(len(row[i]) for row in rows)) for i in range(len(headers))]
    lines = ["  ".join(headers[i].ljust(widths[i]) for i in range(len(headers)))]
    lines.append("  ".join("-" * width for width in widths))
    lines.extend("  ".join(row[i].ljust(widths[i]) for i in range(len(headers))) for row in rows)
    return "\n".join(lines)


def generate_vasp_backend_script(
    root: str | Path,
    *,
    backend: str = "shell",
    output: str | Path | None = None,
    job_name: str = "zstar-vasp-bec",
    nodes: int = 1,
    tasks: int = 1,
    cpus_per_task: int = 1,
    walltime: str = "24:00:00",
    queue: str | None = None,
    account: str | None = None,
    env_script: str | Path | None = None,
    vasp_command: str | None = None,
    min_gap_eV: float = 0.01,
) -> Path:
    """Generate one local, Slurm, or Torque driver for the serial workflow."""

    backend_key = normalize_execution_system(backend)
    if min(nodes, tasks, cpus_per_task) < 1:
        raise ValueError("nodes, tasks, and cpus_per_task must be positive")
    root_path, _manifest = _load_manifest(root)
    if output is None:
        suffix = {"shell": "sh", "slurm": "slurm", "torque": "pbs"}[backend_key]
        target = root_path / f"run_vasp_bec.{suffix}"
    else:
        target = Path(output).resolve()
    if vasp_command is None:
        vasp_command = (
            f"srun --ntasks={tasks} vasp_std"
            if backend_key == "slurm"
            else f"mpirun -np {tasks} vasp_std"
        )

    header = [
        "#!/usr/bin/env bash",
        f"# ZStar VASP BEC execution backend: {backend_key}",
        "# One reference-first driver with gap gating and resumable stage state.",
    ]
    if backend_key == "slurm":
        header.extend(
            [
                f"#SBATCH --job-name={job_name}",
                f"#SBATCH --nodes={nodes}",
                f"#SBATCH --ntasks={tasks}",
                f"#SBATCH --cpus-per-task={cpus_per_task}",
                f"#SBATCH --time={walltime}",
                f"#SBATCH --output={root_path}/.zstar/slurm-%j.out",
            ]
        )
        if queue:
            header.append(f"#SBATCH --partition={queue}")
        if account:
            header.append(f"#SBATCH --account={account}")
    elif backend_key == "torque":
        header.extend(
            [
                f"#PBS -N {job_name}",
                f"#PBS -l nodes={nodes}:ppn={tasks * cpus_per_task}",
                f"#PBS -l walltime={walltime}",
                f"#PBS -o {root_path}/.zstar/torque.out",
                f"#PBS -e {root_path}/.zstar/torque.err",
            ]
        )
        if queue:
            header.append(f"#PBS -q {queue}")
        if account:
            header.append(f"#PBS -A {account}")

    body = [
        "set -euo pipefail",
        f"ROOT={shlex.quote(str(root_path))}",
        'mkdir -p "$ROOT/.zstar"',
        'cd "$ROOT"',
    ]
    if env_script:
        body.append(f"source {shlex.quote(str(Path(env_script).expanduser().resolve()))}")
    body.extend(
        [
            f"export OMP_NUM_THREADS={cpus_per_task}",
            f"zstar vasp-bec run --root \"$ROOT\" --vasp-command {shlex.quote(vasp_command)} "
            f"--min-gap {min_gap_eV:.10g} 2>&1 | tee -a \"$ROOT/.zstar/workflow.log\"",
        ]
    )
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(header + [""] + body) + "\n", encoding="utf-8", newline="\n")
    if target.exists() and hasattr(target, "chmod"):
        target.chmod(target.stat().st_mode | 0o111)
    return target


def _numeric_triplet(line: str) -> list[float] | None:
    fields = line.split()
    values: list[float] = []
    for field in fields:
        if re.fullmatch(_NUMBER, field):
            values.append(_float(field))
    return values[-3:] if len(values) >= 3 else None


def parse_vasp_outcar(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Parse epsilon-infinity and BECs, normalized to ZStar tensor order."""

    source = Path(path)
    lines = source.read_text(encoding="utf-8", errors="ignore").splitlines()
    dielectric_candidates: list[np.ndarray] = []
    born_candidates: list[np.ndarray] = []
    for index, line in enumerate(lines):
        upper = line.upper()
        if (
            "MACROSCOPIC STATIC DIELECTRIC TENSOR" in upper
            and "EXCLUDING" not in upper
            and "IONIC CONTRIBUTION" not in upper
        ):
            rows: list[list[float]] = []
            for raw in lines[index + 1:index + 12]:
                values = _numeric_triplet(raw)
                if values is not None:
                    rows.append(values)
                    if len(rows) == 3:
                        dielectric_candidates.append(np.asarray(rows, dtype=float))
                        break
        if "BORN EFFECTIVE CHARGES" in upper and "EXCLUDING" not in upper:
            tensors: list[np.ndarray] = []
            cursor = index + 1
            while cursor < len(lines):
                if cursor > index + 2 and "BORN EFFECTIVE CHARGES" in lines[cursor].upper():
                    break
                ion_match = re.match(r"\s*ion\s+\d+", lines[cursor], re.I)
                if ion_match:
                    rows = []
                    for raw in lines[cursor + 1:cursor + 8]:
                        fields = raw.split()
                        if len(fields) >= 4 and fields[0] in {"1", "2", "3"}:
                            rows.append([_float(value) for value in fields[1:4]])
                            if len(rows) == 3:
                                break
                    if len(rows) == 3:
                        # VASP: rows=electric-field/polarization, columns=force/displacement.
                        # ZStar: rows=displacement/force, columns=polarization/electric-field.
                        tensors.append(np.asarray(rows, dtype=float).T)
                if tensors and re.match(r"\s*total drift", lines[cursor], re.I):
                    break
                cursor += 1
            if tensors:
                born_candidates.append(np.asarray(tensors))
    if not dielectric_candidates:
        raise ValueError(f"No dielectric tensor including local fields found in {source}")
    if not born_candidates:
        raise ValueError(f"No Born effective charges including local fields found in {source}")
    return dielectric_candidates[-1], born_candidates[-1]


def collect_vasp_bec(
    root: str | Path,
    *,
    output: str | Path = "Z-BORN-all.out",
    born_output: str | Path = "BORN",
    json_output: str | Path = "vasp_bec.json",
    response_output: str | Path | None = "zstar_response.json",
) -> dict:
    root_path, manifest = _load_manifest(root)
    epsilon, tensors = parse_vasp_outcar(root_path / "response" / "OUTCAR")
    labels = list(manifest["labels"])
    if len(tensors) != len(labels):
        raise ValueError(f"OUTCAR has {len(tensors)} BEC tensors but POSCAR has {len(labels)} atoms")
    acoustic_sum = np.sum(tensors, axis=0)
    result = {
        "schema_version": 1,
        "backend": "vasp",
        "method": manifest["method"],
        "natoms_total": len(labels),
        "sum_scope": "all_atoms",
        "tensor_convention": "rows=atomic displacement/force; columns=polarization/electric field",
        "source_tensor_convention": "VASP OUTCAR rows=electric field; columns=force",
        "source_transform": "transpose_each_tensor",
        "epsilon_infinity": epsilon.tolist(),
        "atoms": [
            {"index": index, "label": label, "tensor": tensor.tolist()}
            for index, (label, tensor) in enumerate(zip(labels, tensors), start=1)
        ],
        "acoustic_sum_tensor": acoustic_sum.tolist(),
    }

    output_path = Path(output)
    if not output_path.is_absolute():
        output_path = root_path / output_path
    header = f"{'No. Atom': <8} " + " ".join(
        f"{name:>{BEC_OUTPUT_WIDTH}}"
        for name in ("xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz")
    ) + "\n"
    rows = [header]
    for index, (label, tensor) in enumerate(zip(labels, tensors), start=1):
        values = " ".join(
            f"{value: {BEC_OUTPUT_WIDTH}.{BEC_OUTPUT_PRECISION}f}"
            for value in tensor.reshape(9)
        )
        rows.append(f" {index:>4} {label:<3} {values}\n")
    output_path.write_text("".join(rows), encoding="utf-8", newline="\n")

    born_path = Path(born_output)
    if not born_path.is_absolute():
        born_path = root_path / born_path
    born_rows = [" ".join(f"{value:.10f}" for value in epsilon.reshape(9))]
    born_rows.extend(" ".join(f"{value:.10f}" for value in tensor.reshape(9)) for tensor in tensors)
    born_path.write_text("\n".join(born_rows) + "\n", encoding="utf-8", newline="\n")

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
                "collector": "zstar.vasp_bec.collect_vasp_bec",
                "source": str((root_path / "response" / "OUTCAR").resolve()),
                "legacy_result": str(json_path.resolve()),
            },
        ).write(response_path)
    result.update(
        output=str(output_path),
        born_output=str(born_path),
        json_output=str(json_path),
        response_output=None if response_path is None else str(response_path),
    )
    return result


def compare_vasp_bec(first: str | Path, second: str | Path) -> dict:
    """Compare two normalized VASP BEC JSON files."""

    first_path = Path(first).resolve()
    second_path = Path(second).resolve()
    first_data = json.loads(first_path.read_text(encoding="utf-8"))
    second_data = json.loads(second_path.read_text(encoding="utf-8"))
    first_labels = [str(atom.get("label", "")) for atom in first_data.get("atoms", [])]
    second_labels = [str(atom.get("label", "")) for atom in second_data.get("atoms", [])]
    if first_labels != second_labels:
        raise ValueError(
            f"Atom order differs: first={first_labels}, second={second_labels}"
        )
    first_tensors = np.asarray([atom["tensor"] for atom in first_data["atoms"]], dtype=float)
    second_tensors = np.asarray([atom["tensor"] for atom in second_data["atoms"]], dtype=float)
    if first_tensors.shape != second_tensors.shape:
        raise ValueError(
            f"BEC shapes differ: first={first_tensors.shape}, second={second_tensors.shape}"
        )
    bec_delta = first_tensors - second_tensors
    first_epsilon = np.asarray(first_data["epsilon_infinity"], dtype=float)
    second_epsilon = np.asarray(second_data["epsilon_infinity"], dtype=float)
    epsilon_delta = first_epsilon - second_epsilon
    return {
        "schema_version": 1,
        "first": str(first_path),
        "second": str(second_path),
        "labels": first_labels,
        "bec_max_abs_e": float(np.max(np.abs(bec_delta))),
        "bec_rms_e": float(np.sqrt(np.mean(bec_delta**2))),
        "epsilon_max_abs": float(np.max(np.abs(epsilon_delta))),
        "epsilon_rms": float(np.sqrt(np.mean(epsilon_delta**2))),
        "bec_delta": bec_delta.tolist(),
        "epsilon_delta": epsilon_delta.tolist(),
    }
