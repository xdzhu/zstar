"""Serial, resumable ABACUS/PYATB workflow execution for Born charges."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
from typing import Iterable, Optional, Sequence

from .pyatb_compat import (
    DEFAULT_LEGACY_DOMEGA_EV,
    DEFAULT_LEGACY_OMEGA_MAX_EV,
    configure_optical_input,
    configure_polarization_input,
    detect_pyatb_capabilities,
    read_band_gap,
)
from .configuration import normalize_execution_system


_ATOM_DIR_RE = re.compile(r"^(\d+)\.([^.]+)$")
_DIRECTION_RE = re.compile(r"^([xyz])([+-]?)$")
_DIRECTION_ORDER = {
    "x+": 0,
    "x-": 1,
    "x": 2,
    "y+": 3,
    "y-": 4,
    "y": 5,
    "z+": 6,
    "z-": 7,
    "z": 8,
}


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass(frozen=True)
class WorkflowStage:
    name: str
    path: Path
    reference: bool = False


@dataclass
class StageState:
    name: str
    path: str
    status: str = "pending"
    scf: str = "pending"
    band: str = "pending"
    band_gap_eV: Optional[float] = None
    pyatb: str = "pending"
    dielectric: str = "not-requested"
    started_at: Optional[str] = None
    finished_at: Optional[str] = None
    error: Optional[str] = None


def discover_stages(root: str | Path) -> list[WorkflowStage]:
    """Return reference first, followed by atom/direction stages."""

    root_path = Path(root).resolve()
    reference = root_path / "0.no-move"
    if not reference.is_dir():
        raise FileNotFoundError(f"Reference directory not found: {reference}")

    from .shared_abacus import MANIFEST, load_manifest
    if (root_path / MANIFEST).is_file():
        metadata = load_manifest(root_path)
        return [WorkflowStage("0.no-move", reference, reference=True),
                *(WorkflowStage(item["name"], root_path / item["name"])
                  for item in metadata["stages"])]

    displaced: list[tuple[int, str, int, WorkflowStage]] = []
    for atom_dir in root_path.iterdir():
        if not atom_dir.is_dir():
            continue
        atom_match = _ATOM_DIR_RE.match(atom_dir.name)
        if not atom_match:
            continue
        atom_number = int(atom_match.group(1))
        for direction_dir in atom_dir.iterdir():
            if not direction_dir.is_dir():
                continue
            direction_match = _DIRECTION_RE.match(direction_dir.name)
            if not direction_match:
                continue
            direction_key = direction_dir.name
            order = _DIRECTION_ORDER.get(direction_key, 99)
            stage = WorkflowStage(
                name=f"{atom_dir.name}/{direction_dir.name}",
                path=direction_dir.resolve(),
            )
            displaced.append((atom_number, atom_dir.name, order, stage))

    displaced.sort(key=lambda item: (item[0], item[1], item[2]))
    return [
        WorkflowStage(name="0.no-move", path=reference.resolve(), reference=True),
        *(item[3] for item in displaced),
    ]


class WorkflowStateStore:
    def __init__(self, root: str | Path):
        self.root = Path(root).resolve()
        self.state_dir = self.root / ".zstar"
        self.stage_dir = self.state_dir / "stages"
        self.log_path = self.state_dir / "workflow.jsonl"
        self.stage_dir.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def _safe_name(name: str) -> str:
        return name.replace("/", "__").replace("\\", "__")

    def stage_path(self, name: str) -> Path:
        return self.stage_dir / f"{self._safe_name(name)}.json"

    def load(self, stage: WorkflowStage) -> StageState:
        path = self.stage_path(stage.name)
        if not path.is_file():
            return StageState(name=stage.name, path=str(stage.path))
        data = json.loads(path.read_text(encoding="utf-8"))
        return StageState(**data)

    def save(self, state: StageState) -> None:
        path = self.stage_path(state.name)
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(asdict(state), indent=2), encoding="utf-8")
        temporary.replace(path)

    def event(self, stage: str, event: str, **details) -> None:
        record = {"time": _utc_now(), "stage": stage, "event": event, **details}
        with self.log_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, sort_keys=True) + "\n")


def _find_scf_log(stage_dir: Path) -> Optional[Path]:
    candidates = sorted(stage_dir.glob("OUT.*/running_scf.log"))
    return candidates[-1] if candidates else None


def _contains_completion_marker(path: Path) -> bool:
    try:
        with path.open("rb") as handle:
            handle.seek(max(0, path.stat().st_size - 20000))
            tail = handle.read().decode("utf-8", errors="ignore")
            normalized_tail = " ".join(tail.lower().split())
            if "total time" not in normalized_tail:
                return False

            handle.seek(0)
            for raw_line in handle:
                normalized_line = b" ".join(raw_line.lower().split())
                if b"charge density convergence is achieved" in normalized_line:
                    return True
    except OSError:
        return False
    return False


def scf_is_complete(stage_dir: str | Path) -> bool:
    log = _find_scf_log(Path(stage_dir))
    return bool(log and log.stat().st_size > 0 and _contains_completion_marker(log))


def polarization_is_complete(
    stage_dir: str | Path, subdir: str = "pyatb"
) -> bool:
    path = Path(stage_dir) / subdir / "Out" / "Polarization" / "polarization.dat"
    return path.is_file() and path.stat().st_size > 0


def dielectric_is_complete(stage_dir: str | Path) -> bool:
    optical = Path(stage_dir) / "pyatb" / "Out" / "Optical_Conductivity"
    candidates = (
        optical / "static_dielectric_function.dat",
        optical / "dielectric_function_real_part.dat",
    )
    return any(path.is_file() and path.stat().st_size > 0 for path in candidates)


def band_is_complete(stage_dir: str | Path) -> bool:
    try:
        read_band_gap(Path(stage_dir) / "pyatb-band")
        return True
    except (FileNotFoundError, ValueError):
        return False


def _reference_charge_files(reference_dir: Path) -> list[Path]:
    output_dirs = sorted(reference_dir.glob("OUT.*"))
    files: list[Path] = []
    for output_dir in output_dirs:
        for pattern in ("SPIN*_CHG.cube", "*-CHARGE-DENSITY.restart"):
            files.extend(
                path
                for path in output_dir.glob(pattern)
                if path.is_file() and path.stat().st_size > 0
            )
    return sorted(files)


def _abacus_output_dir(stage_dir: Path) -> tuple[Path, str]:
    input_path = stage_dir / "INPUT"
    if not input_path.is_file():
        input_path = stage_dir / "INPUT-scf"
    suffix = "ABACUS"
    if input_path.is_file():
        match = re.search(
            r"(?mi)^[ \t]*suffix[ \t]+([^#\s]+)",
            input_path.read_text(encoding="utf-8"),
        )
        if match:
            suffix = match.group(1)
    return stage_dir / f"OUT.{suffix}", suffix


def reuse_reference_charge(
    reference_dir: str | Path,
    target_dir: str | Path,
    *,
    overwrite: bool = False,
) -> list[Path]:
    """Copy converged reference charge data into the target ABACUS output folder."""

    source_files = _reference_charge_files(Path(reference_dir))
    if not source_files:
        raise FileNotFoundError(
            "No non-empty charge-density cube or restart file was found under "
            f"{Path(reference_dir) / 'OUT.*'}"
        )
    target, suffix = _abacus_output_dir(Path(target_dir))
    if target.is_symlink():
        raise ValueError(f'Charge output directory is a symlink: {target}. Use a private writable directory.')
    target.mkdir(parents=True, exist_ok=True)
    copied: list[Path] = []
    for source in source_files:
        destination_name = source.name
        if source.name.endswith("-CHARGE-DENSITY.restart"):
            destination_name = f"{suffix}-CHARGE-DENSITY.restart"
        destination = target / destination_name
        if destination.is_symlink():
            raise ValueError(f'Charge destination is a symlink: {destination}. Replace it with a private copy before running ABACUS.')
        if destination.exists() and not overwrite:
            continue
        shutil.copy2(source, destination)
        copied.append(destination)
    return copied


def _run_shell(
    command: str,
    *,
    cwd: Path,
    env: dict[str, str],
    log_path: Path,
    dry_run: bool,
) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"\n[{_utc_now()}] cwd={cwd}\n$ {command}\n")
        log.flush()
        if dry_run:
            return
        result = subprocess.run(
            command,
            cwd=cwd,
            env=env,
            shell=True,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {result.returncode}: {command}"
        )


def _set_abacus_parameter(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"(?m)^([ \t]*){re.escape(key)}(?:[ \t]+.*)?$")
    match = pattern.search(text)
    if match:
        return pattern.sub(
            f"{match.group(1)}{key:<20}{value}", text, count=1
        )
    lines = text.splitlines()
    insert_at = 1 if lines and lines[0].strip() == "INPUT_PARAMETERS" else 0
    lines.insert(insert_at, f"{key:<20}{value}")
    return "\n".join(lines) + ("\n" if text.endswith("\n") else "")


def _prepare_abacus_input(
    stage_dir: Path,
    *,
    reference: bool,
    dimensionality: int = 3,
) -> None:
    source = stage_dir / "INPUT-scf"
    destination = stage_dir / "INPUT"
    if destination.is_symlink():
        raise ValueError(f'INPUT is a symlink: {destination}. Use a private input copy.')
    if not source.is_file():
        if destination.is_file():
            text = destination.read_text(encoding="utf-8")
        else:
            raise FileNotFoundError(
                f"Neither INPUT-scf nor INPUT exists in {stage_dir}"
            )
    else:
        text = source.read_text(encoding="utf-8")
    if not reference:
        text = _set_abacus_parameter(text, "init_chg", "file")
    charge_output = "1 10" if int(dimensionality) in (1, 2) else "1"
    text = _set_abacus_parameter(text, "out_chg", charge_output)
    text = _set_abacus_parameter(text, "out_mat_hs2", "1")
    text = _set_abacus_parameter(text, "out_mat_r", "1")
    destination.write_bytes(text.replace("\r\n", "\n").encode("utf-8"))


def _pyatb_input_command(
    executable: str,
    *,
    reference: bool,
    dimensionality: int,
    mp_density: float,
    legacy_omega_max: float,
    polarization: bool = True,
    output: Optional[str] = None,
) -> str:
    parts = [shlex.quote(executable)]
    if polarization:
        parts.append("--polar")
    parts.extend(["--mp", f"{mp_density:g}"])
    if int(dimensionality) == 1:
        # PYATB's MP-grid parser uses a compact periodicity mask.
        parts.extend(["--dim", "001"])
    elif int(dimensionality) == 2:
        parts.extend(["--dim", "2"])
    if reference:
        parts.extend(
            ["--optical", "--orange", "0.0", f"{legacy_omega_max:g}"]
        )
    if output:
        parts.extend(["--output", shlex.quote(output)])
    return " ".join(parts)


def _pyatb_band_input_command(
    executable: str,
    *,
    gap_mode: str,
    dimensionality: int,
    mp_density: float,
    output: str = "pyatb-band",
) -> str:
    parts = [shlex.quote(executable), "--band"]
    if gap_mode == "mp":
        parts.extend(["--kmode", "mp", "--mp", f"{mp_density:g}"])
    elif gap_mode != "path":
        raise ValueError("gap_mode must be path or mp")
    if int(dimensionality) == 1:
        # PYATB's ASE band-path helper expects a space-separated PBC mask.
        parts.extend(["--dim", shlex.quote("0 0 1")])
    elif int(dimensionality) == 2:
        parts.extend(["--dim", "2"])
    parts.extend(["--output", shlex.quote(output)])
    return " ".join(parts)


def prepare_pyatb_assets(stage_dir: str | Path, pyatb_dir: str | Path) -> list[Path]:
    """Copy STRU-referenced pseudopotentials and orbitals beside PYATB Input."""

    stage = Path(stage_dir).resolve()
    destination = Path(pyatb_dir).resolve()
    stru = stage / "STRU"
    if not stru.is_file() or not destination.is_dir():
        return []

    section_names = {
        "ATOMIC_SPECIES",
        "NUMERICAL_ORBITAL",
        "ABFS_ORBITAL",
        "LATTICE_CONSTANT",
        "LATTICE_VECTORS",
        "LATTICE_PARAMETERS",
        "ATOMIC_POSITIONS",
        "NUMERICAL_DESCRIPTOR",
    }
    assets: list[str] = []
    section: Optional[str] = None
    for raw_line in stru.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        first = line.split()[0].upper()
        if first in section_names:
            section = first
            continue
        fields = line.split()
        if section == "ATOMIC_SPECIES" and len(fields) >= 3:
            assets.append(fields[2])
        elif section in {"NUMERICAL_ORBITAL", "ABFS_ORBITAL"}:
            assets.append(fields[0])

    copied: list[Path] = []
    for asset in dict.fromkeys(assets):
        source = Path(asset).expanduser()
        if not source.is_absolute():
            source = stage / source
        if not source.is_file():
            continue
        target = destination / source.name
        if source.resolve() != target.resolve():
            shutil.copy2(source, target)
        copied.append(target)
    return copied


def _run_insulation_gate(
    reference_dir: Path,
    *,
    pyatb_input: str,
    pyatb_command: str,
    gap_mode: str,
    dimensionality: int,
    mp_density: float,
    min_gap_eV: float,
    env: dict[str, str],
    log_path: Path,
    dry_run: bool,
):
    """Run the one-time band-gap gate for the undisplaced reference."""

    band_dir = reference_dir / "pyatb-band"
    if not band_is_complete(reference_dir):
        if not (band_dir / "Input").is_file():
            _run_shell(
                _pyatb_band_input_command(
                    pyatb_input,
                    gap_mode=gap_mode,
                    dimensionality=dimensionality,
                    mp_density=mp_density,
                ),
                cwd=reference_dir,
                env=env,
                log_path=log_path,
                dry_run=dry_run,
            )
        if not dry_run:
            prepare_pyatb_assets(reference_dir, band_dir)
        _run_shell(
            pyatb_command,
            cwd=band_dir,
            env=env,
            log_path=log_path,
            dry_run=dry_run,
        )
    if dry_run:
        return None

    gap = read_band_gap(band_dir, threshold_eV=min_gap_eV)
    gap_report = {
        **gap.to_dict(),
        "sampling": gap_mode,
        "dimensionality": int(dimensionality),
        "scope": "0.no-move reference gate",
    }
    (reference_dir / "zstar_insulation.json").write_text(
        json.dumps(gap_report, indent=2), encoding="utf-8"
    )
    if not gap.insulating:
        raise RuntimeError(
            "Insulating-path check failed in 0.no-move: "
            f"PYATB gap {gap.gap_eV:.6f} eV is below "
            f"{min_gap_eV:.6f} eV. Displaced calculations were not started."
        )
    return gap


def run_serial_workflow(
    root: str | Path = ".",
    *,
    abacus_command: str = "mpirun -np 1 abacus",
    pyatb_input: str = "pyatb_input",
    pyatb_command: str = "mpirun -np 1 pyatb",
    pyatb_executable: str = "pyatb",
    mp_density: float = 0.08,
    electronic_dielectric: bool = True,
    check_insulating: bool = True,
    gap_mode: str = "path",
    dimensionality: int = 3,
    min_gap_eV: float = 0.01,
    legacy_omega_max: float = DEFAULT_LEGACY_OMEGA_MAX_EV,
    legacy_domega: float = DEFAULT_LEGACY_DOMEGA_EV,
    omp_threads: int = 1,
    dry_run: bool = False,
    stop_after: Optional[int] = None,
) -> list[StageState]:
    """Execute all Born stages in deterministic serial order."""

    root_path = Path(root).resolve()
    stages = discover_stages(root_path)
    from .shared_abacus import MANIFEST, load_manifest, read_forces
    shared = (root_path / MANIFEST).is_file()
    if shared:
        dimensionality = load_manifest(root_path)["dimension"]
        from .pyatb_precision import precision_command
        pyatb_command = precision_command(pyatb_command, pyatb_executable)
    store = WorkflowStateStore(root_path)
    reference_dir = stages[0].path
    caps = detect_pyatb_capabilities(pyatb_executable)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(int(omp_threads))
    env.setdefault("MKL_NUM_THREADS", str(int(omp_threads)))
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    results: list[StageState] = []

    store.event(
        "workflow",
        "start",
        root=str(root_path),
        stages=len(stages),
        dry_run=dry_run,
        pyatb=caps.to_dict(),
    )

    for index, stage in enumerate(stages):
        if stop_after is not None and index >= stop_after:
            break
        state = store.load(stage)
        state.started_at = state.started_at or _utc_now()
        state.status = "running"
        state.error = None
        store.save(state)
        store.event(stage.name, "stage-start", path=str(stage.path))
        log_path = store.state_dir / "logs" / f"{store._safe_name(stage.name)}.log"

        try:
            if scf_is_complete(stage.path):
                state.scf = "completed"
                store.event(stage.name, "scf-skip", reason="completion marker found")
            else:
                _prepare_abacus_input(
                    stage.path,
                    reference=stage.reference,
                    dimensionality=dimensionality,
                )
                if not stage.reference and not dry_run:
                    copied = reuse_reference_charge(reference_dir, stage.path)
                    store.event(
                        stage.name,
                        "charge-reused",
                        files=[str(path) for path in copied],
                    )
                state.scf = "running"
                store.save(state)
                _run_shell(
                    abacus_command,
                    cwd=stage.path,
                    env=env,
                    log_path=log_path,
                    dry_run=dry_run,
                )
                if not dry_run and not scf_is_complete(stage.path):
                    raise RuntimeError(
                        "ABACUS returned successfully but running_scf.log has no completion marker"
                    )
                state.scf = "dry-run" if dry_run else "completed"

            if shared and not dry_run:
                read_forces(stage.path)

            if check_insulating and stage.reference:
                state.band = "running"
                store.save(state)
                gap = _run_insulation_gate(
                    reference_dir,
                    pyatb_input=pyatb_input,
                    pyatb_command=pyatb_command,
                    gap_mode=gap_mode,
                    dimensionality=dimensionality,
                    mp_density=mp_density,
                    min_gap_eV=min_gap_eV,
                    env=env,
                    log_path=log_path,
                    dry_run=dry_run,
                )
                if dry_run:
                    state.band = "dry-run"
                else:
                    state.band_gap_eV = gap.gap_eV
                    state.band = "insulating"
                    store.event(
                        stage.name,
                        "band-gap",
                        gap_eV=gap.gap_eV,
                        sampling=gap_mode,
                        threshold_eV=min_gap_eV,
                        insulating=True,
                    )
            elif check_insulating:
                state.band = "reference-gated"
            else:
                state.band = "not-requested"

            need_dielectric = bool(stage.reference and electronic_dielectric)
            polar_done = polarization_is_complete(stage.path)
            dielectric_done = dielectric_is_complete(stage.path) if need_dielectric else True
            if polar_done and dielectric_done:
                state.pyatb = "completed"
                state.dielectric = "completed" if need_dielectric else "not-requested"
                store.event(stage.name, "pyatb-skip", reason="requested outputs found")
            else:
                pyatb_dir = stage.path / "pyatb"
                pyatb_input_path = pyatb_dir / "Input"
                if not pyatb_input_path.is_file():
                    command = _pyatb_input_command(
                        pyatb_input,
                        reference=need_dielectric,
                        dimensionality=dimensionality,
                        mp_density=mp_density,
                        legacy_omega_max=legacy_omega_max,
                    )
                    _run_shell(
                        command,
                        cwd=stage.path,
                        env=env,
                        log_path=log_path,
                        dry_run=dry_run,
                    )
                if need_dielectric and not dry_run:
                    configure_optical_input(
                        pyatb_input_path,
                        capabilities=caps,
                        static_only=True,
                        legacy_omega_max=legacy_omega_max,
                        legacy_domega=legacy_domega,
                    )
                if not dry_run:
                    configure_polarization_input(
                        pyatb_input_path,
                        dimensionality=dimensionality,
                    )
                    prepare_pyatb_assets(stage.path, pyatb_dir)
                state.pyatb = "running"
                state.dielectric = "running" if need_dielectric else "not-requested"
                store.save(state)
                _run_shell(
                    pyatb_command,
                    cwd=pyatb_dir,
                    env=env,
                    log_path=log_path,
                    dry_run=dry_run,
                )
                if not dry_run and not polarization_is_complete(stage.path):
                    raise RuntimeError("PYATB did not produce Polarization/polarization.dat")
                if not dry_run and need_dielectric and not dielectric_is_complete(stage.path):
                    raise RuntimeError("PYATB did not produce a static dielectric output")
                state.pyatb = "dry-run" if dry_run else "completed"
                if need_dielectric:
                    state.dielectric = "dry-run" if dry_run else "completed"

            state.status = "dry-run" if dry_run else "completed"
            state.finished_at = _utc_now()
            store.event(stage.name, "stage-complete", status=state.status)
        except Exception as exc:
            state.status = "failed"
            state.error = str(exc)
            state.finished_at = _utc_now()
            store.event(stage.name, "stage-failed", error=str(exc))
            store.save(state)
            results.append(state)
            raise
        store.save(state)
        results.append(state)

    store.event("workflow", "finish", completed=len(results))
    return results


def discover_raman_stages(raman_dir: str | Path) -> list[WorkflowStage]:
    """Return the +/- Raman finite-difference stages in manifest order."""

    root = Path(raman_dir).resolve()
    manifest_path = root / "raman_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Raman manifest not found: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    stages: list[WorkflowStage] = []
    for entry in manifest.get("modes", []):
        mode = int(entry["mode"])
        for sign in ("plus", "minus"):
            path = Path(entry[sign]).resolve()
            if not path.is_dir():
                raise FileNotFoundError(
                    f"Raman mode {mode} {sign} directory not found: {path}"
                )
            stages.append(
                WorkflowStage(
                    name=f"mode-{mode:04d}/{sign}",
                    path=path,
                )
            )
    if not stages:
        raise ValueError(f"No Raman stages found in {manifest_path}")
    return stages


def run_raman_workflow(
    raman_dir: str | Path = "raman",
    *,
    reference_dir: str | Path = "0.no-move",
    abacus_command: str = "mpirun -np 1 abacus",
    pyatb_input: str = "pyatb_input",
    pyatb_command: str = "mpirun -np 1 pyatb",
    pyatb_executable: str = "pyatb",
    mp_density: float = 0.08,
    check_insulating: bool = True,
    gap_mode: str = "path",
    dimensionality: int = 3,
    molecular_ir: bool = False,
    min_gap_eV: float = 0.01,
    legacy_omega_max: float = DEFAULT_LEGACY_OMEGA_MAX_EV,
    legacy_domega: float = DEFAULT_LEGACY_DOMEGA_EV,
    omp_threads: int = 1,
    dry_run: bool = False,
    stop_after: Optional[int] = None,
) -> list[StageState]:
    """Run Raman +/- structures serially and calculate response tensors."""

    root = Path(raman_dir).resolve()
    reference = Path(reference_dir).resolve()
    stages = discover_raman_stages(root)
    manifest_path = root / "raman_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["reference_dir"] = str(reference)
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    store = WorkflowStateStore(root)
    caps = detect_pyatb_capabilities(pyatb_executable)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(int(omp_threads))
    env.setdefault("MKL_NUM_THREADS", str(int(omp_threads)))
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    results: list[StageState] = []

    if not dry_run and not _reference_charge_files(reference):
        raise FileNotFoundError(
            f"No reusable reference charge-density cube found under {reference / 'OUT.*'}"
        )

    store.event(
        "raman-workflow",
        "start",
        root=str(root),
        reference=str(reference),
        stages=len(stages),
        dry_run=dry_run,
        pyatb=caps.to_dict(),
    )

    if check_insulating:
        reference_log = store.state_dir / "logs" / "0.no-move.log"
        gap = _run_insulation_gate(
            reference,
            pyatb_input=pyatb_input,
            pyatb_command=pyatb_command,
            gap_mode=gap_mode,
            dimensionality=dimensionality,
            mp_density=mp_density,
            min_gap_eV=min_gap_eV,
            env=env,
            log_path=reference_log,
            dry_run=dry_run,
        )
        if gap is not None:
            store.event(
                "0.no-move",
                "band-gap",
                gap_eV=gap.gap_eV,
                sampling=gap_mode,
                threshold_eV=min_gap_eV,
                insulating=True,
            )

    for index, stage in enumerate(stages):
        if stop_after is not None and index >= stop_after:
            break
        state = store.load(stage)
        state.started_at = state.started_at or _utc_now()
        state.status = "running"
        state.error = None
        state.dielectric = "pending"
        store.save(state)
        store.event(stage.name, "stage-start", path=str(stage.path))
        log_path = (
            store.state_dir / "logs" / f"{store._safe_name(stage.name)}.log"
        )

        try:
            if scf_is_complete(stage.path):
                state.scf = "completed"
                store.event(stage.name, "scf-skip", reason="completion marker found")
            else:
                _prepare_abacus_input(stage.path, reference=False)
                if not dry_run:
                    copied = reuse_reference_charge(reference, stage.path)
                    store.event(
                        stage.name,
                        "charge-reused",
                        files=[str(path) for path in copied],
                    )
                state.scf = "running"
                store.save(state)
                _run_shell(
                    abacus_command,
                    cwd=stage.path,
                    env=env,
                    log_path=log_path,
                    dry_run=dry_run,
                )
                if not dry_run and not scf_is_complete(stage.path):
                    raise RuntimeError(
                        "ABACUS returned successfully but running_scf.log has no "
                        "completion marker"
                    )
                state.scf = "dry-run" if dry_run else "completed"

            if check_insulating:
                state.band = "reference-gated"
            else:
                state.band = "not-requested"

            if dielectric_is_complete(stage.path):
                store.event(
                    stage.name,
                    "pyatb-optical-skip",
                    reason="static dielectric output found",
                )
            else:
                pyatb_dir = stage.path / "pyatb"
                input_path = pyatb_dir / "Input"
                if not input_path.is_file():
                    _run_shell(
                        _pyatb_input_command(
                            pyatb_input,
                            reference=True,
                            dimensionality=dimensionality,
                            mp_density=mp_density,
                            legacy_omega_max=legacy_omega_max,
                            polarization=False,
                        ),
                        cwd=stage.path,
                        env=env,
                        log_path=log_path,
                        dry_run=dry_run,
                    )
                if not dry_run:
                    configure_optical_input(
                        input_path,
                        capabilities=caps,
                        static_only=True,
                        legacy_omega_max=legacy_omega_max,
                        legacy_domega=legacy_domega,
                    )
                    prepare_pyatb_assets(stage.path, pyatb_dir)
                state.pyatb = "running"
                state.dielectric = "running"
                store.save(state)
                _run_shell(
                    pyatb_command,
                    cwd=pyatb_dir,
                    env=env,
                    log_path=log_path,
                    dry_run=dry_run,
                )
                if not dry_run and not dielectric_is_complete(stage.path):
                    raise RuntimeError(
                        "PYATB did not produce a static dielectric output"
                    )
            state.dielectric = "dry-run" if dry_run else "completed"

            if molecular_ir and not polarization_is_complete(
                stage.path, "pyatb-polar"
            ):
                polar_dir = stage.path / "pyatb-polar"
                polar_input = polar_dir / "Input"
                if not polar_input.is_file():
                    _run_shell(
                        _pyatb_input_command(
                            pyatb_input,
                            reference=False,
                            dimensionality=dimensionality,
                            mp_density=mp_density,
                            legacy_omega_max=legacy_omega_max,
                            polarization=True,
                            output="pyatb-polar",
                        ),
                        cwd=stage.path,
                        env=env,
                        log_path=log_path,
                        dry_run=dry_run,
                    )
                if not dry_run:
                    prepare_pyatb_assets(stage.path, polar_dir)
                state.pyatb = "running"
                store.save(state)
                _run_shell(
                    pyatb_command,
                    cwd=polar_dir,
                    env=env,
                    log_path=log_path,
                    dry_run=dry_run,
                )
                if not dry_run and not polarization_is_complete(
                    stage.path, "pyatb-polar"
                ):
                    raise RuntimeError(
                        "PYATB did not produce molecular Polarization/polarization.dat"
                    )
            elif molecular_ir:
                store.event(
                    stage.name,
                    "pyatb-polarization-skip",
                    reason="molecular polarization output found",
                )
            state.pyatb = "dry-run" if dry_run else "completed"

            state.status = "dry-run" if dry_run else "completed"
            state.finished_at = _utc_now()
            store.event(stage.name, "stage-complete", status=state.status)
        except Exception as exc:
            state.status = "failed"
            state.error = str(exc)
            state.finished_at = _utc_now()
            store.event(stage.name, "stage-failed", error=str(exc))
            store.save(state)
            results.append(state)
            raise
        store.save(state)
        results.append(state)

    store.event("raman-workflow", "finish", completed=len(results))
    return results


def workflow_status(root: str | Path = ".") -> list[StageState]:
    root_path = Path(root).resolve()
    from .shared_abacus import MANIFEST, read_forces
    shared = (root_path / MANIFEST).is_file()
    store = WorkflowStateStore(root_path)
    reference_gated = (root_path / "0.no-move" / "zstar_insulation.json").is_file()
    states: list[StageState] = []
    for stage in discover_stages(root_path):
        state = store.load(stage)
        if state.status == "pending":
            if scf_is_complete(stage.path):
                state.scf = "completed"
            if stage.reference and band_is_complete(stage.path):
                try:
                    gap = read_band_gap(Path(stage.path) / "pyatb-band")
                    state.band_gap_eV = gap.gap_eV
                    state.band = "insulating" if gap.insulating else "metallic"
                except (FileNotFoundError, ValueError):
                    pass
            elif not stage.reference and reference_gated:
                state.band = "reference-gated"
            if polarization_is_complete(stage.path):
                state.pyatb = "completed"
            if stage.reference and dielectric_is_complete(stage.path):
                state.dielectric = "completed"
            if state.scf == "completed" and state.pyatb == "completed":
                state.status = "completed"
        if shared and state.scf == 'completed':
            try:
                read_forces(stage.path)
            except ValueError as exc:
                state.status = 'failed'
                state.error = str(exc)
        states.append(state)
    return states


def raman_workflow_status(raman_dir: str | Path = "raman") -> list[StageState]:
    root = Path(raman_dir).resolve()
    store = WorkflowStateStore(root)
    reference_value = json.loads(
        (root / "raman_manifest.json").read_text(encoding="utf-8")
    ).get("reference_dir")
    reference = (
        Path(reference_value)
        if reference_value
        else root.parent / "0.no-move"
    )
    reference_gated = (reference / "zstar_insulation.json").is_file()
    states: list[StageState] = []
    for stage in discover_raman_stages(root):
        state = store.load(stage)
        if state.status == "pending":
            if scf_is_complete(stage.path):
                state.scf = "completed"
            if reference_gated:
                state.band = "reference-gated"
            if dielectric_is_complete(stage.path):
                state.pyatb = "completed"
                state.dielectric = "completed"
            if state.scf == "completed" and state.dielectric == "completed":
                state.status = "completed"
        states.append(state)
    return states


def write_workflow_manifest(root: str | Path = ".") -> Path:
    root_path = Path(root).resolve()
    stages = discover_stages(root_path)
    state_dir = root_path / ".zstar"
    state_dir.mkdir(parents=True, exist_ok=True)
    path = state_dir / "workflow_manifest.json"
    data = {
        "schema": 1,
        "root": str(root_path),
        "created_at": _utc_now(),
        "execution": "serial",
        "charge_reference": (
            "0.no-move/OUT.*/SPIN*_CHG.cube or "
            "0.no-move/OUT.*/*-CHARGE-DENSITY.restart"
        ),
        "insulation_gate": "0.no-move only",
        "stages": [
            {
                "index": index,
                "name": stage.name,
                "path": str(stage.path),
                "reference": stage.reference,
            }
            for index, stage in enumerate(stages)
        ],
    }
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return path


def _quote_bash(value: str) -> str:
    return shlex.quote(str(value))


def generate_backend_script(
    root: str | Path = ".",
    *,
    backend: str = "shell",
    output: Optional[str | Path] = None,
    job_name: str = "zstar-born",
    nodes: int = 1,
    tasks: Optional[int] = None,
    cpus_per_task: Optional[int] = None,
    walltime: str = "24:00:00",
    queue: Optional[str] = None,
    account: Optional[str] = None,
    env_script: Optional[str] = None,
    header_file: Optional[str | Path] = None,
    abacus_command: Optional[str] = None,
    pyatb_command: Optional[str] = None,
    mp_density: float = 0.08,
    check_insulating: bool = True,
    gap_mode: str = "path",
    dimensionality: int = 3,
    min_gap_eV: float = 0.01,
    legacy_omega_max: float = DEFAULT_LEGACY_OMEGA_MAX_EV,
    dry_run: bool = False,
) -> Path:
    """Generate a single backend script that runs the complete serial chain."""
    from .configuration import launcher_command, resolve_parallelism
    from .job_headers import compose_job_script, torque_ppn
    tasks, cpus_per_task = resolve_parallelism(root, tasks=tasks, cpus_per_task=cpus_per_task)
    backend_key = normalize_execution_system(backend)
    if int(nodes) < 1 or int(tasks) < 1 or int(cpus_per_task) < 1:
        raise ValueError("nodes, tasks, and cpus_per_task must be positive")
    root_path = Path(root).resolve()
    if output is None:
        suffix = {"shell": "sh", "slurm": "slurm", "torque": "pbs"}[backend_key]
        output_path = root_path / f"run_zstar_born.{suffix}"
    else:
        output_path = Path(output).resolve()

    if abacus_command is None:
        abacus_command = launcher_command('abacus', root=root, system=backend_key, tasks=tasks)
    if pyatb_command is None:
        pyatb_command = launcher_command('pyatb', root=root, system=backend_key, tasks=tasks)
    if env_script is not None:
        env_script = str(Path(env_script).expanduser().resolve())

    insulation_option = (
        f" --min-gap {float(min_gap_eV):g}"
        if check_insulating
        else " --no-insulation-check"
    )
    run_line = (
        "zstar workflow run"
        f" --root {_quote_bash(str(root_path))}"
        f" --abacus-command {_quote_bash(abacus_command)}"
        f" --pyatb-command {_quote_bash(pyatb_command)}"
        f" --mp-density {float(mp_density):g}"
        f" --gap-mode {_quote_bash(gap_mode)}"
        f" --dimensionality {int(dimensionality)}"
        f"{insulation_option}"
        f" --legacy-omega-max {float(legacy_omega_max):g}"
        f" --omp-threads {int(cpus_per_task)}"
        f"{' --dry-run' if dry_run else ''}"
        ' 2>&1 | tee -a "$ROOT/.zstar/workflow.log"'
    )

    header = [
        "#!/usr/bin/env bash",
        f"# ZStar execution backend: {backend_key}",
        "# Displacement stages run serially and resume from .zstar/stages/*.json.",
    ]
    if backend_key == "slurm":
        header.extend(
            [
                f"#SBATCH --job-name={job_name}",
                f"#SBATCH --nodes={int(nodes)}",
                f"#SBATCH --ntasks={int(tasks)}",
                f"#SBATCH --cpus-per-task={int(cpus_per_task)}",
                f"#SBATCH --time={walltime}",
                f"#SBATCH --output={root_path}/.zstar/slurm-%j.out",
            ]
        )
        if queue:
            header.append(f"#SBATCH --partition={queue}")
        if account:
            header.append(f"#SBATCH --account={account}")
    elif backend_key == "torque":
        ppn = torque_ppn(nodes, tasks, cpus_per_task)
        header.extend(
            [
                f"#PBS -N {job_name}",
                f"#PBS -l nodes={int(nodes)}:ppn={ppn}",
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
        f"ROOT={_quote_bash(str(root_path))}",
        'mkdir -p "$ROOT/.zstar"',
    ]
    body.append('cd "$ROOT"')
    if env_script:
        body.append(f"source {_quote_bash(env_script)}")
    body.extend(
        [
            f"export OMP_NUM_THREADS={int(cpus_per_task)}",
            f"export MKL_NUM_THREADS={int(cpus_per_task)}",
            "export OPENBLAS_NUM_THREADS=1",
            'printf "[%s] ZStar serial Born workflow starts on %s\\n" "$(date -Is)" "$(hostname)" '
            '| tee -a "$ROOT/.zstar/workflow.log"',
            run_line,
            'printf "[%s] ZStar serial Born workflow finished\\n" "$(date -Is)" '
            '| tee -a "$ROOT/.zstar/workflow.log"',
        ]
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    script_text = compose_job_script(root_path, backend_key, header, body, specified=header_file)
    output_path.write_bytes(script_text.encode("utf-8"))
    if os.name != "nt":
        output_path.chmod(output_path.stat().st_mode | 0o111)
    write_workflow_manifest(root_path)
    header_record = json.loads((root_path / '.zstar' / 'job_header.json').read_text(encoding='utf-8'))
    backend_manifest = {
        "schema": 1,
        "backend": backend_key,
        "script": str(output_path),
        "execution": "serial",
        "resume_state": str(root_path / ".zstar" / "stages"),
        "resources": {
            "nodes": int(nodes),
            "tasks": int(tasks),
            "cpus_per_task": int(cpus_per_task),
            "walltime": walltime,
            "queue": queue,
            "account": account,
        },
        "commands": {
            "abacus": abacus_command,
            "pyatb": pyatb_command,
        },
        "environment_script": env_script,
        "job_header": header_record,
        "resources_source": "generated defaults" if header_record['level'] == 'Default' else "selected header",
        "dry_run": bool(dry_run),
    }
    if header_record['level'] != 'Default':
        backend_manifest['resources'] = {'tasks': int(tasks), 'cpus_per_task': int(cpus_per_task)}
    (root_path / ".zstar" / "backend_manifest.json").write_text(
        json.dumps(backend_manifest, indent=2), encoding="utf-8"
    )
    return output_path


def submit_backend_script(script: str | Path, backend: str) -> str:
    backend_key = normalize_execution_system(backend)
    command = {
        "slurm": ["sbatch", str(Path(script).resolve())],
        "torque": ["qsub", str(Path(script).resolve())],
        "shell": ["bash", str(Path(script).resolve())],
    }.get(backend_key)
    if command is None:
        raise ValueError("backend must be one of: shell, slurm, torque")
    result = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        cwd=Path(script).resolve().parent,
    )
    return result.stdout.strip()


def format_status_table(states: Iterable[StageState]) -> str:
    rows = [
        (
            state.name,
            state.status,
            state.scf,
            state.band,
            "" if state.band_gap_eV is None else f"{state.band_gap_eV:.6f}",
            state.pyatb,
            state.dielectric,
            state.error or "",
        )
        for state in states
    ]
    headers = (
        "stage",
        "status",
        "scf",
        "band",
        "gap_eV",
        "pyatb",
        "dielectric",
        "error",
    )
    widths = [
        max([len(headers[i]), *(len(str(row[i])) for row in rows)])
        for i in range(len(headers))
    ]
    lines = [
        "  ".join(headers[i].ljust(widths[i]) for i in range(len(headers))),
        "  ".join("-" * width for width in widths),
    ]
    lines.extend(
        "  ".join(str(row[i]).ljust(widths[i]) for i in range(len(headers)))
        for row in rows
    )
    return "\n".join(lines)
