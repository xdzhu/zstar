"""Calculator-native IR and Raman workflows for VASP and CP2K."""

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
from typing import Iterable, Sequence

import numpy as np

from .cp2k_bec import (
    _find_section,
    _has_section,
    _insert_before_section_end,
    _referenced_assets,
    _set_section_keyword,
)
from .spectra import (
    BornData,
    GammaModes,
    calculate_ir_spectrum,
    calculate_molecular_ir_spectrum,
    calculate_native_line_spectrum,
    calculate_raman_spectrum,
    mode_effective_charges,
    validate_frequencies_stable,
    validate_gamma_stability,
    write_ir_outputs,
    write_molecular_ir_outputs,
    write_native_line_spectrum_outputs,
    write_raman_outputs,
)
from .vasp_bec import parse_vasp_gap, parse_vasp_outcar, render_incar, vasp_output_complete


DEBYE_PER_E_ANGSTROM = 4.80320471257
BOHR3_TO_ANGSTROM3 = 0.14818471147216278
_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


@dataclass
class SpectraStageState:
    name: str
    path: str
    status: str = "pending"
    started_at: str | None = None
    finished_at: str | None = None
    output: str | None = None
    error: str | None = None
    gap_eV: float | None = None


def _mode_indices(
    frequencies_cm1: np.ndarray,
    mode_numbers: Sequence[int] | None,
    acoustic_cutoff_cm1: float,
) -> np.ndarray:
    if mode_numbers:
        indices = np.asarray([int(number) - 1 for number in mode_numbers], dtype=int)
        if np.any(indices < 0) or np.any(indices >= len(frequencies_cm1)):
            raise IndexError("Requested mode number is outside the eigensystem")
        return indices
    return np.flatnonzero(frequencies_cm1 > float(acoustic_cutoff_cm1))


def _real_modes(eigenvectors: np.ndarray) -> np.ndarray:
    output = np.empty(eigenvectors.shape, dtype=float)
    for index, vector in enumerate(eigenvectors):
        flat = vector.reshape(-1)
        phase = np.angle(flat[int(np.argmax(np.abs(flat)))])
        phased = vector * np.exp(-1j * phase)
        if np.linalg.norm(phased.imag) > max(1.0e-7, 1.0e-5 * np.linalg.norm(phased.real)):
            raise ValueError(f"Mode {index + 1} cannot be represented by a real Gamma vector")
        output[index] = phased.real
    return output


def load_vasp_gamma_modes(path: str | Path) -> GammaModes:
    """Load Gamma modes from a VASP vibrational ``vasprun.xml``."""

    source = Path(path).resolve()
    try:
        from pymatgen.io.vasp.outputs import Vasprun

        run = Vasprun(
            source,
            parse_dos=False,
            parse_eigen=False,
            parse_projected_eigen=False,
            parse_potcar_file=False,
            exception_on_bad_xml=True,
        )
    except Exception as exc:
        raise ValueError(f"Cannot read VASP vibrational data from {source}: {exc}") from exc
    eigenvalues = np.asarray(run.normalmode_eigenvals, dtype=float)
    eigenvectors = np.asarray(run.normalmode_eigenvecs, dtype=complex)
    if eigenvalues.ndim != 1 or eigenvectors.shape[0] != len(eigenvalues):
        raise ValueError(f"No complete VASP normal modes found in {source}")
    # pymatgen stores VASP's stable modes as negative omega^2 in THz^2.
    frequencies = np.sign(-eigenvalues) * np.sqrt(np.abs(eigenvalues))
    structure = run.final_structure
    symbols = tuple(str(site.specie.symbol) for site in structure)
    masses = np.asarray([float(site.specie.atomic_mass) for site in structure])
    return GammaModes(
        frequencies_thz=frequencies,
        eigenvectors=eigenvectors,
        masses_amu=masses,
        lattice_angstrom=np.asarray(structure.lattice.matrix, dtype=float),
        symbols=symbols,
        positions_fractional=np.asarray(structure.frac_coords, dtype=float),
    )


def _write_vasp_poscar(
    destination: Path,
    modes: GammaModes,
    displacement_angstrom: np.ndarray | None = None,
) -> None:
    from pymatgen.core import Lattice, Structure
    from pymatgen.io.vasp.inputs import Poscar

    structure = Structure(
        Lattice(modes.lattice_angstrom),
        list(modes.symbols),
        modes.positions_fractional,
        coords_are_cartesian=False,
    )
    if displacement_angstrom is not None:
        for index, vector in enumerate(displacement_angstrom):
            structure.translate_sites(index, vector, frac_coords=False, to_unit_cell=True)
    Poscar(structure).write_file(destination)


def prepare_vasp_spectra(
    input_dir: str | Path,
    modes_xml: str | Path,
    root: str | Path = "vasp_spectra",
    *,
    amplitude: float = 0.02,
    mode_numbers: Sequence[int] | None = None,
    acoustic_cutoff_cm1: float = 5.0,
    imaginary_tolerance_cm1: float = 20.0,
    allow_imaginary: bool = False,
    method: str = "dfpt",
    field_strength: float = 0.001,
    dimensionality: int = 3,
    force: bool = False,
) -> Path:
    """Prepare one reference and +/- VASP dielectric calculations per mode."""

    if dimensionality == 2:
        raise ValueError(
            "VASP 2D out-of-plane IR/Raman is not enabled: use the ZStar real-space "
            "slab polarization workflow for the vacuum direction"
        )
    if dimensionality not in {0, 1, 3}:
        raise ValueError("VASP spectroscopy dimensionality must be 0, 1, or 3")
    if amplitude <= 0.0:
        raise ValueError("amplitude must be positive")
    method_key = method.lower()
    if method_key not in {"dfpt", "finite-field"}:
        raise ValueError("method must be dfpt or finite-field")
    if dimensionality == 1 and method_key != "dfpt":
        raise ValueError("VASP 1D spectroscopy currently requires method=dfpt")
    source = Path(input_dir).resolve()
    required = ("INCAR", "POSCAR", "KPOINTS", "POTCAR")
    missing = [name for name in required if not (source / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing VASP input files: {', '.join(missing)}")
    modes_source = Path(modes_xml).resolve()
    modes = load_vasp_gamma_modes(modes_source)
    validate_gamma_stability(
        modes,
        imaginary_tolerance_cm1=imaginary_tolerance_cm1,
        allow_imaginary=allow_imaginary,
    )
    indices = _mode_indices(modes.frequencies_cm1, mode_numbers, acoustic_cutoff_cm1)
    if not len(indices):
        raise ValueError("No optical modes selected")

    target = Path(root).resolve()
    if target.exists() and force:
        shutil.rmtree(target)
    if target.exists() and any(target.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {target}")
    target.mkdir(parents=True, exist_ok=True)

    original = (source / "INCAR").read_text(encoding="utf-8", errors="ignore")
    remove = ("LEPSILON", "LCALCEPS", "EFIELD_PEAD", "LOPTICS")
    common = {
        "NSW": "0",
        "IBRION": "-1",
        "ISYM": "0",
        "LREAL": ".FALSE.",
        "LCHARG": ".TRUE.",
        "LWAVE": ".TRUE.",
        "EDIFF": "1E-8",
        "LRPA": ".FALSE.",
    }
    if method_key == "dfpt":
        common["LEPSILON"] = ".TRUE."
    else:
        value = f"{field_strength:.10g}"
        common.update(
            {
                "LCALCEPS": ".TRUE.",
                "EFIELD_PEAD": f"{value} {value} {value}",
                "ISMEAR": "0",
                "SIGMA": "0.05",
            }
        )
    reference_incar = render_incar(original, updates=common, remove=remove)
    displaced_incar = render_incar(
        reference_incar,
        updates={"ISTART": "1", "ICHARG": "1"},
    )

    def write_stage(relative: Path, incar: str, displacement: np.ndarray | None) -> None:
        directory = target / relative
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "INCAR").write_text(incar, encoding="utf-8", newline="\n")
        _write_vasp_poscar(directory / "POSCAR", modes, displacement)
        for name in ("KPOINTS", "POTCAR"):
            shutil.copy2(source / name, directory / name)

    write_stage(Path("reference"), reference_incar, None)
    real_vectors = _real_modes(modes.eigenvectors)
    mass_weighted = real_vectors / np.sqrt(modes.masses_amu)[None, :, None]
    entries: list[dict] = []
    stages = [{"name": "reference", "path": "reference", "reference": True}]
    for index in indices:
        entry = {
            "mode": int(index + 1),
            "frequency_cm-1": float(modes.frequencies_cm1[index]),
            "amplitude_A_sqrt_amu": float(amplitude),
        }
        for sign_name, sign in (("plus", 1.0), ("minus", -1.0)):
            relative = Path(f"mode-{index + 1:04d}") / sign_name
            displacement = sign * float(amplitude) * mass_weighted[index]
            write_stage(relative, displaced_incar, displacement)
            entry[sign_name] = relative.as_posix()
            stages.append({"name": relative.as_posix(), "path": relative.as_posix()})
        entries.append(entry)
    manifest = {
        "schema_version": 1,
        "calculator": "vasp",
        "created_at": _utc_now(),
        "source_directory": str(source),
        "modes_source": str(modes_source),
        "dimensionality": dimensionality,
        "periodic_axes": "z" if dimensionality == 1 else ("xyz" if dimensionality == 3 else ""),
        "nac_model": "none" if dimensionality == 1 else "bulk",
        "method": method_key,
        "field_strength_eV_per_angstrom": (
            field_strength if method_key == "finite-field" else None
        ),
        "amplitude_A_sqrt_amu": amplitude,
        "imaginary_tolerance_cm-1": float(imaginary_tolerance_cm1),
        "allow_imaginary": bool(allow_imaginary),
        "modes": entries,
        "stages": stages,
    }
    (target / "spectra_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8", newline="\n"
    )
    return target


def _load_manifest(root: str | Path) -> tuple[Path, dict]:
    root_path = Path(root).resolve()
    path = root_path / "spectra_manifest.json"
    if not path.is_file():
        raise FileNotFoundError(path)
    return root_path, json.loads(path.read_text(encoding="utf-8"))


def _state_path(root: Path) -> Path:
    return root / ".zstar" / "spectra_state.json"


def _write_states(root: Path, states: Sequence[SpectraStageState]) -> None:
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


def _cp2k_output_complete(path: Path) -> bool:
    if not path.is_file() or not path.stat().st_size:
        return False
    tail = path.read_text(encoding="utf-8", errors="ignore")[-100000:]
    return "PROGRAM ENDED AT" in tail and "ABORT" not in tail


def _copy_vasp_restart(reference: Path, target: Path) -> None:
    for name in ("WAVECAR", "CHGCAR"):
        source = reference / name
        if not source.is_file() or not source.stat().st_size:
            raise FileNotFoundError(f"Required VASP restart file is missing: {source}")
        shutil.copy2(source, target / name)


def run_calculator_spectra(
    root: str | Path,
    *,
    command: str | None = None,
    min_gap_eV: float = 0.01,
    omp_threads: int = 1,
    extra_env: dict[str, str] | None = None,
    dry_run: bool = False,
    stop_after: int | None = None,
) -> list[SpectraStageState]:
    """Run a prepared VASP or CP2K spectroscopy workflow with restart state."""

    root_path, manifest = _load_manifest(root)
    calculator = manifest["calculator"]
    if command is None:
        command = "vasp_std" if calculator == "vasp" else "cp2k.ssmp -i input.inp -o output.log"
    old: dict[str, dict] = {}
    path = _state_path(root_path)
    if path.is_file():
        saved = json.loads(path.read_text(encoding="utf-8"))
        old = {item["name"]: item for item in saved.get("stages", [])}
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(int(omp_threads))
    environment.setdefault("OMP_PROC_BIND", "close")
    environment.setdefault("OMP_PLACES", "cores")
    if extra_env:
        environment.update({str(key): str(value) for key, value in extra_env.items()})
    states: list[SpectraStageState] = []

    for stage_data in manifest["stages"]:
        name = stage_data["name"]
        directory = root_path / stage_data["path"]
        previous = old.get(name, {})
        output = directory / ("OUTCAR" if calculator == "vasp" else "output.log")
        complete = vasp_output_complete(output) if calculator == "vasp" else _cp2k_output_complete(output)
        state = SpectraStageState(
            name=name,
            path=str(directory),
            status=previous.get("status", "pending"),
            started_at=previous.get("started_at"),
            finished_at=previous.get("finished_at"),
            output=str(output),
            error=previous.get("error"),
            gap_eV=previous.get("gap_eV"),
        )
        if complete:
            try:
                if calculator == "cp2k":
                    parse_cp2k_native_spectra(output)
                state.status = "completed"
                if calculator == "vasp" and name == "reference":
                    state.gap_eV = parse_vasp_gap(directory / "vasprun.xml")
                    if state.gap_eV < min_gap_eV:
                        state.status = "rejected_metal"
                        state.error = (
                            f"Reference gap {state.gap_eV:.6g} eV is below "
                            f"the {min_gap_eV:.6g} eV threshold"
                        )
            except Exception as exc:
                state.status = "failed"
                state.error = str(exc)
            states.append(state)
            _write_states(root_path, states)
            if state.status != "completed":
                break
            continue

        if calculator == "vasp" and name != "reference":
            reference = next((item for item in states if item.name == "reference"), None)
            accepted = {"completed", "dry-run"} if dry_run else {"completed"}
            if reference is None or reference.status not in accepted:
                state.status = "blocked"
                state.error = "Reference VASP response is not complete and insulating"
                states.append(state)
                _write_states(root_path, states)
                break
            if not dry_run:
                _copy_vasp_restart(root_path / "reference", directory)

        state.started_at = _utc_now()
        state.error = None
        if dry_run:
            state.status = "dry-run"
        else:
            state.status = "running"
            _write_states(root_path, [*states, state])
            try:
                log_path = directory / ("vasp.log" if calculator == "vasp" else "driver.log")
                with log_path.open("w", encoding="utf-8") as log:
                    subprocess.run(
                        command,
                        cwd=directory,
                        env=environment,
                        shell=True,
                        check=True,
                        stdout=log,
                        stderr=subprocess.STDOUT,
                    )
                complete = (
                    vasp_output_complete(output)
                    if calculator == "vasp"
                    else _cp2k_output_complete(output)
                )
                if not complete:
                    raise RuntimeError(f"{calculator.upper()} output has no normal footer")
                if calculator == "cp2k":
                    parse_cp2k_native_spectra(output)
                state.status = "completed"
                state.finished_at = _utc_now()
                if calculator == "vasp" and name == "reference":
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
        _write_states(root_path, states)
        if state.status not in {"completed", "dry-run"}:
            break
        if stop_after is not None and len(states) >= int(stop_after):
            break
    return states


def calculator_spectra_status(root: str | Path) -> list[SpectraStageState]:
    root_path, manifest = _load_manifest(root)
    saved: dict[str, dict] = {}
    path = _state_path(root_path)
    if path.is_file():
        data = json.loads(path.read_text(encoding="utf-8"))
        saved = {item["name"]: item for item in data.get("stages", [])}
    states: list[SpectraStageState] = []
    for item in manifest["stages"]:
        previous = saved.get(item["name"], {})
        states.append(
            SpectraStageState(
                name=item["name"],
                path=str(root_path / item["path"]),
                status=previous.get("status", "pending"),
                started_at=previous.get("started_at"),
                finished_at=previous.get("finished_at"),
                output=previous.get("output"),
                error=previous.get("error"),
                gap_eV=previous.get("gap_eV"),
            )
        )
    return states


def format_calculator_spectra_status(states: Iterable[SpectraStageState]) -> str:
    rows = [
        (item.name, item.status, "" if item.gap_eV is None else f"{item.gap_eV:.6f}", item.error or "")
        for item in states
    ]
    headers = ("stage", "status", "gap_eV", "error")
    widths = [max(len(headers[i]), *(len(row[i]) for row in rows)) for i in range(4)]
    lines = ["  ".join(headers[i].ljust(widths[i]) for i in range(4))]
    lines.append("  ".join("-" * width for width in widths))
    lines.extend("  ".join(row[i].ljust(widths[i]) for i in range(4)) for row in rows)
    return "\n".join(lines)


def generate_calculator_spectra_script(
    root: str | Path,
    *,
    backend: str = "shell",
    output: str | Path | None = None,
    job_name: str = "zstar-spectra",
    nodes: int = 1,
    tasks: int = 1,
    cpus_per_task: int = 1,
    walltime: str = "24:00:00",
    queue: str | None = None,
    account: str | None = None,
    env_script: str | Path | None = None,
    calculator_command: str | None = None,
    min_gap_eV: float = 0.01,
) -> Path:
    """Generate one shell, Slurm, or Torque driver for the serial workflow."""

    backend_key = backend.lower()
    if backend_key not in {"shell", "slurm", "torque"}:
        raise ValueError("backend must be shell, slurm, or torque")
    if min(nodes, tasks, cpus_per_task) < 1:
        raise ValueError("nodes, tasks, and cpus_per_task must be positive")
    root_path, manifest = _load_manifest(root)
    calculator = manifest["calculator"]
    if calculator_command is None:
        launcher = f"srun --ntasks={tasks}" if backend_key == "slurm" else f"mpirun -np {tasks}"
        calculator_command = (
            f"{launcher} vasp_std"
            if calculator == "vasp"
            else f"{launcher} cp2k.psmp -i input.inp -o output.log"
        )
    if output is None:
        suffix = {"shell": "sh", "slurm": "slurm", "torque": "pbs"}[backend_key]
        target = root_path / f"run_spectra.{suffix}"
    else:
        target = Path(output).resolve()
    header = [
        "#!/usr/bin/env bash",
        f"# ZStar {calculator.upper()} IR/Raman backend: {backend_key}",
        "# One serial, resumable driver for all spectroscopy stages.",
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
            f"zstar spectra run --root \"$ROOT\" --command {shlex.quote(calculator_command)} "
            f"--min-gap {min_gap_eV:.10g} --omp-threads {cpus_per_task} "
            '2>&1 | tee -a "$ROOT/.zstar/workflow.log"',
            'zstar spectra collect --root "$ROOT" 2>&1 | tee -a "$ROOT/.zstar/workflow.log"',
        ]
    )
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(header + [""] + body) + "\n", encoding="utf-8", newline="\n")
    target.chmod(target.stat().st_mode | 0o111)
    return target


def _collect_vasp_spectra(
    root: Path,
    manifest: dict,
    *,
    broadening_cm1: float,
    laser_nm: float,
    temperature_K: float,
    points: int,
    plot: bool,
    imaginary_tolerance_cm1: float,
    allow_imaginary: bool,
) -> dict:
    modes = load_vasp_gamma_modes(manifest["modes_source"])
    epsilon, tensors = parse_vasp_outcar(root / "reference" / "OUTCAR")
    selected = [int(entry["mode"]) for entry in manifest["modes"]]
    dimensionality = int(manifest["dimensionality"])
    born = BornData(tensors=tensors, electronic_dielectric=epsilon, source="VASP OUTCAR")
    if dimensionality in {1, 3}:
        ir = calculate_ir_spectrum(
            modes,
            born,
            dimensionality=dimensionality,
            mode_numbers=selected,
            broadening_cm1=broadening_cm1,
            points=points,
            periodic_axis=2,
            imaginary_tolerance_cm1=imaginary_tolerance_cm1,
            allow_imaginary=allow_imaginary,
        )
        ir_summary = write_ir_outputs(root / "ir_spectrum", ir, plot=plot)
    else:
        derivatives = mode_effective_charges(modes, tensors)[np.asarray(selected) - 1]
        derivatives *= DEBYE_PER_E_ANGSTROM
        ir = calculate_molecular_ir_spectrum(
            modes,
            selected,
            derivatives,
            broadening_cm1=broadening_cm1,
            points=points,
            imaginary_tolerance_cm1=imaginary_tolerance_cm1,
            allow_imaginary=allow_imaginary,
        )
        ir_summary = write_molecular_ir_outputs(root / "ir_spectrum", ir, plot=plot)

    raman_tensors: list[np.ndarray] = []
    sources: list[dict] = []
    amplitude = float(manifest["amplitude_A_sqrt_amu"])
    for entry in manifest["modes"]:
        plus, _plus_bec = parse_vasp_outcar(root / entry["plus"] / "OUTCAR")
        minus, _minus_bec = parse_vasp_outcar(root / entry["minus"] / "OUTCAR")
        derivative = (plus - minus) / (2.0 * amplitude)
        if dimensionality == 0:
            derivative *= modes.volume_angstrom3 / (4.0 * math.pi)
        elif dimensionality == 1:
            derivative *= modes.cross_section_angstrom2(2) / (4.0 * math.pi)
        raman_tensors.append(0.5 * (derivative + derivative.T))
        sources.append({"mode": entry["mode"], "plus": entry["plus"], "minus": entry["minus"]})
    tensor_kind = {
        0: "molecular polarizability derivative (Angstrom^3 per Angstrom sqrt(amu))",
        1: "1D line polarizability derivative (Angstrom^2 per Angstrom sqrt(amu))",
        3: "dielectric tensor derivative",
    }[dimensionality]
    raman = calculate_raman_spectrum(
        modes,
        selected,
        np.asarray(raman_tensors),
        tensor_kind=tensor_kind,
        temperature_K=temperature_K,
        laser_nm=laser_nm,
        broadening_cm1=broadening_cm1,
        points=points,
        imaginary_tolerance_cm1=imaginary_tolerance_cm1,
        allow_imaginary=allow_imaginary,
    )
    raman_summary = write_raman_outputs(root / "raman_spectrum", raman, plot=plot)
    result = {
        "schema_version": 1,
        "calculator": "vasp",
        "dimensionality": dimensionality,
        "periodic_axes": manifest.get("periodic_axes"),
        "nac_model": manifest.get("nac_model"),
        "mode_numbers": selected,
        "frequencies_cm-1": modes.frequencies_cm1[np.asarray(selected) - 1].tolist(),
        "electronic_dielectric": epsilon.tolist(),
        "born_tensors": tensors.tolist(),
        "raman_tensors": np.asarray(raman_tensors).tolist(),
        "raman_sources": sources,
        "ir_summary": ir_summary,
        "raman_summary": raman_summary,
    }
    (root / "spectra_results.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8", newline="\n"
    )
    return result


def _ensure_cp2k_moments(text: str, periodic: bool) -> str:
    path = ("FORCE_EVAL", "DFT", "PRINT", "MOMENTS")
    value = "TRUE" if periodic else "FALSE"
    lines = text.splitlines()
    if _has_section(lines, path):
        text = _set_section_keyword(text, path, "PERIODIC", value)
        if not periodic:
            text = _set_section_keyword(text, path, "REFERENCE", "COM")
        return text
    block = ["      &MOMENTS", f"        PERIODIC {value}"]
    if not periodic:
        block.append("        REFERENCE COM")
    block.append("      &END MOMENTS")
    print_path = ("FORCE_EVAL", "DFT", "PRINT")
    if _has_section(lines, print_path):
        return _insert_before_section_end(text, print_path, block)
    return _insert_before_section_end(
        text, ("FORCE_EVAL", "DFT"), ["    &PRINT", *block, "    &END PRINT"]
    )


def _ensure_cp2k_centered_molecule(text: str) -> str:
    """Center isolated coordinates so their density decays at every box face."""

    lines = text.splitlines()
    topology = ("FORCE_EVAL", "SUBSYS", "TOPOLOGY")
    centered = (*topology, "CENTER_COORDINATES")
    if _has_section(lines, centered):
        return text
    block = ["      &CENTER_COORDINATES TRUE", "      &END CENTER_COORDINATES"]
    if _has_section(lines, topology):
        return _insert_before_section_end(text, topology, block)
    return _insert_before_section_end(
        text,
        ("FORCE_EVAL", "SUBSYS"),
        ["    &TOPOLOGY", *block, "    &END TOPOLOGY"],
    )


def _ensure_cp2k_polar(text: str, periodic: bool) -> str:
    lines = text.splitlines()
    properties = ("FORCE_EVAL", "PROPERTIES")
    linres = (*properties, "LINRES")
    polar = (*linres, "POLAR")
    polar_block = [
        "      &POLAR TRUE",
        "        DO_RAMAN TRUE",
        f"        PERIODIC_DIPOLE_OPERATOR {'TRUE' if periodic else 'FALSE'}",
        "      &END POLAR",
    ]
    linres_block = [
        "    &LINRES",
        "      EPS 1.0E-6",
        "      MAX_ITER 50",
        "      PRECONDITIONER FULL_SINGLE_INVERSE",
        *polar_block,
        "    &END LINRES",
    ]
    if not _has_section(lines, properties):
        text = _insert_before_section_end(
            text,
            ("FORCE_EVAL",),
            ["  &PROPERTIES", *linres_block, "  &END PROPERTIES"],
        )
    elif not _has_section(text.splitlines(), linres):
        text = _insert_before_section_end(text, properties, linres_block)
    elif not _has_section(text.splitlines(), polar):
        text = _insert_before_section_end(text, linres, polar_block)
    else:
        text = _set_section_keyword(text, polar, "DO_RAMAN", "TRUE")
        text = _set_section_keyword(
            text,
            polar,
            "PERIODIC_DIPOLE_OPERATOR",
            "TRUE" if periodic else "FALSE",
        )
    text = _set_section_keyword(text, linres, "EPS", "1.0E-6")
    text = _set_section_keyword(text, linres, "MAX_ITER", "50")
    text = _set_section_keyword(
        text, linres, "PRECONDITIONER", "FULL_SINGLE_INVERSE"
    )
    return text


def prepare_cp2k_spectra(
    input_path: str | Path,
    root: str | Path = "cp2k_spectra",
    *,
    displacement_bohr: float = 0.01,
    dimensionality: int = 0,
    imaginary_tolerance_cm1: float = 20.0,
    allow_imaginary: bool = False,
    force: bool = False,
) -> Path:
    """Prepare CP2K native vibrational analysis with IR and Raman intensities."""

    if dimensionality == 2:
        raise ValueError(
            "CP2K 2D out-of-plane spectroscopy is not enabled; use the ZStar "
            "real-space slab polarization workflow"
        )
    if dimensionality not in {0, 3}:
        raise ValueError("CP2K spectroscopy dimensionality must be 0 or 3")
    if displacement_bohr <= 0.0:
        raise ValueError("displacement_bohr must be positive")
    source = Path(input_path).resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    text = source.read_text(encoding="utf-8")
    text = _set_section_keyword(text, ("GLOBAL",), "RUN_TYPE", "VIBRATIONAL_ANALYSIS")
    vib = ("VIBRATIONAL_ANALYSIS",)
    if _has_section(text.splitlines(), vib):
        text = _set_section_keyword(text, vib, "DX", f"{displacement_bohr:.10g}")
        text = _set_section_keyword(text, vib, "INTENSITIES", "TRUE")
        if dimensionality == 3:
            text = _set_section_keyword(text, vib, "FULLY_PERIODIC", "TRUE")
    else:
        text += (
            "\n&VIBRATIONAL_ANALYSIS\n"
            f"  DX {displacement_bohr:.10g}\n"
            "  INTENSITIES TRUE\n"
            + ("  FULLY_PERIODIC TRUE\n" if dimensionality == 3 else "")
            + "&END VIBRATIONAL_ANALYSIS\n"
        )
    text = _ensure_cp2k_moments(text, periodic=dimensionality == 3)
    text = _ensure_cp2k_polar(text, periodic=dimensionality == 3)
    if dimensionality == 0:
        text = _ensure_cp2k_centered_molecule(text)

    target = Path(root).resolve()
    if target.exists() and force:
        shutil.rmtree(target)
    if target.exists() and any(target.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {target}")
    calculation = target / "calculation"
    calculation.mkdir(parents=True, exist_ok=True)
    (calculation / "input.inp").write_text(text, encoding="utf-8", newline="\n")
    for asset in _referenced_assets(text, source.parent):
        shutil.copy2(asset, calculation / asset.name)
    manifest = {
        "schema_version": 1,
        "calculator": "cp2k",
        "created_at": _utc_now(),
        "source_input": str(source),
        "dimensionality": dimensionality,
        "displacement_bohr": displacement_bohr,
        "imaginary_tolerance_cm-1": float(imaginary_tolerance_cm1),
        "allow_imaginary": bool(allow_imaginary),
        "stages": [{"name": "calculation", "path": "calculation"}],
    }
    (target / "spectra_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8", newline="\n"
    )
    return target


def parse_cp2k_native_spectra(path: str | Path) -> dict:
    """Parse CP2K's final native vibrational frequency and intensity blocks."""

    source = Path(path)
    frequencies: list[float] = []
    ir: list[float] = []
    raman: list[float] = []
    text = source.read_text(encoding="utf-8", errors="ignore")
    if re.search(r"\b(?:NAN|INF)\b", text, re.IGNORECASE):
        raise ValueError(f"CP2K spectroscopy output contains non-finite values: {source}")
    for line in text.splitlines():
        upper = line.upper()
        if "VIB|FREQUENCY (CM^-1)" in upper:
            frequencies.extend(_float(value) for value in re.findall(_NUMBER, line.split(")", 1)[-1]))
        elif "VIB|IR INT (KM/MOLE)" in upper:
            ir.extend(_float(value) for value in re.findall(_NUMBER, line.split(")", 1)[-1]))
        elif "VIB|RAMAN (A^4/AMU)" in upper:
            raman.extend(_float(value) for value in re.findall(_NUMBER, line.split(")", 1)[-1]))
    if not frequencies:
        raise ValueError(f"No CP2K vibrational frequencies found in {source}")
    if len(ir) != len(frequencies):
        raise ValueError(
            f"CP2K IR intensity count {len(ir)} does not match {len(frequencies)} frequencies"
        )
    if len(raman) != len(frequencies):
        raise ValueError(
            f"CP2K Raman intensity count {len(raman)} does not match {len(frequencies)} frequencies"
        )
    return {
        "mode_numbers": list(range(1, len(frequencies) + 1)),
        "frequencies_cm-1": frequencies,
        "ir_intensities_km_mol": ir,
        "raman_activities_A4_amu": raman,
    }


def _collect_cp2k_spectra(
    root: Path,
    manifest: dict,
    *,
    broadening_cm1: float,
    points: int,
    plot: bool,
    imaginary_tolerance_cm1: float,
    allow_imaginary: bool,
) -> dict:
    parsed = parse_cp2k_native_spectra(root / "calculation" / "output.log")
    validate_frequencies_stable(
        parsed["frequencies_cm-1"],
        imaginary_tolerance_cm1=imaginary_tolerance_cm1,
        allow_imaginary=allow_imaginary,
    )
    common = {
        "mode_numbers": parsed["mode_numbers"],
        "broadening_cm1": broadening_cm1,
        "points": points,
    }
    ir = calculate_native_line_spectrum(
        parsed["frequencies_cm-1"],
        parsed["ir_intensities_km_mol"],
        activity_kind="IR_intensity",
        activity_unit="km/mol",
        **common,
    )
    raman = calculate_native_line_spectrum(
        parsed["frequencies_cm-1"],
        parsed["raman_activities_A4_amu"],
        activity_kind="Raman_activity",
        activity_unit="Angstrom^4/amu",
        **common,
    )
    parsed.update(
        {
            "schema_version": 1,
            "calculator": "cp2k",
            "dimensionality": manifest["dimensionality"],
            "ir_summary": write_native_line_spectrum_outputs(
                root / "ir_spectrum", ir, stem="ir", plot=plot
            ),
            "raman_summary": write_native_line_spectrum_outputs(
                root / "raman_spectrum", raman, stem="raman", plot=plot
            ),
        }
    )
    (root / "spectra_results.json").write_text(
        json.dumps(parsed, indent=2), encoding="utf-8", newline="\n"
    )
    return parsed


def collect_calculator_spectra(
    root: str | Path,
    *,
    broadening_cm1: float = 8.0,
    laser_nm: float = 532.0,
    temperature_K: float = 300.0,
    points: int = 2001,
    plot: bool = True,
    imaginary_tolerance_cm1: float | None = None,
    allow_imaginary: bool | None = None,
) -> dict:
    root_path, manifest = _load_manifest(root)
    if imaginary_tolerance_cm1 is None:
        imaginary_tolerance_cm1 = float(
            manifest.get("imaginary_tolerance_cm-1", 20.0)
        )
    if allow_imaginary is None:
        allow_imaginary = bool(manifest.get("allow_imaginary", False))
    if manifest["calculator"] == "vasp":
        return _collect_vasp_spectra(
            root_path,
            manifest,
            broadening_cm1=broadening_cm1,
            laser_nm=laser_nm,
            temperature_K=temperature_K,
            points=points,
            plot=plot,
            imaginary_tolerance_cm1=imaginary_tolerance_cm1,
            allow_imaginary=allow_imaginary,
        )
    if manifest["calculator"] == "cp2k":
        return _collect_cp2k_spectra(
            root_path,
            manifest,
            broadening_cm1=broadening_cm1,
            points=points,
            plot=plot,
            imaginary_tolerance_cm1=imaginary_tolerance_cm1,
            allow_imaginary=allow_imaginary,
        )
    raise ValueError(f"Unsupported calculator: {manifest['calculator']}")
