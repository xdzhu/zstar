"""Matched cubic-BTO unified versus separate BEC/phonon experiment.

This research driver does not change the public CLI. All calculations are
new; the source archive is only read. Timings use a monotonic clock.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import time

import numpy as np
import spglib

from zstar import workflow
from zstar.shared_abacus import prepare_shared_abacus, collect_shared_abacus
from zstar.shared_response import DEFAULT_DISTANCE, actual_displacement, make_phonopy, read_structure, write_structure


@contextmanager
def working_directory(path):
    old = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def prepare(root, source):
    from zstar.gen_polar import gen_polar

    if (root / "plan.json").exists():
        raise FileExistsError("Experiment already prepared; use run to resume")
    seed = root / "seed"
    shutil.copytree(source, seed)
    # Local files override any unrelated environment-wide PP/ORB directories.
    for path in (seed / "assets").iterdir():
        shutil.copy2(path, seed / path.name)
    text = (seed / "INPUT").read_text()
    for name, value in {"kspacing": "0", "symmetry": "1", "cal_force": "0",
                        "dft_functional": "pbesol", "scf_thr": "1e-7"}.items():
        text = workflow._set_abacus_parameter(text, name, value)
    (seed / "INPUT").write_text(text)
    atoms = read_structure(seed / "STRU")
    ds = spglib.get_symmetry_dataset(atoms.totuple(), symprec=1e-5)
    if ds.number != 221 or atoms.symbols != ["Ba", "Ti", "O", "O", "O"]:
        raise ValueError("Expected the five-atom cubic Pm-3m BTO archive")
    unified = root / "unified"
    metadata = prepare_shared_abacus(seed / "STRU", root=unified,
        scf_input=seed / "INPUT", dimension=3, symprec=1e-5,
        method="auto", displacement_angstrom=DEFAULT_DISTANCE)

    legacy = root / "legacy_bec"
    legacy.mkdir()
    shutil.copy2(seed / "STRU", legacy / "STRU")
    with working_directory(legacy):
        gen_polar("STRU", scf_input=seed / "INPUT", xc=None, dimension=3,
                  k_grid=0, input_sets=[str(p) for p in (seed / "assets").iterdir()],
                  method="forward", displacement_angstrom=DEFAULT_DISTANCE,
                  extract_starred_atoms_only=True, force_delete=False)
    old_stages = workflow.discover_stages(legacy)
    for stage in old_stages:
        shutil.copy2(seed / "KPT", stage.path / "KPT")
    entries = []
    ref = read_structure(legacy / "0.no-move/STRU")
    for stage in old_stages[1:]:
        atom, vector = actual_displacement(ref, read_structure(stage.path / "STRU"))
        entries.append({"name": stage.name, "atom": atom,
                        "displacement_A": vector.tolist(),
                        "distance_A": float(np.linalg.norm(vector))})
    write_json(legacy / "displacements.json", entries)

    # Old phonons are separate force-only SCFs, not reused BEC force outputs.
    phonon_dir = root / "legacy_phonon"
    phonon_dir.mkdir()
    phonon = make_phonopy(atoms, symprec=1e-5)
    phonon.generate_displacements(distance=DEFAULT_DISTANCE, is_plusminus="auto")
    phonon.save(filename=str(phonon_dir / "phonopy_disp.yaml"))
    shutil.copy2(seed / "STRU", phonon_dir / "STRU")
    force_input = text
    for name, value in {"cal_force": "1", "init_chg": "auto", "out_chg": "0",
                        "out_mat_hs2": "0", "out_mat_r": "0"}.items():
        force_input = workflow._set_abacus_parameter(force_input, name, value)
    force_entries = []
    for i, cell in enumerate(phonon.supercells_with_displacements, 1):
        stage = phonon_dir / f"disp-{i:03d}"
        stage.mkdir()
        write_structure(seed / "STRU", stage / "STRU", cell)
        (stage / "INPUT").write_text(force_input)
        shutil.copy2(seed / "KPT", stage / "KPT")
        for asset in (seed / "assets").iterdir():
            shutil.copy2(asset, stage / asset.name)
        atom, vector = actual_displacement(atoms, read_structure(stage / "STRU"))
        force_entries.append({"name": stage.name, "atom": atom,
                              "displacement_A": vector.tolist(),
                              "distance_A": float(np.linalg.norm(vector))})
    write_json(phonon_dir / "displacements.json", force_entries)
    hashes = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
              for sub in (seed, unified, legacy, phonon_dir)
              for p in sub.rglob("*") if p.is_file()}
    plan = {"system": "cubic BaTiO3", "space_group": ds.international,
            "source": str(source), "cell_A": atoms.cell.tolist(),
            "xc": "PBEsol", "ecutwfc_Ry": 100, "scf_thr": 1e-7,
            "DFT_mesh": [9, 9, 9], "PYATB_mp_density": 0.08,
            "nominal_distance_A": DEFAULT_DISTANCE,
            "unified_displacements": len(metadata["stages"]),
            "legacy_BEC_displacements": len(entries),
            "legacy_phonon_displacements": len(force_entries),
            "legacy_BEC_scheme": "symmetry-inequivalent atoms, Cartesian forward",
            "legacy_phonon_start": "independent atomic-charge SCFs, no matrix/density export",
            "unified_and_legacy_BEC_start": "reference first, private reference-charge copies",
            "precision": "same output-only full-precision PYATB adapter on both BEC routes",
            "timing": "completed command monotonic seconds times 40 reserved physical cores; components reported separately",
            "exclusions": "no relaxation, Raman, or finite-q dispersion; no stable static dielectric claim if soft modes remain",
            "input_sha256": hashes}
    write_json(root / "plan.json", plan)
    print(json.dumps({k: v for k, v in plan.items() if k != "input_sha256"}, indent=2))


def run(root, route, binary):
    from zstar.pyatb_precision import precision_command
    from zstar.shared_abacus import read_forces

    output = root / route
    if not (root / "plan.json").exists():
        raise ValueError("Prepare and audit the experiment first")
    lock = output / ".worker.lock"
    with lock.open("x") as handle:
        handle.write(f"{socket.gethostname()} {os.getpid()}\n")
    records = output / "component_times.jsonl"
    original = workflow._run_shell
    abacus = f"mpirun -np 1 {binary}"
    bindir = Path(sys.executable).parent

    def timed(command, **kwargs):
        is_pyatb = "-m zstar.pyatb_precision" in command
        environment = dict(kwargs["env"])
        environment.update(OMP_NUM_THREADS="1" if is_pyatb else "40",
                           MKL_NUM_THREADS="1" if is_pyatb else "40",
                           OPENBLAS_NUM_THREADS="1", I_MPI_PIN_DOMAIN="omp")
        kwargs["env"] = environment
        record = {"command": command, "cwd": str(kwargs["cwd"]),
                  "host": socket.gethostname(), "mpi": 40 if is_pyatb else 1,
                  "omp": 1 if is_pyatb else 40, "reserved_cores": 40,
                  "kind": "SCF" if command == abacus else
                          "input_preparation" if "pyatb_input" in command else
                          "band_gate" if str(kwargs["cwd"]).endswith("pyatb-band") else
                          "polarization_and_electronic_response"}
        start = time.monotonic()
        try:
            value = original(command, **kwargs)
            record["success"] = True
            return value
        except Exception as exc:
            record.update(success=False, error=repr(exc))
            raise
        finally:
            record["wall_seconds"] = time.monotonic() - start
            record["reserved_core_hours"] = record["wall_seconds"] * 40 / 3600
            with records.open("a") as handle:
                handle.write(json.dumps(record) + "\n")
            print(json.dumps(record), flush=True)

    workflow._run_shell = timed
    worker = {"host": socket.gethostname(), "pid": os.getpid(),
               "python": sys.executable, "abacus": binary,
               "cpu": subprocess.run(["lscpu"], capture_output=True, text=True).stdout,
               "status": "running"}
    write_json(output / "worker.json", worker)
    try:
        if route == "legacy_phonon":
            for item in json.loads((output / "displacements.json").read_text()):
                stage = output / item["name"]
                if not workflow.scf_is_complete(stage):
                    timed(abacus, cwd=stage, env=dict(os.environ),
                          log_path=stage / "driver.log", dry_run=False)
                read_forces(stage)
        else:
            workflow.run_serial_workflow(output, abacus_command=abacus,
                pyatb_input=str(bindir / "pyatb_input"), pyatb_executable=str(bindir / "pyatb"),
                pyatb_command=precision_command(f"mpirun -np 40 {bindir / 'pyatb'}"),
                omp_threads=40, dimensionality=3, mp_density=0.08)
            if route == "unified":
                collect_shared_abacus(output)
        write_json(output / "completed.json", {"host": socket.gethostname(), "status": "completed"})
        worker["status"] = "completed"
    except Exception as exc:
        worker.update(status="failed", error=repr(exc))
        raise
    finally:
        write_json(output / "worker.json", worker)
        workflow._run_shell = original
        lock.unlink()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["prepare", "run"])
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--route", choices=["unified", "legacy_bec", "legacy_phonon"])
    parser.add_argument("--abacus", default="/home/zhuxd/Software/abacus/INSTALL/3.10.0-LTS/bin/abacus")
    args = parser.parse_args()
    if args.action == "prepare":
        prepare(args.root.resolve(), args.source.resolve())
    else:
        run(args.root.resolve(), args.route, args.abacus)
