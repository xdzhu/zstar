"""Audit and export completed BTO runs without inventing missing results."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from phonopy.file_IO import write_FORCE_CONSTANTS, write_FORCE_SETS

from zstar import workflow
from zstar.pyatb_compat import read_static_dielectric
from zstar.shared_abacus import _dipole_changes, load_manifest, read_forces
from zstar.shared_response import (SharedResponse, actual_displacement,
    make_phonopy, project_response, read_structure, reconstruct_responses, symmetry_operations)


def load(path):
    return json.loads(path.read_text())


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def ledger(route, expected_scf):
    if not (route / "completed.json").is_file():
        raise ValueError(f"Not completed: {route}")
    items = [json.loads(line) for line in (route / "component_times.jsonl").read_text().splitlines()]
    if not all(item["success"] for item in items):
        raise ValueError(f"Failed commands require separate accounting: {route}")
    if sum(item["kind"] == "SCF" for item in items) != expected_scf:
        raise ValueError(f"Unexpected SCF count, possibly restarts: {route}")
    groups = {}
    for item in items:
        group = groups.setdefault(item["kind"], {"calls": 0, "wall_seconds": 0., "core_hours": 0.})
        group["calls"] += 1
        group["wall_seconds"] += item["wall_seconds"]
        group["core_hours"] += item["reserved_core_hours"]
    worker = load(route / "worker.json")
    # Earlier ledgers retained worker.json as an immutable launch snapshot.
    worker["status"] = load(route / "completed.json")["status"]
    return {"components": groups, "core_hours": sum(g["core_hours"] for g in groups.values()),
            "hosts": sorted({item["host"] for item in items}),
            "worker": worker}


def phonopy_from_observations(atoms, entries, force_reference):
    p = make_phonopy(atoms, symprec=1e-5)
    p.dataset = {"natom": len(atoms), "first_atoms": [
        {"number": s["atom"], "displacement": np.array(s["displacement_A"]),
         "forces": np.array(s["forces_eV_A"]) - force_reference} for s in entries]}
    p.produce_force_constants(fc_calculator="traditional")
    return p


def direct_cartesian_born(root, entries, atoms):
    dipoles = _dipole_changes(root, {"dimension": 3, "stages": entries})
    measured = {}
    for item, dipole in zip(entries, dipoles):
        u = np.array(item["displacement_A"])
        axis = np.flatnonzero(np.abs(u) > 1e-8)
        if len(axis) != 1:
            raise ValueError("Legacy baseline is not Cartesian")
        tensor = measured.setdefault(item["atom"], np.full((3, 3), np.nan))
        tensor[:, axis[0]] = dipole / u[axis[0]]
    if not all(np.isfinite(t).all() for t in measured.values()):
        raise ValueError("Incomplete Cartesian response")
    # No unified least-squares fit is used for the reference BEC calculation.
    expanded = np.zeros((len(atoms), 3, 3))
    counts = np.zeros(len(atoms))
    for atom, tensor in measured.items():
        for rotation, permutation in symmetry_operations(make_phonopy(atoms)):
            expanded[permutation[atom]] += rotation @ tensor @ rotation.T
            counts[permutation[atom]] += 1
    if np.any(counts == 0):
        raise ValueError("Uncovered atomic orbit")
    return expanded / counts[:, None, None], {str(k): v.tolist() for k, v in measured.items()}


def report(root, *, quiet=False):
    plan = load(root / "plan.json")
    for name, expected in plan["input_sha256"].items():
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Prepared input changed: {name}")
    unified, old_bec, old_ph = (root / n for n in ("unified", "legacy_bec", "legacy_phonon"))
    costs = {name: ledger(root / name, count) for name, count in
             (("unified", 4), ("legacy_bec", 10), ("legacy_phonon", 3))}
    result = load(unified / "shared_response_result.json")
    atoms = read_structure(unified / "0.no-move/STRU")
    for path in (old_bec / "0.no-move/STRU", old_ph / "STRU"):
        other = read_structure(path)
        np.testing.assert_allclose(other.cell, atoms.cell, atol=1e-10, rtol=0)
        np.testing.assert_allclose(other.positions, atoms.positions, atol=1e-10, rtol=0)
    for route, entries in ((old_bec, load(old_bec / "displacements.json")),
                           (old_ph, load(old_ph / "displacements.json"))):
        for item in entries:
            index, u = actual_displacement(atoms, read_structure(route / item["name"] / "STRU"))
            assert index == item["atom"]
            np.testing.assert_allclose(u, item["displacement_A"], atol=1e-10, rtol=0)
            if not workflow.scf_is_complete(route / item["name"]):
                raise ValueError("Incomplete electronic SCF")

    zold, direct = direct_cartesian_born(old_bec, load(old_bec / "displacements.json"), atoms)
    fentries = load(old_ph / "displacements.json")
    for item in fentries:
        item["forces_eV_A"] = read_forces(old_ph / item["name"]).tolist()
    pold = phonopy_from_observations(atoms, fentries, np.zeros((len(atoms), 3)))
    fc_old = pold.force_constants.transpose(1, 0, 3, 2).copy()
    write_FORCE_SETS(pold.dataset, filename=str(old_ph / "FORCE_SETS"))
    write_FORCE_CONSTANTS(pold.force_constants, filename=str(old_ph / "FORCE_CONSTANTS.raw"))
    projected_old = project_response(SharedResponse(zold, fc_old, {}))
    write_FORCE_CONSTANTS(projected_old.force_constants, filename=str(old_ph / "FORCE_CONSTANTS"))
    pold.force_constants = projected_old.force_constants
    pold.run_qpoints([[0, 0, 0]], with_eigenvectors=True)
    pold.save(filename=str(old_ph / "phonopy.yaml"), settings={"force_constants": True})
    f_old = pold.qpoints.frequencies[0] * 33.3564095198152
    raw = reconstruct_responses(len(atoms), result["observations"], symmetry_operations(make_phonopy(atoms)),
                                reference_forces=result["reference_forces_eV_A"])
    independent = phonopy_from_observations(atoms, result["observations"], result["reference_forces_eV_A"])
    independent_error = float(np.max(np.abs(independent.force_constants.transpose(1, 0, 3, 2) - raw.force_constants)))
    f_new = np.array(result["frequencies_THz"]) * 33.3564095198152
    epsilon_new, _ = read_static_dielectric(unified / "0.no-move/pyatb")
    epsilon_old, _ = read_static_dielectric(old_bec / "0.no-move/pyatb")
    precision_count = {}
    for name, route in (("unified", unified), ("legacy_bec", old_bec)):
        stages = workflow.discover_stages(route)
        precision_count[name] = sum((s.path / "pyatb/Out/Polarization/zstar_precision.json").is_file() for s in stages)
        if precision_count[name] != len(stages):
            raise ValueError("Missing full-precision polarization record")
    total_old = costs["legacy_bec"]["core_hours"] + costs["legacy_phonon"]["core_hours"]
    total_new = costs["unified"]["core_hours"]
    output = {"status": "completed_and_compared", "system": plan["system"],
        "xc": plan["xc"], "cell_A": atoms.cell.tolist(), "timing": costs,
        "counts": {"legacy_BEC_displacements": 9, "legacy_phonon_displacements": 3,
                   "unified_displacements": 3, "legacy_total_SCFs": 13, "unified_total_SCFs": 4},
        "legacy_total_core_hours": total_old, "unified_total_core_hours": total_new,
        "core_hours_saved": total_old - total_new, "core_hours_reduction_percent": 100 * (1 - total_new / total_old),
        "speedup": total_old / total_new,
        "legacy_born_cartesian_representatives_e": direct,
        "legacy_born_raw_e": zold.tolist(), "legacy_born_projected_e": projected_old.born.tolist(),
        "unified_born_raw_e": raw.born.tolist(), "unified_born_projected_e": result["born_projected_e"],
        "max_raw_BEC_difference_e": float(np.max(np.abs(raw.born - zold))),
        "max_projected_BEC_difference_e": float(np.max(np.abs(np.array(result["born_projected_e"]) - projected_old.born))),
        "max_raw_force_constants_difference_eV_A2": float(np.max(np.abs(raw.force_constants - fc_old))),
        "relative_raw_force_constants_difference": float(np.linalg.norm(raw.force_constants - fc_old) / np.linalg.norm(fc_old)),
        "independent_Phonopy_unified_raw_difference_eV_A2": independent_error,
        "legacy_frequencies_cm1": f_old.tolist(), "unified_frequencies_cm1": f_new.tolist(),
        "max_frequency_difference_cm1": float(np.max(np.abs(f_new - f_old))),
        "epsilon_infinity_legacy": epsilon_old.tolist(), "epsilon_infinity_unified": epsilon_new.tolist(),
        "static_phonon_response_status": "not_validated_unstable_reference" if min(f_old.min(), f_new.min()) < -1 else "requires_static_closure_check",
        "diagnostics": result["diagnostics"], "projected_diagnostics": result["projected_diagnostics"],
        "full_precision_polarization_stages": precision_count,
        "insulation": {"unified": load(unified / "0.no-move/zstar_insulation.json"),
                       "legacy": load(old_bec / "0.no-move/zstar_insulation.json")}}
    dump(root / "benchmark_summary.json", output)
    if not quiet:
        print(json.dumps({k: output[k] for k in (
            "status", "system", "counts", "legacy_total_core_hours", "unified_total_core_hours",
            "core_hours_reduction_percent", "speedup", "max_raw_BEC_difference_e",
            "max_frequency_difference_cm1", "static_phonon_response_status")}, indent=2))
    return output


def export(root, dest):
    def copy(path, target):
        if path.is_file():
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)

    if dest.exists():
        raise FileExistsError("Export to a new directory")
    dest.mkdir(parents=True)
    for name in ("INPUT", "STRU", "KPT"):
        copy(root / "seed" / name, dest / "run" / name)
    shutil.copytree(root / "seed/assets", dest / "run/assets")
    for name in ("plan.json", "benchmark_summary.json"):
        copy(root / name, dest / "results" / name)
    # Keep every prepared input, including resolver sidecars not required by
    # numerical post-processing but needed for a complete provenance audit.
    for name in load(root / "plan.json")["input_sha256"]:
        if Path(name).parts[0] != "seed":
            copy(root / name, dest / "results" / name)
    for route in ("unified", "legacy_bec", "legacy_phonon"):
        source, target = root / route, dest / "results" / route
        for path in source.iterdir():
            if path.is_file() and path.suffix != ".lock":
                copy(path, target / path.name)
        if route == "unified":
            meta = load_manifest(source)
            for name in meta["input_hashes"]:
                copy(source / name, target / name)
            names = ["0.no-move"] + [s["name"] for s in meta["stages"]]
        else:
            names = [s["name"] for s in load(source / "displacements.json")]
            if route == "legacy_bec":
                names.insert(0, "0.no-move")
        for name in names:
            stage = source / name
            for pattern in ("STRU", "INPUT", "INPUT-scf", "KPT", "*.upf", "*.orb",
                "OUT.*/running_scf.log", "OUT.*/kpoints", "OUT.*/INPUT", "OUT.*/warning.log",
                "pyatb/Input", "pyatb/Out/input.json", "pyatb/Out/Polarization/*",
                "pyatb/Out/Optical_Conductivity/*", "pyatb/Out/Optical_Conductivity/**/*",
                "pyatb-band/band_gap.json", "zstar_insulation.json"):
                for path in stage.glob(pattern):
                    copy(path, target / path.relative_to(source))
        for path in (source / ".zstar").rglob("*"):
            if path.suffix in (".json", ".jsonl"):
                copy(path, target / path.relative_to(source))
    dump(dest / "checksums.json", {p.relative_to(dest).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                                    for p in sorted(dest.rglob("*")) if p.is_file()})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--export", type=Path)
    args = parser.parse_args()
    report(args.root.resolve())
    if args.export:
        export(args.root.resolve(), args.export.resolve())
