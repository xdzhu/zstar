"""Read existing ABACUS response archives; never run DFT or alter source files."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
import platform
from datetime import datetime, timezone

import numpy as np
import phonopy
import spglib
from phonopy.interface.abacus import read_abacus_output

BOHR = 0.529177210903
E = 1.602176634e-19
ROOT = Path("/home/zhuxd/abacus")
OLD = ROOT / "agent-runs/20260723-zstar-validation"


def cases():
    benchmark = ROOT / "agent-runs/20260903-spectra-backend-benchmark/sic_abacus"
    hf = ROOT / "zstar_validation/hfo2_pbesol_tzdp9_20260901"
    yield "SiC", 3, benchmark / "workflow", benchmark / "phonon"
    yield "HfO2_t_TZDP", 3, hf / "bec", hf / "phonon_gamma"
    for name in ("bto_tet", "pto", "hfo2", "in2se3"):
        yield name + "_legacy", 2 if name == "in2se3" else 3, OLD / "cases" / name, OLD / "phonon-fresh" / name
    yield "hBN", 2, OLD / "cases/hbn", None
    yield "MoS2_D3BJ", 2, ROOT / "agent-runs/zstar_mos2_pbed3bj_response_20260831", ROOT / "agent-runs/zstar_mos2_pbed3bj_20260831"
    yield "In2Se3_PBEsol", 2, ROOT / "agent-runs/20260827-in2se3-pbesol-bec/bec", None
    for name in ("h2o", "ch4"):
        yield name.upper() + "_PBE", 0, ROOT / "agent-runs/20260826-zstar-molecular-apt" / name, None


def structure(path):
    from zstar.stru_analyzer import stru_analyzer

    a0, vectors, symbols, counts, mode, coords, _, _, masses, _, _ = stru_analyzer(str(path))
    lattice = np.asarray(vectors) * a0 * BOHR
    xyz = np.array([v for s in symbols for v in coords[s]])
    if mode.lower().startswith("cart"):
        xyz = xyz * a0 * BOHR @ np.linalg.inv(lattice)
    labels = [s for s in symbols for _ in range(counts[s])]
    mass = [float(m) for s, m in zip(symbols, masses) for _ in range(counts[s])]
    return {"cell_A": lattice.tolist(), "fractional_positions": xyz.tolist(), "symbols": labels, "masses_amu": mass}


def evidence(path):
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def read_electronic_dielectric(path):
    rows = [line.split("#", 1)[0].split() for line in path.read_text().splitlines()]
    rows = [row for row in rows if row]
    if len(rows[0]) != 9:
        rows = rows[1:]  # Optional Phonopy conversion-factor header.
    if not rows or len(rows[0]) != 9:
        raise ValueError(f"Missing dielectric tensor in {path}")
    return np.array([float(v) for v in rows[0]]).reshape(3, 3)


def inputs(stage):
    candidates = sorted(stage.glob("OUT.*/INPUT")) or [stage / "INPUT"]
    if not candidates[0].exists():
        return {}
    result = {}
    for line in candidates[0].read_text(errors="replace").splitlines():
        parts = line.split("#", 1)[0].split()
        if len(parts) > 1:
            result[parts[0]] = " ".join(parts[1:])
    keys = ("dft_functional", "ecutwfc", "scf_thr", "basis_type", "cal_force", "out_mat_hs2", "out_mat_r", "nspin", "noncolin", "lspinorb", "efield_flag", "vdw_method", "kspacing", "init_chg")
    return {key: result[key] for key in keys if key in result}


def forces(stage):
    logs = sorted(stage.glob("OUT.*/running_scf.log"))
    for path in logs:
        text = path.read_text(errors="replace")
        if "TOTAL-FORCE (eV/Angstrom)" in text:
            return np.asarray(read_abacus_output(str(path))).tolist(), evidence(path)
    return None, None


def collect_phonon(root):
    if root is None:
        return None
    disp, force_sets = root / "phonopy_disp.yaml", root / "FORCE_SETS"
    if not disp.exists() or not force_sets.exists():
        return {"root": str(root), "status": "missing_displacements_or_forces"}
    import yaml

    meta = yaml.safe_load(disp.read_text())
    unit = meta["physical_unit"]["length"]
    length = BOHR if unit.lower() in ("au", "bohr") else 1.0
    p = phonopy.load(str(disp), force_sets_filename=str(force_sets), is_nac=False, produce_fc=False)
    unitcell = p.unitcell
    data = [{"atom": int(d["number"]), "displacement_A": (d["displacement"] * length).tolist(), "forces_eV_A": d["forces"].tolist()} for d in p.dataset["first_atoms"]]
    settings_path = next(iter(sorted(root.glob("disp-*/OUT.*/INPUT"))), None)
    settings = inputs(settings_path.parent.parent) if settings_path else {}
    result = {"root": str(root), "status": "available", "native_length_unit": unit,
            "generation_symprec_A": float(meta["phonopy"].get("symmetry_tolerance", 1e-5)) * length,
            "structure": {"cell_A": (unitcell.cell * length).tolist(), "fractional_positions": unitcell.scaled_positions.tolist(), "symbols": unitcell.symbols, "masses_amu": unitcell.masses.tolist()},
            "supercell_matrix": p.supercell_matrix.tolist(), "stages": data, "settings": settings,
            "sources": [evidence(disp), evidence(force_sets)]}
    qpoints = root / "qpoints.yaml"
    if qpoints.exists():
        result["archived_qpoints"] = yaml.safe_load(qpoints.read_text())
        result["sources"].append(evidence(qpoints))
    return result


def collect_case(name, dim, root, phroot):
    from zstar.deal_polar import _parse_pyatb_polar_file, _read_pyatb_geom
    from zstar.spectra import read_pyatb_polarization

    reference = root / "0.no-move"
    geom = structure(reference / "STRU")
    cell = np.array(geom["cell_A"])
    pos = np.array(geom["fractional_positions"])
    polpath = Path("pyatb/Out/Polarization/polarization.dat")

    def read(stage):
        if dim == 0:
            p, q, _ = read_pyatb_polarization(stage / polpath)
            return p, q
        values = _parse_pyatb_polar_file(stage / polpath)
        return np.array(values[:3]), np.array(values[3:])

    p0, q0 = read(reference)
    transform, vol_m3 = _read_pyatb_geom(reference / "pyatb/Out/input.json")
    assert np.isclose(vol_m3 / 1e-30, abs(np.linalg.det(cell)), rtol=2e-6)
    f0, f0_source = forces(reference)
    result = {"name": name, "dimension": dim, "root": str(root), "structure": geom,
              "settings": inputs(reference), "reference_forces_eV_A": f0,
              "sources": [evidence(reference / "STRU"), evidence(reference / polpath)], "stages": [],
              "phonon": collect_phonon(phroot), "tensor_convention": "Z[polarization, displacement] in e"}
    if f0_source:
        result["sources"].append(f0_source)
    for atomdir in sorted(root.iterdir()):
        match = re.fullmatch(r"([1-9]\d*)\.([A-Za-z]+)", atomdir.name)
        if not match:
            continue
        atom = int(match[1]) - 1
        hybrid = None
        if dim == 2:
            hybrid_path = atomdir / "zstar_2d_bec.json"
            hybrid = json.loads(hybrid_path.read_text())
            result["sources"].append(evidence(hybrid_path))
        for stage in sorted(atomdir.iterdir()):
            if not re.fullmatch(r"[xyz][+-]?", stage.name) or not stage.is_dir():
                continue
            moved = structure(stage / "STRU")
            delta = np.array(moved["fractional_positions"]) - pos
            delta -= np.rint(delta)
            u = delta @ cell
            others = np.delete(u, atom, axis=0)
            assert np.max(np.abs(others), initial=0) < 1e-7, (name, stage, "multiatom displacement")
            p, q = read(stage)
            dp = p - p0
            dp -= np.rint(dp / q0) * q0
            response = dp @ transform * vol_m3 / E / 1e-10
            if hybrid:
                sign = "minus" if stage.name.endswith("-") else "plus"
                diag = hybrid["diagnostics"]
                ref = diag["reference_dipole"]
                value = diag["directions"][stage.name[0]][sign]
                dmu = value["dipole_e_bohr"] - ref["dipole_e_bohr"]
                dmu -= round(dmu / ref["height_bohr"]) * ref["height_bohr"]
                response[2] = dmu * BOHR
            force, source = forces(stage)
            item = {"name": str(stage.relative_to(root)), "atom": atom, "displacement_A": u[atom].tolist(),
                    "dipole_change_e_A": response.tolist(), "forces_eV_A": force,
                    "sources": [evidence(stage / "STRU"), evidence(stage / polpath)], "settings": inputs(stage)}
            if source:
                item["sources"].append(source)
            result["stages"].append(item)
    born = root / "BORN"
    if born.exists():
        result["electronic_dielectric"] = read_electronic_dielectric(born).tolist()
        result["sources"].append(evidence(born))
    stored = root / "Z-BORN-symm.out"
    if stored.exists():
        from zstar.molecular_bec import _read_zborn
        result["stored_symmetrized_BEC"] = _read_zborn(stored)
        result["sources"].append(evidence(stored))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--code-root", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, args.code_root)
    rows, errors = [], []
    for spec in cases():
        try:
            row = collect_case(*spec)
            rows.append(row)
            print(spec[0], len(row["stages"]), "stages", flush=True)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            errors.append({"case": spec[0], "error": str(exc)})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    metadata = {"collected_UTC": datetime.now(timezone.utc).isoformat(), "python": platform.python_version(),
                "numpy": np.__version__, "phonopy": phonopy.__version__, "spglib": spglib.__version__,
                "collector": evidence(Path(__file__)), "readers": [evidence(Path(args.code_root) / "zstar" / name) for name in ("stru_analyzer.py", "deal_polar.py", "spectra.py")]}
    args.output.write_text(json.dumps({"schema": 1, "metadata": metadata, "cases": rows, "errors": errors}, indent=2) + "\n")


if __name__ == "__main__":
    main()
