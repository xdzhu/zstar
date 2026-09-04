"""Prepare isolated, provenance-tracked DFT controls without editing archives."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np

from zstar.shared_response import read_structure, write_structure

BASE = Path("/home/zhuxd/abacus/agent-runs")
ROOT = BASE / "20260904-shared-response-benchmark"
SEED = BASE / "20260827-in2se3-pbesol-bec/relax-explicit-d3/STRU"


def prepare_relax():
    atoms = read_structure(SEED)
    cart = atoms.positions
    cart[:, 2] -= (cart[:, 2].max() + cart[:, 2].min()) / 2
    a, height = 4.106, 30.0
    fractional = atoms.scaled_positions.copy()
    atoms.cell = [[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, height]]
    fractional[:, 2] = (cart[:, 2] + height / 2) / height
    atoms.scaled_positions = fractional
    # ABBCA, read bottom to top, is FE-ZB-prime up to a cyclic site relabeling.
    layers = np.argsort(atoms.positions[:, 2])
    for mode, folder in (("cell-relax", "in2se3_nc2017/relax"),
                         ("relax", "in2se3_nc2017/fixed_a_control")):
        stage = ROOT / folder
        stage.mkdir(parents=True, exist_ok=False)
        write_structure(SEED, stage / "STRU", atoms)
        (stage / "KPT").write_text("K_POINTS\n0\nGamma\n12 12 1 0 0 0\n")
        params = {"suffix": "SHARED_IN2SE3", "calculation": mode,
                  "basis_type": "lcao", "dft_functional": "pbe", "ecutwfc": 100,
                  "ks_solver": "genelpa", "scf_thr": "1e-8", "scf_nmax": 200,
                  "smearing_method": "gauss", "smearing_sigma": 0.005,
                  "mixing_type": "pulay", "mixing_beta": 0.3,
                  "symmetry": 1, "cal_force": 1, "cal_stress": int(mode == "cell-relax"),
                  "force_thr_ev": 0.005, "stress_thr": 0.5, "relax_nmax": 150,
                  "fixed_axes": "c", "out_stru": 1, "out_chg": "1 10",
                  "efield_flag": 1, "dip_cor_flag": 1, "efield_dir": 2,
                  "efield_pos_max": 0.95, "efield_pos_dec": 0.1, "efield_amp": 0}
        (stage / "INPUT").write_text("INPUT_PARAMETERS\n" + "".join(
            f"{key:<24}{value}\n" for key, value in params.items()))
        provenance = {"reference_doi": "10.1038/ncomms14956", "seed": str(SEED),
                      "seed_sha256": hashlib.sha256(SEED.read_bytes()).hexdigest(),
                      "in_plane_lattice_A": a, "cell_height_A": height,
                      "vacuum_gap_A": float(height - np.ptp(cart[:, 2])),
                      "stacking_bottom_to_top": fractional[layers, :2].tolist(),
                      "relaxation_complete": False,
                      "matched": ["PBE", "Gamma 12x12x1", "vacuum >15 A", "dipole correction", "force <0.005 eV/A"],
                      "different_from_paper": ["ABACUS LCAO rather than VASP PAW", "ONCV potentials and localized orbitals", "100 Ry density cutoff", "smearing 0.005 Ry (not specified in paper)"]}
        (stage / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
        print(stage)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["relax", "symmetry-relax"])
    args = parser.parse_args()
    if args.action == "relax":
        prepare_relax()
    else:
        from zstar.workflow import _set_abacus_parameter
        old = ROOT / "in2se3_nc2017/relax"
        log = old / "OUT.SHARED_IN2SE3/running_cell-relax.log"
        if "Relaxation is converged!" not in log.read_text():
            raise RuntimeError("Cell relaxation is not converged")
        source = old / "OUT.SHARED_IN2SE3/STRU_ION_D"
        original = read_structure(source)
        atoms = original.copy()
        area = np.linalg.norm(np.cross(atoms.cell[0], atoms.cell[1]))
        a = np.sqrt(area * 2 / np.sqrt(3))
        atoms.cell = [[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, 30]]
        frac = atoms.scaled_positions.copy()
        frac[:, :2] = [[0, 0], [2/3, 1/3], [1/3, 2/3], [2/3, 1/3], [1/3, 2/3]]
        atoms.scaled_positions = frac
        movement = atoms.positions - original.positions
        if np.max(np.abs(movement)) > 1e-4 or np.max(np.abs(atoms.cell - original.cell)) > 1e-4:
            raise ValueError("Large symmetry idealization requires explicit investigation")
        stage = ROOT / "in2se3_nc2017/relax_symmetry_verified"
        stage.mkdir(exist_ok=False)
        write_structure(source, stage / "STRU", atoms)
        shutil.copy2(old / "KPT", stage / "KPT")
        text = _set_abacus_parameter((old / "INPUT").read_text(), "calculation", "relax")
        text = _set_abacus_parameter(text, "cal_stress", "0")
        (stage / "INPUT").write_text(text)
        (stage / "provenance.json").write_text(json.dumps({
            "source": str(source), "reason": "Remove numerical cell shear after converged cell relaxation, then recheck forces with a fresh fixed-cell relaxation",
            "maximum_cartesian_atom_change_A": float(np.max(np.abs(movement))),
            "maximum_cell_component_change_A": float(np.max(np.abs(atoms.cell - original.cell))),
            "lattice_a_A": a, "fresh_relaxation_required": True}, indent=2) + "\n")
        print(stage)


if __name__ == "__main__":
    main()
