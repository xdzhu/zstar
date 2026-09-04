"""Geometry and accounting guards for the cubic-BTO experiment."""

import json
from pathlib import Path

import numpy as np
import pytest

from tools.shared_response.cubic_bto_benchmark import prepare
from tools.shared_response.cubic_bto_report import direct_cartesian_born, ledger
from zstar.shared_abacus import load_manifest
from zstar.shared_response import DEFAULT_DISTANCE, make_phonopy, read_structure, symmetry_operations
from zstar.workflow import discover_stages


@pytest.fixture(scope="module")
def experiment(tmp_path_factory):
    root = tmp_path_factory.mktemp("cubic-bto")
    source = Path(__file__).resolve().parents[2] / "examples/3d_bulk/BaTiO3_cubic/run"
    prepare(root, source)
    return root, source


def test_real_stage_counts_and_rank(experiment):
    root, _ = experiment
    meta = load_manifest(root / "unified")
    assert len(meta["stages"]) == 3
    assert len(discover_stages(root / "legacy_bec")) == 10
    assert len(list((root / "legacy_phonon").glob("disp-*"))) == 3
    atoms = read_structure(root / "unified/0.no-move/STRU")
    ops = symmetry_operations(make_phonopy(atoms))
    for s in meta["stages"]:
        orbit = [r @ s["displacement_A"] for r, perm in ops if perm[s["atom"]] == s["atom"]]
        assert np.linalg.matrix_rank(orbit) == 3
        assert np.isclose(s["distance_A"], DEFAULT_DISTANCE, atol=1e-12, rtol=0)


def test_matched_inputs_and_private_assets(experiment):
    root, source = experiment
    assert (root / "seed/assets/O.upf").read_bytes() == (source / "assets/O.upf").read_bytes()
    assert not any(p.is_symlink() for p in root.rglob("*"))
    for stage in discover_stages(root / "legacy_bec"):
        assert (stage.path / "KPT").read_bytes() == (root / "seed/KPT").read_bytes()
        text = (stage.path / "INPUT-scf").read_text()
        values = {v[0]: v[1] for line in text.splitlines() if len(v := line.split()) >= 2 and not line.startswith("#")}
        assert values["cal_force"] == "0"
        assert values["kspacing"] == "0"
        assert values["dft_functional"] == "pbesol"
        assert float(values["scf_thr"]) == 1e-7
    assert (root / "legacy_phonon/disp-001/INPUT").read_text().split("out_mat_hs2", 1)[1].split()[0] == "0"


def test_existing_experiment_is_not_overwritten(experiment):
    root, source = experiment
    with pytest.raises(FileExistsError):
        prepare(root, source)


def test_independent_cartesian_reference(monkeypatch, experiment):
    root, _ = experiment
    route = root / "legacy_bec"
    entries = json.loads((route / "displacements.json").read_text())
    known = {0: np.eye(3) * 3, 1: np.eye(3) * 7, 2: np.diag([-6., -2., -2.])}
    monkeypatch.setattr("tools.shared_response.cubic_bto_report._dipole_changes",
                        lambda *_: [known[s["atom"]] @ s["displacement_A"] for s in entries])
    tensor, direct = direct_cartesian_born(route, entries, read_structure(route / "STRU"))
    np.testing.assert_allclose(tensor[0], known[0], atol=1e-12)
    np.testing.assert_allclose(tensor[1], known[1], atol=1e-12)
    np.testing.assert_allclose(np.diagonal(tensor[2:], axis1=1, axis2=2),
                               [[-6, -2, -2], [-2, -6, -2], [-2, -2, -6]], atol=1e-12)
    np.testing.assert_allclose(tensor.sum(axis=0), 0, atol=1e-12)


def test_missing_completion_cannot_be_reported(tmp_path):
    with pytest.raises(ValueError, match="Not completed"):
        ledger(tmp_path, 4)
