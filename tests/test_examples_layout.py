import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"


def _manifest_cases():
    manifest = json.loads((EXAMPLES / "manifest.json").read_text(encoding="utf-8"))
    return manifest["cases"]


def test_manifest_cases_are_self_describing():
    for record in _manifest_cases():
        case = EXAMPLES / record["path"]
        assert case.is_dir(), record["id"]
        assert (case / "run.sh").is_file(), record["id"]
        assert (case / "README.md").is_file(), record["id"]
        assert (case / "README.zh-CN.md").is_file(), record["id"]
        assert (case / "run").is_dir(), record["id"]
        assert (case / "results").is_dir(), record["id"]
        assert any((case / "results").rglob("*")), record["id"]


def test_abacus_cases_ship_matching_assets():
    for record in _manifest_cases():
        if not record["calculator"].startswith("ABACUS"):
            continue
        case = EXAMPLES / record["path"]
        if record.get("assets_required") is False:
            continue
        assets = case / "run" / "assets"
        asset_root = assets if assets.is_dir() else case / "run"
        assert list(asset_root.rglob("*.upf")), record["id"]
        assert list(asset_root.rglob("*.orb")), record["id"]


def test_backend_cases_declare_external_asset_boundary():
    for record in _manifest_cases():
        if record["calculator"].startswith("ABACUS"):
            continue
        case = EXAMPLES / record["path"]
        metadata = json.loads((case / "case.json").read_text(encoding="utf-8"))
        provenance = case / "ASSET_PROVENANCE.md"
        assert provenance.is_file(), record["id"]
        if record["calculator"] == "VASP":
            assert metadata["licensed_inputs_required"]
        if record["calculator"] == "CP2K":
            assert metadata["data_dir_required"] is True


def test_legacy_case_layout_is_not_reintroduced():
    legacy_names = {"input", "reference_results", "reference_spectroscopy"}
    for path in EXAMPLES.rglob("*"):
        if path.is_dir():
            assert path.name not in legacy_names, path


def test_run_directories_contain_inputs_only():
    allowed_asset_dirs = {"assets", "pp", "orb"}
    for record in _manifest_cases():
        run_dir = EXAMPLES / record["path"] / "run"
        unexpected = [
            path for path in run_dir.rglob("*")
            if path.is_dir() and path.name not in allowed_asset_dirs
        ]
        assert not unexpected, (record["id"], unexpected)
        assert not (run_dir / "work").exists(), record["id"]
        assert not (run_dir / "native").exists(), record["id"]
