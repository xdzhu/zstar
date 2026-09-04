"""Recompute the BTO comparison offline using immutable archived outputs."""

import hashlib
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import shutil
import sys
import tempfile

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from tools.shared_response.cubic_bto_report import report
from zstar.shared_abacus import collect_shared_abacus, prepare_shared_abacus


def verify():
    case = Path(__file__).resolve().parent
    for name, expected in json.loads((case / "checksums.json").read_text()).items():
        if hashlib.sha256((case / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Archive checksum mismatch: {name}")
    reference = json.loads((case / "results/benchmark_summary.json").read_text())
    with tempfile.TemporaryDirectory(prefix="zstar-cubic-bto-") as temporary:
        root = Path(temporary) / "experiment"
        shutil.copytree(case / "results", root)
        shutil.copytree(case / "run", root / "seed")
        for asset in (root / "seed/assets").iterdir():
            shutil.copy2(asset, root / "seed" / asset.name)
        collected = collect_shared_abacus(root / "unified")
        np.testing.assert_allclose(collected["born_raw_e"], reference["unified_born_raw_e"], atol=1e-8, rtol=0)
        value = report(root, quiet=True)
        for key in ("max_raw_BEC_difference_e", "max_frequency_difference_cm1",
                    "legacy_total_core_hours", "unified_total_core_hours"):
            np.testing.assert_allclose(value[key], reference[key], atol=1e-5, rtol=1e-6)
        assert value["static_phonon_response_status"] == "not_validated_unstable_reference"
        with redirect_stdout(io.StringIO()):
            prepared = prepare_shared_abacus(root / "seed/STRU", root=Path(temporary) / "fresh",
                                             scf_input=root / "seed/INPUT")
        assert len(prepared["stages"]) == 3
    print("PASS: checksums, exact inputs, cubic symmetry, raw/projected BEC, independent Phonopy, frequencies, and measured costs")


if __name__ == "__main__":
    verify()
