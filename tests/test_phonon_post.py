import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from zstar.phonon_post import get_phonopy_params, run_eigen_irrep


class PhononPostTests(unittest.TestCase):
    def test_reads_phonopy_configuration_without_dumping_yaml(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "phonopy_disp.yaml"
            path.write_text(
                "phonopy:\n"
                "  configuration:\n"
                "    dim: [2, 2, 1]\n"
                "    symmetry_tolerance: 0.0005\n"
                "space_group:\n"
                "  type: P6_3/mmc\n",
                encoding="utf-8",
            )
            dim, tolerance, group = get_phonopy_params(path)
            self.assertEqual(dim, "2 2 1")
            self.assertAlmostEqual(tolerance, 5.0e-4)
            self.assertEqual(group, "P6_3/mmc")

    def test_postprocess_raises_when_force_logs_are_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "STRU").write_text("structure", encoding="utf-8")
            (root / "phonopy_disp.yaml").write_text(
                "phonopy:\n  configuration:\n    dim: \"1 1 1\"\n",
                encoding="utf-8",
            )
            previous = Path.cwd()
            os.chdir(root)
            try:
                with self.assertRaisesRegex(FileNotFoundError, "running"):
                    run_eigen_irrep()
            finally:
                os.chdir(previous)

    def test_postprocess_runs_force_qpoint_and_irrep_stages(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "STRU").write_text("structure", encoding="utf-8")
            (root / "phonopy_disp.yaml").write_text(
                "phonopy:\n  configuration:\n    dim: \"1 1 1\"\n",
                encoding="utf-8",
            )
            log = root / "disp-001" / "OUT.TEST" / "running_scf.log"
            log.parent.mkdir(parents=True)
            log.write_text("done", encoding="utf-8")
            previous = Path.cwd()
            os.chdir(root)
            try:
                with (
                    patch("zstar.phonon_post._run_phonopy") as run,
                    patch(
                        "zstar.phonon_post.write_phonopy_compatible_stru"
                    ) as normalize,
                ):
                    normalize.side_effect = (
                        lambda source, destination: Path(destination).write_text(
                            Path(source).read_text(encoding="utf-8"),
                            encoding="utf-8",
                        )
                    )
                    report = run_eigen_irrep()
            finally:
                os.chdir(previous)

            self.assertEqual(run.call_count, 3)
            self.assertEqual(report["force_logs"], 1)
