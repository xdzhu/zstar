import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from zstar.phonon_post import get_phonopy_params, run_eigen_irrep
from zstar.phonopy_stru import write_phonopy_compatible_stru


class PhononPostTests(unittest.TestCase):
    def test_irrep_view_can_omit_magnetism_without_changing_geometry(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "STRU"
            source.write_text(
                "ATOMIC_SPECIES\nFe 55.845 Fe.upf\nO 15.999 O.upf\n"
                "\nLATTICE_CONSTANT\n1.8897261255\n"
                "\nLATTICE_VECTORS\n2 0 0\n0 2 0\n0 0 2\n"
                "\nATOMIC_POSITIONS\nDirect\n"
                "\nFe\n0.0\n1\n0 0 0 m 1 1 1 mag 5\n"
                "\nO\n0.0\n1\n0.5 0.5 0.5 m 1 1 1 mag 0\n",
                encoding="utf-8",
            )
            magnetic = root / "STRU.magnetic"
            irrep = root / "STRU.irrep"
            write_phonopy_compatible_stru(source, magnetic)
            write_phonopy_compatible_stru(
                source, irrep, include_magnetism=False
            )

            magnetic_text = magnetic.read_text(encoding="utf-8")
            irrep_text = irrep.read_text(encoding="utf-8")
            self.assertIn("mag 5", magnetic_text)
            self.assertIn("mag 0", magnetic_text)
            self.assertNotIn(" mag ", irrep_text)
            for coordinate in ("0.000000000000", "0.500000000000"):
                self.assertIn(coordinate, magnetic_text)
                self.assertIn(coordinate, irrep_text)

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
                        lambda source, destination, **kwargs: Path(destination).write_text(
                            Path(source).read_text(encoding="utf-8"),
                            encoding="utf-8",
                        )
                    )
                    report = run_eigen_irrep()
            finally:
                os.chdir(previous)

            self.assertEqual(run.call_count, 3)
            self.assertEqual(normalize.call_count, 2)
            self.assertFalse(normalize.call_args_list[1].kwargs["include_magnetism"])
            self.assertEqual(report["force_logs"], 1)

    def test_bulk_nac_passes_explicit_q_direction(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "STRU").write_text("structure", encoding="utf-8")
            (root / "phonopy_disp.yaml").write_text(
                "phonopy:\n  configuration:\n    dim: '1 1 1'\n",
                encoding="utf-8",
            )
            (root / "BORN").write_text("born\n", encoding="utf-8")
            log = root / "disp-001" / "OUT.TEST" / "running_scf.log"
            log.parent.mkdir(parents=True)
            log.write_text("done", encoding="utf-8")
            previous = Path.cwd()
            os.chdir(root)
            try:
                with (
                    patch("zstar.phonon_post._run_phonopy") as run,
                    patch("zstar.phonon_post.write_phonopy_compatible_stru") as normalize,
                ):
                    normalize.side_effect = lambda source, destination, **kwargs: Path(destination).write_text("structure")
                    report = run_eigen_irrep(
                        nac=True, physical_dim=3, q_direction=(1.0, 0.0, 0.0)
                    )
            finally:
                os.chdir(previous)
            qpoint_command = run.call_args_list[1].args[0]
            self.assertIn("--nac", qpoint_command)
            self.assertIn("--q-direction", qpoint_command)
            self.assertEqual(report["q_direction"], [1.0, 0.0, 0.0])

    def test_bulk_nac_requires_born_file_before_phonopy(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "STRU").write_text("structure", encoding="utf-8")
            (root / "phonopy_disp.yaml").write_text(
                "phonopy:\n  configuration:\n    dim: '1 1 1'\n",
                encoding="utf-8",
            )
            previous = Path.cwd()
            os.chdir(root)
            try:
                with self.assertRaisesRegex(FileNotFoundError, "Copy the BEC workflow's BORN"):
                    run_eigen_irrep(nac=True, physical_dim=3)
            finally:
                os.chdir(previous)

    def test_low_dimensional_bulk_nac_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "STRU").write_text("structure", encoding="utf-8")
            (root / "phonopy_disp.yaml").write_text(
                "phonopy:\n  configuration:\n    dim: '1 1 1'\n",
                encoding="utf-8",
            )
            previous = Path.cwd()
            os.chdir(root)
            try:
                with self.assertRaisesRegex(ValueError, "2d-cutoff"):
                    run_eigen_irrep(nac=True, physical_dim=2, nac_model="gonze")
                with self.assertRaisesRegex(NotImplementedError, "true low-dimensional"):
                    run_eigen_irrep(nac=True, physical_dim=2, nac_model="2d-cutoff")
            finally:
                os.chdir(previous)

    def test_vasp_postprocess_collects_vasprun_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "POSCAR").write_text("vasp structure", encoding="utf-8")
            (root / "phonopy_disp.yaml").write_text(
                "phonopy:\n  configuration:\n    dim: '1 1 1'\n",
                encoding="utf-8",
            )
            output = root / "disp-001" / "vasprun.xml"
            output.parent.mkdir(parents=True)
            output.write_text("placeholder", encoding="utf-8")
            previous = Path.cwd()
            os.chdir(root)
            try:
                with patch("zstar.phonon_post._run_phonopy") as run:
                    report = run_eigen_irrep(f_stru="POSCAR")
            finally:
                os.chdir(previous)
            commands = [call.args[0] for call in run.call_args_list]
            self.assertEqual(Path(commands[0][2]), Path("disp-001") / "vasprun.xml")
            self.assertTrue(all("--vasp" in command for command in commands[1:]))
            self.assertEqual(report["calculator"], "vasp")
