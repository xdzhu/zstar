import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from zstar.phonon_gen import run_phonopy_and_process_files


class PhononGenerationTests(unittest.TestCase):
    def test_missing_optional_job_script_does_not_block_generation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "STRU").write_text(
                "ATOMIC_SPECIES\nH 1 H.upf\n",
                encoding="utf-8",
            )
            (root / "INPUT").write_text(
                "INPUT_PARAMETERS\ncal_force 1\n",
                encoding="utf-8",
            )
            for name in ("KPT", "STRU-001"):
                (root / name).write_text(name, encoding="utf-8")
            previous = Path.cwd()
            os.chdir(root)
            try:
                with (
                    patch("zstar.phonon_gen.subprocess.run") as run,
                    patch("zstar.phonon_gen.create_symlink") as symlink,
                    patch(
                        "zstar.phonopy_stru.write_phonopy_compatible_stru"
                    ) as normalize,
                ):
                    normalize.side_effect = (
                        lambda source, destination: Path(destination).write_text(
                            Path(source).read_text(encoding="utf-8"),
                            encoding="utf-8",
                        )
                    )
                    run.return_value.returncode = 0
                    run.return_value.stdout = ""
                    run.return_value.stderr = ""
                    generated = run_phonopy_and_process_files(
                        abacus_sub="missing.sh"
                    )
            finally:
                os.chdir(previous)

            self.assertEqual(generated, [str(root / "disp-001")])
            self.assertEqual(symlink.call_count, 3)
            self.assertFalse((root / "disp-001" / "missing.sh").exists())

    def test_abacus_input_must_request_forces(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "STRU").write_text(
                "ATOMIC_SPECIES\nH 1 H.upf\n",
                encoding="utf-8",
            )
            (root / "INPUT").write_text(
                "INPUT_PARAMETERS\ncal_force 0\n",
                encoding="utf-8",
            )
            (root / "KPT").write_text("K_POINTS\n", encoding="utf-8")
            previous = Path.cwd()
            os.chdir(root)
            try:
                with self.assertRaisesRegex(ValueError, "cal_force 1"):
                    run_phonopy_and_process_files()
            finally:
                os.chdir(previous)

    def test_vasp_generation_uses_vasp_phonopy_mode(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ("POSCAR", "INCAR", "KPOINTS", "POTCAR", "POSCAR-001"):
                (root / name).write_text(name, encoding="utf-8")
            previous = Path.cwd()
            os.chdir(root)
            try:
                with (
                    patch("zstar.phonon_gen.subprocess.run") as run,
                    patch("zstar.phonon_gen.create_symlink") as symlink,
                ):
                    run.return_value.returncode = 0
                    run.return_value.stdout = ""
                    run.return_value.stderr = ""
                    generated = run_phonopy_and_process_files(
                        f_stru="POSCAR", vasp_sub="missing.sh"
                    )
            finally:
                os.chdir(previous)
            command = run.call_args.args[0]
            self.assertIn("--vasp", command)
            self.assertNotIn("--abacus", command)
            self.assertEqual(generated, [str(root / "disp-001")])
            self.assertEqual(symlink.call_count, 4)
