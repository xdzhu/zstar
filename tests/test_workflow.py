import tempfile
from pathlib import Path
import unittest

from zstar.workflow import (
    _prepare_abacus_input,
    _pyatb_band_input_command,
    _pyatb_input_command,
    discover_stages,
    generate_backend_script,
    prepare_pyatb_assets,
    reuse_reference_charge,
    run_raman_workflow,
)


class WorkflowTests(unittest.TestCase):
    def _make_tree(self, root: Path):
        for relative in (
            "0.no-move",
            "2.O/z-",
            "1.Hf/y+",
            "1.Hf/x+",
            "1.Hf/x-",
        ):
            directory = root / relative
            directory.mkdir(parents=True)
            (directory / "INPUT-scf").write_text(
                "INPUT_PARAMETERS\n", encoding="utf-8"
            )

    def test_reference_is_first_and_displacements_are_deterministic(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._make_tree(root)
            names = [stage.name for stage in discover_stages(root)]
            self.assertEqual(
                names,
                [
                    "0.no-move",
                    "1.Hf/x+",
                    "1.Hf/x-",
                    "1.Hf/y+",
                    "2.O/z-",
                ],
            )

    def test_backend_scripts_are_single_root_drivers(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._make_tree(root)
            slurm = generate_backend_script(
                root,
                backend="slurm",
                queue="compute",
                min_gap_eV=0.02,
                dry_run=True,
            )
            text = slurm.read_text(encoding="utf-8")
            self.assertIn("#SBATCH --partition=compute", text)
            self.assertIn(
                f"#SBATCH --output={root.resolve()}/.zstar/slurm-%j.out",
                text,
            )
            self.assertIn("zstar workflow run", text)
            self.assertIn("--min-gap 0.02", text)
            self.assertIn("--gap-mode path", text)
            self.assertIn("--legacy-omega-max 30", text)
            self.assertIn("--dry-run", text)
            self.assertIn("# ZStar execution backend: slurm", text)
            self.assertIn(
                "--abacus-command 'srun --ntasks=1 abacus'", text
            )
            self.assertEqual(text.count("zstar workflow run"), 1)
            self.assertNotIn(b"\r\n", slurm.read_bytes())
            self.assertTrue((root / ".zstar" / "workflow_manifest.json").is_file())
            self.assertTrue((root / ".zstar" / "backend_manifest.json").is_file())

            shell = generate_backend_script(
                root,
                backend="shell",
                tasks=2,
                cpus_per_task=4,
            )
            shell_text = shell.read_text(encoding="utf-8")
            self.assertIn("# ZStar execution backend: shell", shell_text)
            self.assertIn(
                "--abacus-command 'mpirun -np 2 abacus'", shell_text
            )

            torque = generate_backend_script(
                root,
                backend="torque",
                tasks=2,
                cpus_per_task=4,
            )
            torque_text = torque.read_text(encoding="utf-8")
            self.assertIn("#PBS -l nodes=1:ppn=8", torque_text)
            self.assertIn("# ZStar execution backend: torque", torque_text)

            divided_torque = generate_backend_script(
                root,
                backend="torque",
                output=root / "divided.pbs",
                nodes=2,
                tasks=3,
                cpus_per_task=4,
            )
            self.assertIn(
                "#PBS -l nodes=2:ppn=6",
                divided_torque.read_text(encoding="utf-8"),
            )

    def test_band_path_is_default_lightweight_gate(self):
        command = _pyatb_band_input_command(
            "pyatb_input",
            gap_mode="path",
            dimensionality=3,
            mp_density=0.08,
        )
        self.assertIn("--band", command)
        self.assertNotIn("--kmode", command)
        self.assertNotIn("--mp", command)

    def test_mp_gap_gate_is_explicit_opt_in(self):
        command = _pyatb_band_input_command(
            "pyatb_input",
            gap_mode="mp",
            dimensionality=3,
            mp_density=0.08,
        )
        self.assertIn("--band", command)
        self.assertIn("--kmode mp", command)
        self.assertIn("--mp 0.08", command)

    def test_raman_static_input_does_not_request_polarization(self):
        command = _pyatb_input_command(
            "pyatb_input",
            reference=True,
            mp_density=0.08,
            legacy_omega_max=0.1,
            polarization=False,
        )
        self.assertIn("--optical", command)
        self.assertNotIn("--polar", command)

    def test_pyatb_assets_are_copied_next_to_generated_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pyatb = root / "pyatb"
            pyatb.mkdir()
            (root / "H.upf").write_text("pseudo", encoding="utf-8")
            (root / "H.orb").write_text("orbital", encoding="utf-8")
            (root / "STRU").write_text(
                "ATOMIC_SPECIES\n"
                "H 1.0 H.upf\n\n"
                "NUMERICAL_ORBITAL\n"
                "H.orb\n\n"
                "LATTICE_CONSTANT\n"
                "1.0\n",
                encoding="utf-8",
            )

            copied = prepare_pyatb_assets(root, pyatb)

            self.assertEqual(
                {path.name for path in copied},
                {"H.upf", "H.orb"},
            )
            self.assertEqual((pyatb / "H.upf").read_text(), "pseudo")
            self.assertEqual((pyatb / "H.orb").read_text(), "orbital")

    def test_raman_workflow_dry_run_is_serial_and_resumable(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raman = root / "raman"
            stages = [
                raman / "mode-0004" / "plus",
                raman / "mode-0004" / "minus",
            ]
            for stage in stages:
                stage.mkdir(parents=True)
                (stage / "INPUT-scf").write_text(
                    "INPUT_PARAMETERS\n", encoding="utf-8"
                )
            manifest = {
                "schema": 1,
                "modes": [
                    {
                        "mode": 4,
                        "plus": str(stages[0]),
                        "minus": str(stages[1]),
                    }
                ],
            }
            (raman / "raman_manifest.json").write_text(
                __import__("json").dumps(manifest), encoding="utf-8"
            )
            states = run_raman_workflow(
                raman,
                reference_dir=root / "0.no-move",
                dry_run=True,
            )
            self.assertEqual(
                [state.name for state in states],
                ["mode-0004/plus", "mode-0004/minus"],
            )
            self.assertTrue(all(state.status == "dry-run" for state in states))
            log_text = "\n".join(
                path.read_text(encoding="utf-8")
                for path in (raman / ".zstar" / "logs").glob("*.log")
            )
            self.assertIn("pyatb_input --band --output pyatb-band", log_text)
            self.assertEqual(log_text.count("pyatb_input --band"), 1)
            self.assertNotIn("--kmode", log_text)

    def test_born_gap_gate_runs_only_for_reference(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._make_tree(root)
            from zstar.workflow import run_serial_workflow

            states = run_serial_workflow(root, dry_run=True, stop_after=2)
            log_text = "\n".join(
                path.read_text(encoding="utf-8")
                for path in (root / ".zstar" / "logs").glob("*.log")
            )
            self.assertEqual(log_text.count("pyatb_input --band"), 1)
            self.assertEqual(states[0].band, "dry-run")
            self.assertEqual(states[1].band, "reference-gated")

    def test_displaced_input_reuses_reference_charge(self):
        with tempfile.TemporaryDirectory() as tmp:
            stage = Path(tmp)
            (stage / "INPUT-scf").write_text(
                "INPUT_PARAMETERS\ninit_chg auto\nout_chg 0\n",
                encoding="utf-8",
            )
            _prepare_abacus_input(stage, reference=False)
            text = (stage / "INPUT").read_text(encoding="utf-8")
            self.assertIn("init_chg            file", text)
            self.assertIn("out_chg             1", text)
            self.assertIn("out_mat_hs2         1", text)
            self.assertIn("out_mat_r           1", text)

    def test_reference_charge_is_copied_into_abacus_output_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference = root / "0.no-move"
            target = root / "1.Hf" / "x+"
            source_output = reference / "OUT.POLAR"
            source_output.mkdir(parents=True)
            target.mkdir(parents=True)
            (source_output / "SPIN1_CHG.cube").write_text(
                "charge", encoding="utf-8"
            )
            (source_output / "POLAR-CHARGE-DENSITY.restart").write_text(
                "restart", encoding="utf-8"
            )
            (target / "INPUT").write_text(
                "INPUT_PARAMETERS\nsuffix TARGET\n", encoding="utf-8"
            )

            copied = reuse_reference_charge(reference, target)

            self.assertEqual(
                {path.name for path in copied},
                {"SPIN1_CHG.cube", "TARGET-CHARGE-DENSITY.restart"},
            )
            self.assertTrue((target / "OUT.TARGET" / "SPIN1_CHG.cube").is_file())
            self.assertTrue(
                (target / "OUT.TARGET" / "TARGET-CHARGE-DENSITY.restart").is_file()
            )
            self.assertFalse((target / "SPIN1_CHG.cube").exists())


if __name__ == "__main__":
    unittest.main()
