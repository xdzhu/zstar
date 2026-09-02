from pathlib import Path
import tempfile
import unittest

from zstar.phonon_workflow import (
    discover_phonon_stages,
    generate_phonon_script,
    phonon_workflow_status,
    run_phonon_workflow,
)


class PhononWorkflowTests(unittest.TestCase):
    def _folders(self, root: Path):
        for number in (2, 1, 10):
            stage = root / f"disp-{number:03d}"
            stage.mkdir()
            (stage / "STRU").write_text("placeholder\n", encoding="utf-8")

    def test_discovery_is_numeric_and_dry_run_is_serial(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._folders(root)
            stages = discover_phonon_stages(root)
            self.assertEqual([item[0].name for item in stages], ["disp-001", "disp-002", "disp-010"])
            states = run_phonon_workflow(root, dry_run=True, stop_after=2)
            self.assertEqual([item.status for item in states], ["dry-run", "dry-run", "pending"])
            self.assertTrue((root / ".zstar" / "phonon_state.json").is_file())

    def test_status_does_not_treat_stale_dry_run_as_complete(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._folders(root)
            run_phonon_workflow(root, dry_run=True, stop_after=1)
            states = phonon_workflow_status(root)
            self.assertEqual(states[0].status, "dry-run")
            self.assertEqual(states[1].status, "pending")

    def test_slurm_script_uses_one_canonical_driver(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._folders(root)
            script = generate_phonon_script(root, backend="slurm", tasks=20)
            text = script.read_text(encoding="utf-8")
            self.assertIn("#SBATCH --ntasks=20", text)
            self.assertIn("srun --ntasks=20 abacus", text)
            self.assertEqual(text.count("zstar phonon run"), 1)


if __name__ == "__main__":
    unittest.main()
