import json
from pathlib import Path
import tempfile
import unittest

from zstar.agent_skill import (
    SKILL_NAME,
    install_agent_skill,
    packaged_skill_path,
    preflight_report,
)


class AgentSkillTests(unittest.TestCase):
    def test_skill_name_and_frontmatter_match_directory(self):
        skill = packaged_skill_path()
        text = (skill / "SKILL.md").read_text(encoding="utf-8")
        self.assertEqual(skill.name, SKILL_NAME)
        self.assertIn(f"name: {SKILL_NAME}\n", text)
        self.assertNotIn("TODO", text)
        self.assertTrue((skill / "references" / "spectroscopy.md").is_file())
        self.assertFalse(
            (skill / "references" / "spectroscopy-and-md.md").exists()
        )

    def test_preflight_reports_required_input_and_existing_state(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            missing = preflight_report(root, lane="bec", dimensionality="bulk")
            self.assertFalse(missing["ready"])
            self.assertIn("STRU is required", missing["blockers"][0])

            (root / "STRU").write_text("ATOMIC_SPECIES\n", encoding="utf-8")
            stage_dir = root / ".zstar" / "stages"
            stage_dir.mkdir(parents=True)
            (stage_dir / "0.no-move.json").write_text(
                json.dumps({"status": "completed"}), encoding="utf-8"
            )
            ready = preflight_report(root, lane="bec", dimensionality="bulk")
            self.assertTrue(ready["ready"])
            self.assertEqual(ready["state"]["counts"], {"completed": 1})

    def test_install_copies_valid_skill_and_requires_force_to_replace(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            installed = install_agent_skill(root)
            self.assertEqual(installed.name, SKILL_NAME)
            self.assertTrue((installed / "SKILL.md").is_file())
            self.assertFalse(any(installed.rglob("*.pyc")))
            self.assertFalse(any(path.name == "__pycache__" for path in installed.rglob("*")))
            with self.assertRaises(FileExistsError):
                install_agent_skill(root)
            replaced = install_agent_skill(root, force=True)
            self.assertEqual(replaced, installed)

    def test_preflight_accepts_supported_1d_workflow_with_scope_warning(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "STRU").write_text("placeholder\n", encoding="utf-8")
            report = preflight_report(root, lane="bec", dimensionality="1d")

            self.assertTrue(report["ready"])
            self.assertEqual(report["dimensionality"], "1d")
            self.assertFalse(report["blockers"])
            self.assertTrue(any("bulk NAC" in item for item in report["warnings"]))

    def test_preflight_allows_1d_database_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = preflight_report(tmp, lane="database", dimensionality="1d")

            self.assertTrue(report["ready"])
            self.assertTrue(any("bulk NAC" in item for item in report["warnings"]))

    def test_preflight_describes_molecular_bec_as_apt(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "STRU").write_text("placeholder\n", encoding="utf-8")
            report = preflight_report(root, lane="bec", dimensionality="molecule")

            self.assertTrue(report["ready"])
            self.assertTrue(
                any("atomic polar tensors" in item for item in report["warnings"])
            )


if __name__ == "__main__":
    unittest.main()
