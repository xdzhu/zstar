import importlib
import json
import os
from pathlib import Path
import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from unittest.mock import patch

from zstar.cli_frontend import handle_canonical_cli
from zstar.cli import zstar_cli
from zstar.configuration import load_config, set_config_value
from zstar.project_manifest import write_manifest


class CanonicalCliArchitectureTests(unittest.TestCase):
    def test_package_module_entry_point_is_importable(self):
        module = importlib.import_module("zstar.__main__")
        self.assertIs(module.zstar_cli, zstar_cli)

    def test_spectra_family_help_uses_canonical_actions(self):
        output = StringIO()
        with redirect_stdout(output):
            zstar_cli(["spectra", "--help"])
        self.assertIn("pre, run, job, stat, post", output.getvalue())
        self.assertNotIn("prepare,run,status", output.getvalue())

    def test_simple_families_route_to_established_handlers(self):
        calls = []
        runner = lambda argv: calls.append(list(argv))
        cases = [
            (["dielectric", "static", "--born", "BORN"], ["calc", "--born", "BORN"]),
            (["diel", "zero"], ["calc"]),
            (["stru", "convert", "--to", "vasp"], ["vasp"]),
            (["stru", "wyckoff"], ["wyckoff"]),
            (["data", "qnep", "check", "--input", "train.xyz"], ["qnep", "check", "--input", "train.xyz"]),
            (["skill", "path"], ["agent-skill", "path"]),
        ]
        for command, expected in cases:
            calls.clear()
            self.assertTrue(handle_canonical_cli(command, runner))
            self.assertEqual(calls, [expected])

    def test_bec_manifest_drives_later_actions(self):
        calls = []
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            previous = Path.cwd()
            try:
                os.chdir(root)
                self.assertTrue(
                    handle_canonical_cli(
                        ["bec", "pre", "--calculator", "abacus", "--dim", "2", "--method", "central"],
                        lambda argv: calls.append(list(argv)),
                    )
                )
                manifest = json.loads((root / ".zstar" / "bec.json").read_text())
                self.assertEqual(manifest["calculator"], "abacus")
                self.assertEqual(manifest["dimensionality"], 2)
                self.assertEqual(manifest["options"]["method"], "central")

                calls.clear()
                handle_canonical_cli(["bec", "run"], lambda argv: calls.append(list(argv)))
                self.assertEqual(calls[0][:2], ["workflow", "run"])
                self.assertIn("--dimensionality", calls[0])

                calls.clear()
                handle_canonical_cli(["bec", "post"], lambda argv: calls.append(list(argv)))
                self.assertEqual(calls[0][0], "deal")
                self.assertIn("central", calls[0])
            finally:
                os.chdir(previous)

    def test_abacus_bec_pre_and_post_execute_inside_root(self):
        calls = []
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "polar"

            def runner(argv):
                calls.append((list(argv), Path.cwd()))

            handle_canonical_cli(
                ["bec", "pre", "--root", str(root), "--calculator", "abacus"],
                runner,
            )
            self.assertEqual(calls[-1][1], root.resolve())

            calls.clear()
            handle_canonical_cli(["bec", "post", "--root", str(root)], runner)
            self.assertEqual(calls[-1][1], root.resolve())

    def test_old_ph_command_is_not_consumed_as_family_alias(self):
        self.assertFalse(handle_canonical_cli(["ph", "--stru", "STRU"], lambda _: None))

    def test_phonon_manifest_preserves_physical_dimensionality(self):
        calls = []
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "STRU").write_text("ATOMIC_SPECIES\n", encoding="utf-8")
            handle_canonical_cli(
                ["phonon", "pre", "--root", str(root), "--physical-dim", "2"],
                lambda argv: calls.append(list(argv)),
            )
            self.assertEqual(calls[0][0], "ph")
            manifest = json.loads((root / ".zstar" / "phonon.json").read_text())
            self.assertEqual(manifest["dimensionality"], 2)

            calls.clear()
            handle_canonical_cli(
                ["phonon", "post", "--root", str(root)],
                lambda argv: calls.append(list(argv)),
            )
            self.assertEqual(calls[0][0], "postph")
            self.assertIn("--physical-dim", calls[0])
            self.assertEqual(calls[0][calls[0].index("--physical-dim") + 1], "2")
            self.assertIn("--stru", calls[0])

    def test_abacus_spectra_manifest_routes_status_and_job(self):
        calls = []
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "raman"
            handle_canonical_cli(
                [
                    "spectra", "pre", "--calculator", "abacus",
                    "--kind", "raman", "--root", str(root),
                    "--stru", "STRU", "--qpoints", "qpoints.yaml",
                ],
                lambda argv: calls.append(list(argv)),
            )
            self.assertEqual(calls[0][:2], ["raman", "prepare"])
            self.assertTrue((root / ".zstar" / "spectra.json").is_file())

            calls.clear()
            handle_canonical_cli(
                ["spectra", "stat", "--root", str(root)],
                lambda argv: calls.append(list(argv)),
            )
            self.assertEqual(calls, [["raman", "status", "--raman-dir", str(root)]])

            handle_canonical_cli(
                ["spectra", "job", "--root", str(root), "--system", "slurm"],
                lambda argv: calls.append(list(argv)),
            )
            script = root / "run_zstar_spectra.slurm"
            self.assertTrue(script.is_file())
            self.assertIn("zstar spectra run", script.read_text(encoding="utf-8"))

    def test_canonical_bec_job_accepts_pbs_alias_and_writes_pbs_driver(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "0.no-move").mkdir()
            (root / "1.Ba" / "x+").mkdir(parents=True)
            for path in (root / "0.no-move", root / "1.Ba" / "x+"):
                (path / "INPUT-scf").write_text(
                    "INPUT_PARAMETERS\n", encoding="utf-8"
                )
            handle_canonical_cli(
                [
                    "bec", "job", "--root", str(root), "--system", "pbs",
                    "--tasks", "2", "--cpus-per-task", "4",
                ],
                zstar_cli,
            )
            driver = root / "run_zstar_born.pbs"
            self.assertTrue(driver.is_file())
            text = driver.read_text(encoding="utf-8")
            self.assertIn("# ZStar execution backend: torque", text)
            self.assertIn("#PBS -l nodes=1:ppn=8", text)
            self.assertIn("mpirun -np 2 abacus", text)

    def test_canonical_phonon_job_accepts_pbs_alias(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            stage = root / "disp-001"
            stage.mkdir()
            (stage / "STRU").write_text("ATOMIC_SPECIES\n", encoding="utf-8")
            write_manifest("phonon", root=root, calculator="abacus", dimensionality=3)
            handle_canonical_cli(
                [
                    "phonon", "job", "--root", str(root), "--system", "pbs",
                    "--tasks", "2", "--cpus-per-task", "4",
                ],
                zstar_cli,
            )
            driver = root / "run_zstar_phonon.pbs"
            self.assertTrue(driver.is_file())
            text = driver.read_text(encoding="utf-8")
            self.assertIn("# ZStar phonon execution system: torque", text)
            self.assertIn("#PBS -l nodes=1:ppn=8", text)

    def test_canonical_spectra_job_accepts_pbs_alias(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_manifest(
                "spectra", root=root, calculator="abacus", dimensionality=3,
                options={"kind": "all"},
            )
            handle_canonical_cli(
                [
                    "spectra", "job", "--root", str(root), "--system", "pbs",
                    "--tasks", "2", "--cpus-per-task", "4",
                ],
                zstar_cli,
            )
            driver = root / "run_zstar_spectra.pbs"
            self.assertTrue(driver.is_file())
            text = driver.read_text(encoding="utf-8")
            self.assertIn("# ZStar spectroscopy execution system: torque", text)
            self.assertIn("#PBS -l nodes=1:ppn=8", text)

    def test_project_and_environment_configuration_precedence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            set_config_value("executables.abacus", "/project/abacus", root=root)
            self.assertEqual(load_config(root)["executables"]["abacus"], "/project/abacus")
            with patch.dict(os.environ, {"ZSTAR_ABACUS_EXECUTABLE": "/env/abacus"}):
                self.assertEqual(load_config(root)["executables"]["abacus"], "/env/abacus")

    def test_top_level_help_exposes_canonical_families(self):
        output = StringIO()
        with self.assertRaises(SystemExit), redirect_stdout(output):
            zstar_cli(["--help"])
        text = output.getvalue()
        for family in (
            "bec", "phonon", "spectra", "dielectric", "pot", "backend",
            "config", "response", "density", "stru", "data", "skill",
        ):
            self.assertIn(family, text)
        command_help = text.split("options:", maxsplit=1)[0]
        self.assertNotIn("gen                 Generate polarization data", command_help)
        self.assertNotIn("postph", command_help)
        self.assertNotIn("cp2k-bec", command_help)
        self.assertNotIn("md                  Calculate", command_help)

    def test_backend_list_check_reports_configured_executables(self):
        output = StringIO()
        with redirect_stdout(output):
            zstar_cli(["backend", "list", "--check"])
        text = output.getvalue()
        self.assertIn("Configured executables:", text)
        self.assertIn("abacus", text)


if __name__ == "__main__":
    unittest.main()
