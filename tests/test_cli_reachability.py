from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
import unittest

from zstar.cli import zstar_cli


class CliReachabilityTests(unittest.TestCase):
    def test_all_canonical_leaves_expose_help(self):
        commands = [
            ["bec", action, "--help"]
            for action in ("pre", "run", "stat", "post", "job")
        ]
        commands += [
            ["phonon", action, "--help"]
            for action in ("pre", "run", "stat", "post", "irrep", "job")
        ]
        commands += [
            ["spectra", action, "--help"]
            for action in ("pre", "run", "stat", "post", "job")
        ]
        commands += [
            ["dielectric", action, "--help"]
            for action in ("static", "freq", "optics")
        ]
        commands += [
            ["stru", "convert", "--help"],
            ["stru", "wyckoff", "--help"],
            ["data", "qnep", "--help"],
            ["data", "db", "--help"],
            ["skill", "install", "--help"],
            ["skill", "path", "--help"],
            ["skill", "preflight", "--help"],
            ["config", "init", "--help"],
            ["config", "show", "--help"],
            ["config", "set", "--help"],
            ["config", "check", "--help"],
            ["backend", "list", "--help"],
            ["response", "--help"],
            ["density", "--help"],
            ["pot", "--help"],
        ]
        for command in commands:
            with self.subTest(command=" ".join(command)):
                with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
                    try:
                        zstar_cli(command)
                    except SystemExit as exc:
                        self.assertEqual(exc.code, 0)

    def test_action_aliases_remain_reachable(self):
        commands = [
            ["bec", "prepare", "--help"],
            ["bec", "status", "--help"],
            ["bec", "collect", "--help"],
            ["bec", "script", "--help"],
            ["ph", "pre", "--help"],
            ["phonon", "status", "--help"],
            ["spectra", "prepare", "--help"],
            ["spectra", "status", "--help"],
            ["diel", "zero", "--help"],
            ["potential", "--help"],
        ]
        for command in commands:
            with self.subTest(command=" ".join(command)):
                with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
                    try:
                        zstar_cli(command)
                    except SystemExit as exc:
                        self.assertEqual(exc.code, 0)


if __name__ == "__main__":
    unittest.main()
