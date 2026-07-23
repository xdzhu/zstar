import tempfile
from pathlib import Path
import unittest

import numpy as np

from zstar.pyatb_compat import (
    PyATBCapabilities,
    configure_optical_input,
    detect_pyatb_capabilities,
    read_band_gap,
    read_static_dielectric,
)


BASE_INPUT = """INPUT_PARAMETERS
{
    nspin 1
}

POLARIZATION
{
    occ_band 4
    nk1 5
    nk2 6
    nk3 1
}
"""


class PyATBCompatTests(unittest.TestCase):
    def test_detects_adjacent_pip_target_installation(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp)
            executable = target / "bin" / "pyatb"
            source = target / "pyatb" / "io" / "default_input.py"
            executable.parent.mkdir()
            source.parent.mkdir(parents=True)
            executable.write_text("#!/usr/bin/python\n", encoding="utf-8")
            source.write_text(
                "static_dielectric_only = False\n", encoding="utf-8"
            )
            (target / "pyatb-1.2.3.dist-info").mkdir()

            capabilities = detect_pyatb_capabilities(str(executable))

            self.assertTrue(capabilities.static_dielectric_only)
            self.assertEqual(capabilities.version, "1.2.3")
            self.assertEqual(
                capabilities.detection, "adjacent pip-target installation"
            )

    def test_configures_new_direct_static_interface(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "Input"
            path.write_text(BASE_INPUT, encoding="utf-8")
            caps = PyATBCapabilities("new", "pyatb", True, "test")
            report = configure_optical_input(path, capabilities=caps)
            text = path.read_text(encoding="utf-8")
            self.assertIn("OPTICAL_CONDUCTIVITY", text)
            self.assertIn("static_dielectric_only      1", text)
            self.assertEqual(report["mode"], "direct-static")

    def test_configures_old_compact_spectrum_without_unknown_key(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "Input"
            path.write_text(BASE_INPUT, encoding="utf-8")
            caps = PyATBCapabilities("old", "pyatb", False, "test")
            report = configure_optical_input(
                path,
                capabilities=caps,
                legacy_omega_max=0.2,
                legacy_domega=0.2,
            )
            text = path.read_text(encoding="utf-8")
            self.assertIn("omega                       0.0 0.2", text)
            self.assertNotIn("static_dielectric_only", text)
            self.assertEqual(report["mode"], "legacy-compact-spectrum")

    def test_reads_both_dielectric_output_formats(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            static = root / "static_dielectric_function.dat"
            static.write_text(
                "# xx xy xz yx yy yz zx zy zz\n"
                "2 0 0 0 3 0 0 0 4\n",
                encoding="utf-8",
            )
            tensor, source = read_static_dielectric(root)
            np.testing.assert_allclose(tensor, np.diag([2, 3, 4]))
            self.assertEqual(source, static.resolve())

            static.unlink()
            legacy = root / "dielectric_function_real_part.dat"
            legacy.write_text(
                "# omega xx xy xz yx yy yz zx zy zz\n"
                "0.0 5 0 0 0 6 0 0 0 7\n",
                encoding="utf-8",
            )
            tensor, _ = read_static_dielectric(root)
            np.testing.assert_allclose(tensor, np.diag([5, 6, 7]))

    def test_band_gap_gate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            info = root / "band_info.dat"
            info.write_text(
                "Band gap (eV): 1.2500\n"
                "Eigenvalue of VBM (eV): -0.4000\n"
                "Eigenvalue of CBM (eV): 0.8500\n",
                encoding="utf-8",
            )
            result = read_band_gap(root, threshold_eV=0.1)
            self.assertTrue(result.insulating)
            self.assertAlmostEqual(result.gap_eV, 1.25)

            info.write_text("Fermi Energy (eV): 0.0\n", encoding="utf-8")
            result = read_band_gap(root, threshold_eV=0.01)
            self.assertFalse(result.insulating)
            self.assertEqual(result.gap_eV, 0.0)


if __name__ == "__main__":
    unittest.main()
