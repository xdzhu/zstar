import tempfile
from pathlib import Path
import unittest

import numpy as np

from zstar.pyatb_compat import (
    PyATBCapabilities,
    configure_optical_input,
    configure_polarization_input,
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
    def test_one_dimensional_polarization_pads_all_berry_loops(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "Input"
            path.write_text(BASE_INPUT, encoding="utf-8")

            report = configure_polarization_input(
                path,
                dimensionality=1,
                periodic_axis=2,
            )

            text = path.read_text(encoding="utf-8")
            self.assertRegex(text, r"nk1\s+5")
            self.assertRegex(text, r"nk2\s+6")
            self.assertRegex(text, r"nk3\s+2")
            self.assertEqual(report["original_grid"], [5, 6, 1])
            self.assertEqual(report["effective_grid"], [5, 6, 2])
            self.assertEqual(report["periodic_axis"], 2)
            self.assertTrue(
                (path.parent / "zstar_pyatb_polarization_compat.json").is_file()
            )

    def test_three_dimensional_polarization_input_is_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "Input"
            path.write_text(BASE_INPUT, encoding="utf-8")

            report = configure_polarization_input(path, dimensionality=3)

            self.assertFalse(report["modified"])
            self.assertEqual(path.read_text(encoding="utf-8"), BASE_INPUT)

    def test_two_dimensional_polarization_pads_out_of_plane_loop(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "Input"
            path.write_text(BASE_INPUT, encoding="utf-8")

            report = configure_polarization_input(path, dimensionality=2)

            text = path.read_text(encoding="utf-8")
            self.assertRegex(text, r"nk1\s+5")
            self.assertRegex(text, r"nk2\s+6")
            self.assertRegex(text, r"nk3\s+2")
            self.assertEqual(report["effective_grid"], [5, 6, 2])
            self.assertEqual(report["periodic_axes"], [0, 1])

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

    def test_band_gap_recovers_manifold_gap_when_pyatb_uses_same_band(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "band_info.dat").write_text(
                "Fermi Energy (eV): 8.5656\n"
                "Band gap (eV): 0.0020\n"
                "Eigenvalue of VBM (eV): 8.5637\n"
                "Eigenvalue of CBM (eV): 8.5656\n"
                "VBM 1 (band index and k coor): 1 0 0 0\n"
                "CBM 1 (band index and k coor): 1 0.5 0 0\n",
                encoding="utf-8",
            )
            np.savetxt(
                root / "band.dat",
                np.array(
                    [
                        [-5.0, 8.40, 11.20, 15.0],
                        [-4.8, 8.57, 11.09, 15.2],
                        [-4.9, 8.50, 11.30, 15.1],
                    ]
                ),
            )

            result = read_band_gap(root, threshold_eV=0.1)

            self.assertTrue(result.insulating)
            self.assertAlmostEqual(result.vbm_eV, 8.57)
            self.assertAlmostEqual(result.cbm_eV, 11.09)
            self.assertAlmostEqual(result.gap_eV, 2.52)
            self.assertEqual(Path(result.source), (root / "band.dat").resolve())

    def test_band_gap_uses_integer_occupied_count_at_degenerate_vbm(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            band_dir = root / "Out" / "Band_Structure"
            band_dir.mkdir(parents=True)
            (root / "get_Energy.out").write_text(
                "E_FERMI (eV) = 10.689367773\n"
                "Occupied bands = 41\n"
                "NBANDS = 51\n"
                "NELEC = 82\n",
                encoding="utf-8",
            )
            (band_dir / "band_info.dat").write_text(
                "Fermi Energy (eV): 10.689367773\n"
                "Band gap (eV): 0.0000\n"
                "Eigenvalue of VBM (eV): 10.6894\n"
                "Eigenvalue of CBM (eV): 10.6894\n"
                "VBM 1 (band index and k coor): 39 0 0 0\n"
                "CBM 1 (band index and k coor): 40 0 0 0\n",
                encoding="utf-8",
            )
            bands = np.zeros((3, 43))
            bands[:, 39] = [10.68, 10.67, 10.66]
            bands[:, 40] = [10.68936783, 10.68, 10.67]
            bands[:, 41] = [14.15, 14.03508710, 14.20]
            np.savetxt(band_dir / "band.dat", bands)

            result = read_band_gap(band_dir, threshold_eV=0.01)

            self.assertTrue(result.insulating)
            self.assertAlmostEqual(result.vbm_eV, 10.68936783)
            self.assertAlmostEqual(result.cbm_eV, 14.03508710)
            self.assertAlmostEqual(result.gap_eV, 3.34571927)
            self.assertEqual(Path(result.source), (band_dir / "band.dat").resolve())

    def test_spin_band_crossing_is_metal_not_millielectronvolt_gap(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "band_info.dat").write_text(
                "Fermi Energy (eV): 12.3638\n"
                "Band gap (eV): 0.0032\n"
                "Eigenvalue of VBM (eV): 12.3610\n"
                "Eigenvalue of CBM (eV): 12.3642\n"
                "VBM 1 (band index and k coor): 1 0 0 0\n"
                "CBM 1 (band index and k coor): 1 0.5 0 0\n",
                encoding="utf-8",
            )
            bands = np.array(
                [
                    [-5.0, 12.20, 13.0],
                    [-4.9, 12.50, 13.2],
                    [-4.8, 12.36, 13.1],
                ]
            )
            np.savetxt(root / "band_up.dat", bands)
            np.savetxt(root / "band_dn.dat", bands)
            log_dir = root / "OUT.test"
            log_dir.mkdir()
            (log_dir / "running_scf.log").write_text(
                "nelec for spin up = 40.5\n"
                "nelec for spin down = 40.5\n",
                encoding="utf-8",
            )

            result = read_band_gap(root, threshold_eV=0.01)

            self.assertFalse(result.insulating)
            self.assertEqual(result.gap_eV, 0.0)
            self.assertAlmostEqual(result.vbm_eV, 12.3638)
            self.assertEqual(Path(result.source), (root / "band_up.dat").resolve())

    def test_spin_polarized_gap_uses_total_band_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "band_info.dat").write_text(
                "Fermi Energy (eV): 12.1525\n\n"
                "For nspin up:\n"
                "Band gap (eV): 1.3650\n"
                "Eigenvalue of VBM (eV): 11.8979\n"
                "Eigenvalue of CBM (eV): 13.2630\n"
                "VBM 1 (band index and k coor): 42 0 0 0\n"
                "CBM 1 (band index and k coor): 43 0 0 0\n\n"
                "For nspin down:\n"
                "Band gap (eV): 0.0074\n"
                "Eigenvalue of VBM (eV): 12.1480\n"
                "Eigenvalue of CBM (eV): 12.1553\n"
                "VBM 1 (band index and k coor): 40 0 0 0\n"
                "CBM 1 (band index and k coor): 41 0 0 0\n\n"
                "For total band:\n"
                "Band gap (eV): 0.0074\n"
                "Eigenvalue of VBM (eV): 12.1480\n"
                "Eigenvalue of CBM (eV): 12.1553\n"
                "VBM 1 (band index and k coor): 40 0 0 0\n"
                "CBM 1 (band index and k coor): 41 0 0 0\n",
                encoding="utf-8",
            )

            result = read_band_gap(root, threshold_eV=0.01)

            self.assertFalse(result.insulating)
            self.assertAlmostEqual(result.gap_eV, 0.0074)
            self.assertAlmostEqual(result.vbm_eV, 12.1480)
            self.assertAlmostEqual(result.cbm_eV, 12.1553)


if __name__ == "__main__":
    unittest.main()
