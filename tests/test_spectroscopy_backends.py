import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from zstar.cli import zstar_cli
from zstar.spectra import GammaModes, calculate_native_line_spectrum
from zstar.spectroscopy_backends import (
    collect_calculator_spectra,
    generate_calculator_spectra_script,
    parse_cp2k_native_spectra,
    prepare_cp2k_spectra,
    prepare_vasp_spectra,
    run_calculator_spectra,
)


CP2K_INPUT = """&GLOBAL
  PROJECT test
  RUN_TYPE ENERGY
&END GLOBAL
&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    &SCF
      &OT
      &END OT
    &END SCF
  &END DFT
  &SUBSYS
    &CELL
      ABC 10 10 10
      PERIODIC NONE
    &END CELL
    &COORD
      H 0 0 0
      H 0 0 0.74
    &END COORD
  &END SUBSYS
&END FORCE_EVAL
"""


def fake_modes() -> GammaModes:
    eigenvectors = np.zeros((6, 2, 3), dtype=complex)
    eigenvectors[3, 0, 0] = 1.0
    eigenvectors[3, 1, 0] = -1.0
    eigenvectors[4, 0, 1] = 1.0
    eigenvectors[4, 1, 1] = -1.0
    eigenvectors[5, 0, 2] = 1.0
    eigenvectors[5, 1, 2] = -1.0
    return GammaModes(
        frequencies_thz=np.asarray([0, 0, 0, 10, 11, 12], dtype=float),
        eigenvectors=eigenvectors,
        masses_amu=np.asarray([28.0, 12.0]),
        lattice_angstrom=np.eye(3) * 4.0,
        symbols=("Si", "C"),
        positions_fractional=np.asarray([[0, 0, 0], [0.25, 0.25, 0.25]]),
    )


def vasp_response(epsilon: float, charge: float) -> str:
    return (
        " MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)\n"
        " ------------------------------------------------------\n"
        f" {epsilon} 0 0\n 0 {epsilon} 0\n 0 0 {epsilon}\n"
        " BORN EFFECTIVE CHARGES (including local field effects) (in |e|)\n"
        " ------------------------------------------------------\n"
        " ion 1\n"
        f" 1 {charge} 0 0\n 2 0 {charge} 0\n 3 0 0 {charge}\n"
        " ion 2\n"
        f" 1 {-charge} 0 0\n 2 0 {-charge} 0\n 3 0 0 {-charge}\n"
        " total drift: 0 0 0\n"
    )


class SpectroscopyBackendTests(unittest.TestCase):
    @patch("zstar.spectroscopy_backends.format_calculator_spectra_status", return_value="ok")
    @patch("zstar.spectroscopy_backends.run_calculator_spectra", return_value=[])
    def test_cli_calculator_command_does_not_replace_top_level_command(
        self, mocked_run, _mocked_format
    ):
        zstar_cli(
            [
                "spectra",
                "run",
                "--root",
                "work",
                "--command",
                "vasp_std",
                "--dry-run",
            ]
        )
        self.assertEqual(mocked_run.call_args.kwargs["command"], "vasp_std")

    def test_cp2k_prepare_enables_both_native_intensities(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "input.inp"
            source.write_text(CP2K_INPUT, encoding="utf-8")
            root = prepare_cp2k_spectra(source, Path(tmp) / "spectra", dimensionality=0)
            text = (root / "calculation" / "input.inp").read_text()
            self.assertIn("RUN_TYPE VIBRATIONAL_ANALYSIS", text)
            self.assertIn("INTENSITIES TRUE", text)
            self.assertIn("PERIODIC FALSE", text)
            self.assertIn("REFERENCE COM", text)
            self.assertIn("&CENTER_COORDINATES TRUE", text)
            self.assertIn("DO_RAMAN TRUE", text)
            self.assertIn("PERIODIC_DIPOLE_OPERATOR FALSE", text)
            self.assertIn("PRECONDITIONER FULL_SINGLE_INVERSE", text)
            self.assertIn("MAX_ITER 50", text)
            states = run_calculator_spectra(root, dry_run=True)
            self.assertEqual([state.status for state in states], ["dry-run"])
            script = generate_calculator_spectra_script(
                root, backend="slurm", tasks=8, cpus_per_task=2
            )
            rendered = script.read_text()
            self.assertIn("#SBATCH --ntasks=8", rendered)
            self.assertEqual(rendered.count("zstar spectra run"), 1)
            self.assertIn("zstar spectra collect", rendered)

    def test_cp2k_native_parser_reads_three_mode_blocks(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "output.log"
            output.write_text(
                " VIB|Frequency (cm^-1)  1.000000E+03 2.000000E+03\n"
                " VIB|IR int (KM/Mole)   1.500000E+01 2.500000E+01\n"
                " VIB|Raman (A^4/amu)    3.500000E+01 4.500000E+01\n"
                " VIB|Frequency (cm^-1)  3.000000E+03\n"
                " VIB|IR int (KM/Mole)   5.500000E+01\n"
                " VIB|Raman (A^4/amu)    6.500000E+01\n",
                encoding="utf-8",
            )
            result = parse_cp2k_native_spectra(output)
            self.assertEqual(result["frequencies_cm-1"], [1000.0, 2000.0, 3000.0])
            self.assertEqual(result["ir_intensities_km_mol"], [15.0, 25.0, 55.0])
            self.assertEqual(result["raman_activities_A4_amu"], [35.0, 45.0, 65.0])
            output.write_text(
                " VIB|Frequency (cm^-1)  1.000000E+03\n"
                " VIB|IR int (KM/Mole)   1.500000E+01\n"
                " VIB|Raman (A^4/amu)    NaN\n"
            )
            with self.assertRaisesRegex(ValueError, "non-finite"):
                parse_cp2k_native_spectra(output)

    def test_cp2k_collection_writes_native_activity_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "calculation").mkdir()
            (root / "spectra_manifest.json").write_text(
                json.dumps(
                    {
                        "calculator": "cp2k",
                        "dimensionality": 0,
                        "stages": [{"name": "calculation", "path": "calculation"}],
                    }
                )
            )
            (root / "calculation" / "output.log").write_text(
                " VIB|Frequency (cm^-1)  1.000000E+03\n"
                " VIB|IR int (KM/Mole)   1.500000E+01\n"
                " VIB|Raman (A^4/amu)    3.500000E+01\n"
            )
            result = collect_calculator_spectra(root, points=101, plot=False)
            self.assertEqual(result["calculator"], "cp2k")
            self.assertTrue((root / "ir_spectrum" / "ir_modes.csv").is_file())
            self.assertTrue((root / "raman_spectrum" / "raman_modes.csv").is_file())

    def test_completed_cp2k_output_with_nan_is_a_failed_stage(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "calculation").mkdir()
            (root / "spectra_manifest.json").write_text(
                json.dumps(
                    {
                        "calculator": "cp2k",
                        "stages": [{"name": "calculation", "path": "calculation"}],
                    }
                )
            )
            (root / "calculation" / "output.log").write_text(
                " VIB|Frequency (cm^-1)  1.000000E+03\n"
                " VIB|IR int (KM/Mole)   1.500000E+01\n"
                " VIB|Raman (A^4/amu)    NaN\n"
                " PROGRAM ENDED AT test\n"
            )
            states = run_calculator_spectra(root)
            self.assertEqual(states[0].status, "failed")
            self.assertIn("non-finite", states[0].error)

    @patch("zstar.spectroscopy_backends.load_vasp_gamma_modes", return_value=fake_modes())
    def test_vasp_prepare_writes_mode_pairs_and_response_inputs(self, _mock_modes):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "input"
            source.mkdir()
            (source / "INCAR").write_text("ENCUT=520\nIBRION=6\n", encoding="utf-8")
            (source / "POSCAR").write_text("placeholder", encoding="utf-8")
            (source / "KPOINTS").write_text("Gamma\n0\nGamma\n1 1 1\n0 0 0\n")
            (source / "POTCAR").write_text("licensed-placeholder")
            modes_xml = Path(tmp) / "vasprun.xml"
            modes_xml.write_text("placeholder")
            root = prepare_vasp_spectra(
                source,
                modes_xml,
                Path(tmp) / "spectra",
                mode_numbers=[4],
            )
            manifest = json.loads((root / "spectra_manifest.json").read_text())
            self.assertEqual(len(manifest["stages"]), 3)
            for relative in ("reference", "mode-0004/plus", "mode-0004/minus"):
                incar = (root / relative / "INCAR").read_text()
                self.assertIn("LEPSILON = .TRUE.", incar)
                self.assertIn("ISYM = 0", incar)
                self.assertTrue((root / relative / "POSCAR").is_file())

    def test_two_dimensional_native_backends_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "input.inp"
            source.write_text(CP2K_INPUT)
            with self.assertRaisesRegex(ValueError, "real-space"):
                prepare_cp2k_spectra(source, Path(tmp) / "spectra", dimensionality=2)
            with self.assertRaisesRegex(ValueError, "real-space"):
                prepare_vasp_spectra(
                    Path(tmp) / "vasp",
                    Path(tmp) / "vasprun.xml",
                    Path(tmp) / "vasp-spectra",
                    dimensionality=2,
                )

    def test_native_line_spectrum_preserves_tabulated_activities(self):
        result = calculate_native_line_spectrum(
            [100.0, 200.0],
            [2.0, 4.0],
            activity_kind="IR_intensity",
            activity_unit="km/mol",
            points=101,
        )
        self.assertTrue(np.allclose(result.activities, [2.0, 4.0]))
        self.assertGreater(float(np.max(result.spectrum)), 0.0)

    @patch("zstar.spectroscopy_backends.load_vasp_gamma_modes", return_value=fake_modes())
    def test_vasp_collection_writes_both_spectra(self, _mock_modes):
        for dimensionality in (0, 3):
            with self.subTest(dimensionality=dimensionality):
                with tempfile.TemporaryDirectory() as tmp:
                    root = Path(tmp)
                    manifest = {
                        "calculator": "vasp",
                        "modes_source": str(root / "vasprun.xml"),
                        "dimensionality": dimensionality,
                        "amplitude_A_sqrt_amu": 0.02,
                        "modes": [
                            {
                                "mode": 4,
                                "plus": "mode-0004/plus",
                                "minus": "mode-0004/minus",
                            }
                        ],
                        "stages": [],
                    }
                    (root / "spectra_manifest.json").write_text(json.dumps(manifest))
                    for relative, epsilon in (
                        ("reference", 5.0),
                        ("mode-0004/plus", 5.1),
                        ("mode-0004/minus", 4.9),
                    ):
                        directory = root / relative
                        directory.mkdir(parents=True)
                        (directory / "OUTCAR").write_text(vasp_response(epsilon, 2.0))
                    result = collect_calculator_spectra(root, points=101, plot=False)
                    self.assertEqual(result["mode_numbers"], [4])
                    self.assertTrue((root / "ir_spectrum" / "ir_modes.csv").is_file())
                    self.assertTrue((root / "raman_spectrum" / "raman_modes.csv").is_file())


if __name__ == "__main__":
    unittest.main()
