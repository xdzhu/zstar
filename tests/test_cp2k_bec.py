import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from zstar.cp2k_bec import (
    DEBYE_PER_E_ANGSTROM,
    collect_cp2k_bec,
    compare_cp2k_bec,
    ensure_moments,
    ensure_periodic_moments,
    parse_cp2k_moment,
    parse_native_apt,
    prepare_cp2k_bec,
    prepare_native_apt,
    unwrap_dipole_delta,
    validate_cp2k_bec_input,
)


MINIMAL_INPUT = """&GLOBAL
  PROJECT test
  RUN_TYPE ENERGY_FORCE
&END GLOBAL
&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &SCF
      EPS_SCF 1.0E-8
      SCF_GUESS ATOMIC
      &OT
      &END OT
    &END SCF
  &END DFT
  &SUBSYS
    &CELL
      ABC [angstrom] 5.0 5.0 5.0
    &END CELL
    &COORD [angstrom]
      Mg 0.0 0.0 0.0
      O  2.5 2.5 2.5
    &END COORD
  &END SUBSYS
&END FORCE_EVAL
"""


def cp2k_output(dipole, quantum=24.01602136):
    return (
        "Dipole vectors are based on the periodic (Berry phase) operator.\n"
        f"  [X] [ {quantum: .8f} 0.00000000 0.00000000 ] [i]\n"
        f"  [Y]=[ 0.00000000 {quantum: .8f} 0.00000000 ]*[j]\n"
        f"  [Z] [ 0.00000000 0.00000000 {quantum: .8f} ] [k]\n"
        "  Dipole moment [Debye]\n"
        f"    X= {dipole[0]: .10E} Y= {dipole[1]: .10E} "
        f"Z= {dipole[2]: .10E} Total= 0.0\n"
        "  PROGRAM ENDED AT 2026-01-01 00:00:00\n"
    )


def cp2k_molecular_output(dipole):
    return (
        "  Dipole moment [Debye]\n"
        f"    X= {dipole[0]: .10E} Y= {dipole[1]: .10E} "
        f"Z= {dipole[2]: .10E} Total= 0.0\n"
        "  PROGRAM ENDED AT 2026-01-01 00:00:00\n"
    )


class Cp2kBecTests(unittest.TestCase):
    def test_moments_are_injected_under_dft_print(self):
        rendered = ensure_periodic_moments(MINIMAL_INPUT)
        self.assertIn("&MOMENTS", rendered)
        self.assertIn("PERIODIC TRUE", rendered)
        self.assertLess(rendered.index("&MOMENTS"), rendered.index("&END DFT"))

    def test_nonperiodic_moments_use_center_of_mass_reference(self):
        rendered = ensure_moments(MINIMAL_INPUT, periodic=False)
        self.assertIn("PERIODIC FALSE", rendered)
        self.assertIn("REFERENCE COM", rendered)
        self.assertLess(rendered.index("&MOMENTS"), rendered.index("&END DFT"))

    def test_prepare_generates_central_stages_and_restart_guess(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "input.inp"
            source.write_text(MINIMAL_INPUT, encoding="utf-8")
            root = prepare_cp2k_bec(
                source,
                Path(tmp) / "work",
                atoms="1",
                method="central",
                displacement_angstrom=0.01,
            )
            manifest = json.loads((root / "cp2k_bec_manifest.json").read_text())
            self.assertEqual(len(manifest["stages"]), 7)
            plus = (root / "atom-0001-Mg" / "x-plus" / "input.inp").read_text()
            minus = (root / "atom-0001-Mg" / "x-minus" / "input.inp").read_text()
            self.assertIn("Mg      0.010000000000", plus)
            self.assertIn("Mg     -0.010000000000", minus)
            self.assertIn("SCF_GUESS RESTART", plus)
            self.assertIn("WFN_RESTART_FILE_NAME reference-RESTART.wfn", plus)

    def test_parse_and_unwrap_periodic_dipole(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "output.log"
            path.write_text(cp2k_output((23.9, 0.2, -0.3)), encoding="utf-8")
            moment = parse_cp2k_moment(path)
            delta, shifts = unwrap_dipole_delta(
                moment.dipole_debye, (0.1, 0.0, 0.0), moment.quantum_debye
            )
            self.assertTrue(np.allclose(delta, (-0.21602136, 0.2, -0.3)))
            self.assertEqual(shifts.tolist(), [1, 0, 0])

    def test_collect_preserves_displacement_row_convention(self):
        expected = np.array(
            [[2.1, 0.2, -0.3], [0.4, 1.8, 0.5], [-0.6, 0.7, 2.4]]
        )
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "input.inp"
            source.write_text(MINIMAL_INPUT, encoding="utf-8")
            root = prepare_cp2k_bec(source, Path(tmp) / "work", atoms="1")
            (root / "reference" / "output.log").write_text(
                cp2k_output((0.0, 0.0, 0.0)), encoding="utf-8"
            )
            scale = 0.01 * DEBYE_PER_E_ANGSTROM
            for row, direction in zip(expected, ("x", "y", "z")):
                atom_dir = root / "atom-0001-Mg"
                (atom_dir / f"{direction}-plus" / "output.log").write_text(
                    cp2k_output(scale * row), encoding="utf-8"
                )
                (atom_dir / f"{direction}-minus" / "output.log").write_text(
                    cp2k_output(-scale * row), encoding="utf-8"
                )
            result = collect_cp2k_bec(root)
            self.assertEqual(result["sum_scope"], "selected_atoms")
            actual = np.asarray(result["atoms"][0]["tensor"])
            self.assertTrue(np.allclose(actual, expected, atol=1e-9))
            self.assertTrue((root / "Z-BORN-all.out").is_file())
            self.assertIn(
                "2.10000000",
                (root / "Z-BORN-all.out").read_text(encoding="utf-8"),
            )
            self.assertTrue((root / "zstar_response.json").is_file())
            response = json.loads((root / "zstar_response.json").read_text())
            self.assertEqual(response["schema"], "zstar-response")
            self.assertEqual(response["backend"], "cp2k")
            self.assertEqual(response["quantities"][0]["shape"], [1, 3, 3])

    def test_collect_molecular_apt_without_periodic_quantum(self):
        expected = np.array(
            [[-0.7, 0.1, 0.0], [0.2, -0.5, 0.3], [0.0, -0.1, -0.4]]
        )
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "input.inp"
            source.write_text(MINIMAL_INPUT, encoding="utf-8")
            root = prepare_cp2k_bec(
                source,
                Path(tmp) / "molecule",
                atoms="1",
                dimensionality=0,
            )
            manifest = json.loads((root / "cp2k_bec_manifest.json").read_text())
            self.assertEqual(manifest["dimensionality"], 0)
            prepared = (root / "reference" / "input.inp").read_text()
            self.assertIn("PERIODIC FALSE", prepared)
            self.assertIn("REFERENCE COM", prepared)

            (root / "reference" / "output.log").write_text(
                cp2k_molecular_output((0.3, -0.2, 0.1)), encoding="utf-8"
            )
            scale = 0.01 * DEBYE_PER_E_ANGSTROM
            for row, direction in zip(expected, ("x", "y", "z")):
                atom_dir = root / "atom-0001-Mg"
                (atom_dir / f"{direction}-plus" / "output.log").write_text(
                    cp2k_molecular_output(np.array((0.3, -0.2, 0.1)) + scale * row),
                    encoding="utf-8",
                )
                (atom_dir / f"{direction}-minus" / "output.log").write_text(
                    cp2k_molecular_output(np.array((0.3, -0.2, 0.1)) - scale * row),
                    encoding="utf-8",
                )

            result = collect_cp2k_bec(root)
            self.assertEqual(result["dimensionality"], 0)
            self.assertEqual(result["quantity"], "atomic_polar_tensor")
            self.assertAlmostEqual(result["atoms"][0]["gapt"], -1.6 / 3.0)
            self.assertTrue(
                np.allclose(np.asarray(result["atoms"][0]["tensor"]), expected)
            )
            response = json.loads((root / "zstar_response.json").read_text())
            self.assertEqual(response["dimensionality"]["value"], 0)
            self.assertEqual(response["quantities"][0]["name"], "atomic_polar_tensor")

    def test_native_apt_parser_and_comparison(self):
        native_text = """     1    Mg        2.0000000000
        2.0000000000        0.1000000000        0.0000000000
        0.2000000000        2.1000000000        0.0000000000
        0.0000000000        0.0000000000        1.9000000000

Sum of Born charges: 2.0
"""
        with tempfile.TemporaryDirectory() as tmp:
            native = Path(tmp) / "native.data"
            native.write_text(native_text, encoding="utf-8")
            parsed = parse_native_apt(native)
            self.assertEqual(parsed[0]["label"], "Mg")
            self.assertEqual(parsed[0]["tensor"][0][1], 0.2)
            self.assertEqual(parsed[0]["tensor"][1][0], 0.1)
            zstar = Path(tmp) / "cp2k_bec.json"
            zstar.write_text(
                json.dumps({"atoms": [{"index": 1, "label": "Mg", "tensor": parsed[0]["tensor"]}]}),
                encoding="utf-8",
            )
            comparison = compare_cp2k_bec(zstar, native)
            self.assertEqual(comparison["max_abs"], 0.0)
            self.assertEqual(comparison["components"], 9)
            self.assertAlmostEqual(
                comparison["per_atom"][0]["zstar_gapt"], 2.0
            )
            self.assertAlmostEqual(
                comparison["per_atom"][0]["gapt_difference"], 0.0
            )
            self.assertEqual(comparison["zstar_acoustic_sum_max_abs"], 2.1)

    def test_native_input_and_domain_guards(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "input.inp"
            source.write_text(MINIMAL_INPUT, encoding="utf-8")
            root = prepare_native_apt(source, Path(tmp) / "native")
            rendered = (root / "input.inp").read_text()
            self.assertIn("APT_FD TRUE", rendered)
            self.assertIn("APT_FD_DE 0.0003", rendered)
            self.assertIn("RUN_TYPE ENERGY", rendered)
            existing = MINIMAL_INPUT.replace(
                "&END FORCE_EVAL",
                "  &PROPERTIES\n    &LINRES\n    &END LINRES\n  &END PROPERTIES\n&END FORCE_EVAL",
            )
            source.write_text(existing, encoding="utf-8")
            merged = prepare_native_apt(source, Path(tmp) / "merged")
            merged_text = (merged / "input.inp").read_text()
            self.assertEqual(merged_text.count("&PROPERTIES"), 1)
            self.assertIn("&DCDR", merged_text)
        with self.assertRaisesRegex(ValueError, "Gamma-point"):
            validate_cp2k_bec_input(
                MINIMAL_INPUT.replace("&SCF", "&KPOINTS\n&END KPOINTS\n    &SCF", 1)
            )


if __name__ == "__main__":
    unittest.main()
