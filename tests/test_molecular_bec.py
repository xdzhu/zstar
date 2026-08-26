import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from zstar.molecular_bec import _ATOM_DIR_RE, calculate_molecular_apt
from zstar.response_schema import response_record_from_bec_result
from zstar.spectra import ELEMENTARY_CHARGE


def write_polarization(stage: Path, values, quanta=(1.0, 1.0, 1.0)) -> None:
    output = stage / "pyatb" / "Out" / "Polarization"
    output.mkdir(parents=True)
    lines = [
        f"direction is in {axis}, P = {value:.16e} (mod {quantum:.16e}) C/m^2\n"
        for axis, value, quantum in zip("abc", values, quanta)
    ]
    (output / "polarization.dat").write_text("".join(lines), encoding="utf-8")


class MolecularAptTests(unittest.TestCase):
    def test_reference_directory_is_not_an_atom_directory(self):
        self.assertIsNone(_ATOM_DIR_RE.match("0.no-move"))
        self.assertEqual(_ATOM_DIR_RE.match("2.H").groups(), ("2", "H"))

    def test_central_difference_uses_dipole_derivative_in_cartesian_basis(self):
        expected = np.asarray(
            [[-0.52, 0.07, 0.01], [0.03, -0.41, 0.09], [-0.02, 0.05, -0.33]]
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference = root / "0.no-move"
            write_polarization(reference, (0.0, 0.0, 0.0))
            input_json = reference / "pyatb" / "Out" / "input.json"
            input_json.parent.mkdir(parents=True, exist_ok=True)
            input_json.write_text(
                json.dumps({"lattice_vector": (np.eye(3) * 10.0).tolist()}),
                encoding="utf-8",
            )

            atom = root / "1.O"
            displacement_m = 0.01e-10
            volume_m3 = 1000.0e-30
            for row, direction in zip(expected, "xyz"):
                delta = row * (2.0 * displacement_m * ELEMENTARY_CHARGE) / volume_m3
                write_polarization(atom / f"{direction}+", 0.5 * delta)
                write_polarization(atom / f"{direction}-", -0.5 * delta)

            actual, diagnostics = calculate_molecular_apt(
                atom, reference, method="central", displacement_angstrom=0.01
            )
            self.assertTrue(np.allclose(actual, expected, atol=1.0e-10))
            self.assertEqual(diagnostics["method"], "central")
            self.assertEqual(
                diagnostics["directions"]["x"]["branch_shifts"], [0, 0, 0]
            )

    def test_dimensionality_zero_uses_atomic_polar_tensor_name(self):
        data = {
            "backend": "abacus-pyatb",
            "atoms": [{"index": 1, "label": "H", "tensor": np.eye(3).tolist()}],
            "tensor_convention": "rows=displacement; columns=dipole",
        }
        molecular = response_record_from_bec_result(data, dimensionality=0)
        bulk = response_record_from_bec_result(data, dimensionality=3)
        self.assertEqual(molecular.quantities[0].name, "atomic_polar_tensor")
        self.assertEqual(bulk.quantities[0].name, "born_effective_charge")


if __name__ == "__main__":
    unittest.main()
