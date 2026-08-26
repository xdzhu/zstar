from pathlib import Path
import tempfile
import unittest

import numpy as np
import yaml

from zstar.interoperability import (
    intrinsic_polarizability_from_supercell,
    response_record_from_phonopy,
    validate_nac_model,
)


def write_qpoints(path: Path) -> None:
    bands = []
    for frequency, direction in zip((1.0, 2.0, 3.0), range(3)):
        vector = [[[0.0, 0.0] for _ in range(3)]]
        vector[0][direction] = [1.0, 0.0]
        bands.append({"frequency": frequency, "eigenvector": vector})
    data = {
        "primitive_cell": {
            "lattice": [[5.0, 0, 0], [0, 5.0, 0], [0, 0, 20.0]],
            "points": [{"symbol": "X", "coordinates": [0, 0, 0.5], "mass": 10.0}],
        },
        "phonon": [{"q-position": [0, 0, 0], "band": bands}],
    }
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")


class InteroperabilityTests(unittest.TestCase):
    def test_sheet_response_is_vacuum_invariant(self):
        target = np.diag([2.5, 1.5, 0.5])
        first_lattice = np.diag([4.0, 4.0, 20.0])
        second_lattice = np.diag([4.0, 4.0, 40.0])
        first_eps = np.eye(3) + 4.0 * np.pi * target / 20.0
        second_eps = np.eye(3) + 4.0 * np.pi * target / 40.0
        first = intrinsic_polarizability_from_supercell(
            first_eps, first_lattice, dimensionality=2
        )
        second = intrinsic_polarizability_from_supercell(
            second_eps, second_lattice, dimensionality=2
        )
        np.testing.assert_allclose(first.values, target)
        np.testing.assert_allclose(second.values, target)
        self.assertEqual(first.unit, "angstrom")

    def test_line_response_uses_nonperiodic_cross_section(self):
        target = np.diag([3.0, 2.0, 1.0])
        lattice = np.diag([10.0, 12.0, 5.0])
        epsilon = np.eye(3) + 4.0 * np.pi * target / 120.0
        result = intrinsic_polarizability_from_supercell(
            epsilon, lattice, dimensionality=1, periodic_axes="z"
        )
        np.testing.assert_allclose(result.values, target)
        self.assertEqual(result.name, "line_polarizability")
        self.assertEqual(result.unit, "angstrom^2")

    def test_low_dimensional_nac_requires_cutoff_model(self):
        with self.assertRaisesRegex(ValueError, "1d-cutoff"):
            validate_nac_model(1, "gonze")
        with self.assertRaisesRegex(ValueError, "2d-cutoff"):
            validate_nac_model(2, "bulk")
        self.assertEqual(validate_nac_model(1, "1d-cutoff"), "1d-cutoff")
        self.assertEqual(validate_nac_model(3, "gonze"), "gonze")

    def test_phonopy_import_preserves_modes_and_nac_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            qpoints = root / "qpoints.yaml"
            write_qpoints(qpoints)
            born = root / "BORN"
            born.write_text(
                "2 0 0 0 2 0 0 0 1.2\n1 0 0 0 1 0 0 0 1\n",
                encoding="utf-8",
            )
            record = response_record_from_phonopy(
                qpoints, born_path=born, dimensionality=2
            )
            self.assertEqual(record.backend, "phonopy")
            self.assertEqual(record.quantity("gamma_frequency").shape, (3,))
            self.assertEqual(record.quantity("gamma_eigenvector_real").shape, (3, 1, 3))
            self.assertEqual(record.quantity("born_effective_charge").shape, (1, 3, 3))
            self.assertEqual(record.quantity("supercell_electronic_dielectric").shape, (3, 3))


if __name__ == "__main__":
    unittest.main()
