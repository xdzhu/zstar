import tempfile
from pathlib import Path
import unittest

import numpy as np
import yaml

from zstar.spectra import (
    BornData,
    calculate_ir_spectrum,
    calculate_raman_spectrum,
    collect_raman_tensors,
    load_gamma_modes,
    prepare_raman_displacements,
    read_born_data,
)


def write_qpoints(path: Path):
    bands = []
    for frequency, direction in zip((1.0, 2.0, 3.0), range(3)):
        vector = [[[0.0, 0.0] for _ in range(3)]]
        vector[0][direction] = [1.0, 0.0]
        bands.append({"frequency": frequency, "eigenvector": vector})
    data = {
        "primitive_cell": {
            "lattice": [[5.0, 0, 0], [0, 5.0, 0], [0, 0, 20.0]],
            "points": [
                {"symbol": "X", "coordinates": [0.0, 0.0, 0.5], "mass": 1.0}
            ],
        },
        "phonon": [{"q-position": [0, 0, 0], "band": bands}],
    }
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")


def write_split_qpoints(path: Path):
    write_qpoints(path)
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    cell = data.pop("primitive_cell")
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    (path.parent / "phonopy.yaml").write_text(
        yaml.safe_dump({"primitive_cell": cell}, sort_keys=False),
        encoding="utf-8",
    )


class SpectraTests(unittest.TestCase):
    def test_load_gamma_modes_uses_companion_phonopy_yaml(self):
        with tempfile.TemporaryDirectory() as tmp:
            qpoints = Path(tmp) / "qpoints.yaml"
            write_split_qpoints(qpoints)
            modes = load_gamma_modes(qpoints)
            self.assertEqual(modes.symbols, ("X",))
            self.assertAlmostEqual(modes.cell_height_angstrom, 20.0)

    def test_ir_effective_charge_and_2d_response(self):
        with tempfile.TemporaryDirectory() as tmp:
            qpoints = Path(tmp) / "qpoints.yaml"
            write_qpoints(qpoints)
            modes = load_gamma_modes(qpoints)
            born = BornData(
                tensors=np.asarray([np.diag([2.0, 3.0, 4.0])]),
                electronic_dielectric=np.diag([2.0, 2.0, 1.5]),
                source="test",
            )
            result = calculate_ir_spectrum(
                modes,
                born,
                dimensionality=2,
                acoustic_cutoff_cm1=0.0,
                points=101,
            )
            np.testing.assert_allclose(
                result.effective_charges, np.diag([2.0, 3.0, 4.0])
            )
            self.assertEqual(result.response_kind, "2D sheet polarizability (Angstrom)")
            self.assertTrue(np.all(np.isfinite(result.response_real)))

    def test_reduced_phonopy_born_auto_loads_full_sibling(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "BORN").write_text(
                "1.0\n"
                "2 0 0 0 2 0 0 0 2\n"
                "1 0 0 0 1 0 0 0 1\n",
                encoding="utf-8",
            )
            (root / "Z-BORN-all.out").write_text(
                "1 X 1 0 0 0 1 0 0 0 1\n"
                "2 X -1 0 0 0 -1 0 0 0 -1\n",
                encoding="utf-8",
            )
            born = read_born_data(root / "BORN", natoms=2)
            self.assertEqual(born.tensors.shape, (2, 3, 3))
            np.testing.assert_allclose(born.electronic_dielectric, np.eye(3) * 2)
            self.assertIn("Z-BORN-all.out", born.source)

    def test_raman_placzek_isotropic_tensor(self):
        with tempfile.TemporaryDirectory() as tmp:
            qpoints = Path(tmp) / "qpoints.yaml"
            write_qpoints(qpoints)
            modes = load_gamma_modes(qpoints)
            result = calculate_raman_spectrum(
                modes,
                [1],
                np.asarray([np.eye(3)]),
                points=101,
            )
            self.assertAlmostEqual(result.activities[0], 1.0)
            self.assertAlmostEqual(result.depolarization_ratios[0], 0.0)
            self.assertAlmostEqual(float(np.max(result.spectrum)), 1.0)

    def test_prepare_raman_mode_pair(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            qpoints = root / "qpoints.yaml"
            write_qpoints(qpoints)
            stru = root / "STRU"
            stru.write_text(
                """ATOMIC_SPECIES
X 1.0 X.upf

NUMERICAL_ORBITAL
X.orb

LATTICE_CONSTANT
1.889726125

LATTICE_VECTORS
5 0 0
0 5 0
0 0 20

ATOMIC_POSITIONS
Direct

X
0.0
1
0.0 0.0 0.5 m 1 1 1
""",
                encoding="utf-8",
            )
            modes = load_gamma_modes(qpoints)
            manifest = prepare_raman_displacements(
                stru,
                modes,
                root / "raman",
                mode_numbers=[1],
            )
            self.assertEqual(len(manifest["modes"]), 1)
            plus = root / "raman" / "mode-0001" / "plus" / "STRU"
            minus = root / "raman" / "mode-0001" / "minus" / "STRU"
            self.assertTrue(plus.is_file())
            self.assertNotEqual(
                plus.read_text(encoding="utf-8"),
                minus.read_text(encoding="utf-8"),
            )

    def test_collect_raman_tensors_from_workflow_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            qpoints = root / "qpoints.yaml"
            write_qpoints(qpoints)
            stru = root / "STRU"
            stru.write_text(
                """ATOMIC_SPECIES
X 1.0 X.upf

NUMERICAL_ORBITAL
X.orb

LATTICE_CONSTANT
1.889726125

LATTICE_VECTORS
5 0 0
0 5 0
0 0 20

ATOMIC_POSITIONS
Direct

X
0.0
1
0.0 0.0 0.5 m 1 1 1
""",
                encoding="utf-8",
            )
            modes = load_gamma_modes(qpoints)
            prepare_raman_displacements(
                stru,
                modes,
                root / "raman",
                mode_numbers=[1],
                amplitude=0.02,
            )
            for sign, scale in (("plus", 1.04), ("minus", 0.96)):
                output = (
                    root
                    / "raman"
                    / "mode-0001"
                    / sign
                    / "pyatb"
                    / "Out"
                    / "Optical_Conductivity"
                )
                output.mkdir(parents=True)
                tensor = np.eye(3) * scale
                (output / "static_dielectric_function.dat").write_text(
                    " ".join(str(value) for value in tensor.reshape(-1)) + "\n",
                    encoding="utf-8",
                )
            numbers, tensors, kind = collect_raman_tensors(root / "raman")
            np.testing.assert_array_equal(numbers, [1])
            np.testing.assert_allclose(tensors[0], np.eye(3) * 2.0)
            self.assertEqual(kind, "dielectric tensor derivative")


if __name__ == "__main__":
    unittest.main()
