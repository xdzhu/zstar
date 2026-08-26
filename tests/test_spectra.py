import json
import tempfile
from pathlib import Path
import unittest

import numpy as np
import yaml

from zstar.spectra import (
    BornData,
    GammaModes,
    calculate_ir_spectrum,
    calculate_molecular_ir_spectrum,
    calculate_raman_spectrum,
    collect_molecular_dipole_derivatives,
    collect_raman_tensors,
    load_gamma_modes,
    mode_effective_charges,
    prepare_raman_displacements,
    read_born_data,
    read_pyatb_polarization,
    write_ir_outputs,
    write_molecular_ir_outputs,
    write_raman_outputs,
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


def write_polarization(path: Path, values, quanta=(100.0, 100.0, 100.0)):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            f"The calculated polarization direction is in {axis}, "
            f"P = {value:.12e} (mod {quantum:.12e}) C/m^2."
            for axis, value, quantum in zip("abc", values, quanta)
        )
        + "\n",
        encoding="utf-8",
    )


class SpectraTests(unittest.TestCase):
    def test_mode_effective_charge_uses_displacement_rows(self):
        modes = GammaModes(
            frequencies_thz=np.asarray([10.0]),
            eigenvectors=np.asarray([[[1.0, 0.0, 0.0]]]),
            masses_amu=np.asarray([1.0]),
            lattice_angstrom=np.eye(3),
            symbols=("X",),
            positions_fractional=np.zeros((1, 3)),
        )
        canonical = np.asarray(
            [[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]]
        )

        effective = mode_effective_charges(modes, canonical)

        np.testing.assert_allclose(effective, [[1.0, 2.0, 3.0]])

    def test_load_gamma_modes_uses_companion_phonopy_yaml(self):
        with tempfile.TemporaryDirectory() as tmp:
            qpoints = Path(tmp) / "qpoints.yaml"
            write_split_qpoints(qpoints)
            modes = load_gamma_modes(qpoints)
            self.assertEqual(modes.symbols, ("X",))
            self.assertAlmostEqual(modes.cell_height_angstrom, 20.0)

    def test_load_gamma_modes_converts_companion_bohr_lattice(self):
        with tempfile.TemporaryDirectory() as tmp:
            qpoints = Path(tmp) / "qpoints.yaml"
            write_qpoints(qpoints)
            (qpoints.parent / "phonopy.yaml").write_text(
                yaml.safe_dump(
                    {"physical_unit": {"length": "au"}}, sort_keys=False
                ),
                encoding="utf-8",
            )
            modes = load_gamma_modes(qpoints)
            self.assertAlmostEqual(
                modes.cell_height_angstrom, 20.0 * 0.529177210903
            )

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

    def test_ir_1d_line_response_is_vacuum_invariant(self):
        target = np.diag([3.0, 2.0, 1.0])
        responses = []
        for vacuum in (10.0, 20.0):
            lattice = np.diag([vacuum, vacuum, 5.0])
            modes = GammaModes(
                frequencies_thz=np.asarray([10.0]),
                eigenvectors=np.asarray([[[1.0, 0.0, 0.0]]]),
                masses_amu=np.asarray([1.0]),
                lattice_angstrom=lattice,
                symbols=("X",),
                positions_fractional=np.asarray([[0.0, 0.0, 0.0]]),
            )
            born = BornData(
                tensors=np.zeros((1, 3, 3)),
                electronic_dielectric=np.eye(3) + target / (vacuum * vacuum),
                source="test",
            )
            result = calculate_ir_spectrum(
                modes,
                born,
                dimensionality=1,
                acoustic_cutoff_cm1=0.0,
                points=11,
            )
            responses.append(result.response_real[0])
            self.assertEqual(
                result.response_kind,
                "1D line polarizability (Angstrom^2; SI-reduced)",
            )
            self.assertEqual(result.response_unit, "angstrom^2")
            self.assertIn("nonperiodic_cross_section", result.response_convention)
        np.testing.assert_allclose(responses[0], target)
        np.testing.assert_allclose(responses[1], target)

    def test_read_and_collect_molecular_dipole_derivatives(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "modes"
            plus = root / "mode-0001" / "plus"
            minus = root / "mode-0001" / "minus"
            write_polarization(
                plus / "pyatb-polar" / "Out" / "Polarization" / "polarization.dat",
                (0.01, 49.99, 0.03),
            )
            write_polarization(
                minus / "pyatb-polar" / "Out" / "Polarization" / "polarization.dat",
                (-0.01, -49.99, -0.03),
            )
            (root / "raman_manifest.json").write_text(
                json.dumps(
                    {
                        "amplitude_A_sqrt_amu": 0.02,
                        "modes": [
                            {
                                "mode": 1,
                                "plus": "/stale/location/plus",
                                "minus": "/stale/location/minus",
                            }
                        ],
                    },
                    indent=2,
                ),
                encoding="utf-8",
            )
            values, quanta, source = read_pyatb_polarization(
                plus / "pyatb-polar"
            )
            np.testing.assert_allclose(values, [0.01, 49.99, 0.03])
            np.testing.assert_allclose(quanta, [100.0, 100.0, 100.0])
            self.assertEqual(source.name, "polarization.dat")

            numbers, derivatives, kind = collect_molecular_dipole_derivatives(
                root,
                cell_volume_angstrom3=500.0,
            )
            expected_delta = np.asarray([0.02, -0.02, 0.06])
            expected = expected_delta * 500.0e-30 / (0.04 * 3.33564e-30)
            np.testing.assert_array_equal(numbers, [1])
            np.testing.assert_allclose(derivatives[0], expected)
            self.assertIn("molecular dipole derivative", kind)
            self.assertTrue((root / "molecular_ir_derivatives.json").is_file())

            lattice = np.asarray(
                [[2.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 3.0]]
            )
            _, skew_derivatives, _ = collect_molecular_dipole_derivatives(
                root,
                cell_volume_angstrom3=500.0,
                cell_lattice_angstrom=lattice,
            )
            basis_to_cartesian = lattice / np.linalg.norm(
                lattice, axis=1
            )[:, None]
            expected_skew = (
                expected_delta
                @ basis_to_cartesian
                * 500.0e-30
                / (0.04 * 3.33564e-30)
            )
            np.testing.assert_allclose(skew_derivatives[0], expected_skew)

    def test_read_pyatb_reconstructs_small_polarization_from_phases(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "polarization.dat"
            path.write_text(
                "The Ionic Phase      : 0.002000 0.000000 0.000000\n"
                "The Electronic Phase : -0.002011 0.000000 0.000000\n"
                "The calculated polarization direction is in a, "
                "P = -0.000000 (mod 0.040024) C/m^2.\n"
                "The calculated polarization direction is in b, "
                "P = 0.000000 (mod 0.040024) C/m^2.\n"
                "The calculated polarization direction is in c, "
                "P = 0.000000 (mod 0.040024) C/m^2.\n",
                encoding="utf-8",
            )
            values, quanta, _ = read_pyatb_polarization(path)
            np.testing.assert_allclose(quanta, [0.040024] * 3)
            np.testing.assert_allclose(values, [-4.40264e-7, 0.0, 0.0])

            precise = path.read_text(encoding="utf-8").replace(
                "P = -0.000000 (mod", "P = -0.000000410000 (mod", 1
            )
            path.write_text(precise, encoding="utf-8")
            precise_values, _, _ = read_pyatb_polarization(path)
            self.assertEqual(precise_values[0], -4.1e-7)

            bulk_like = path.read_text(encoding="utf-8").replace(
                "-0.000000410000", "-0.000000", 1
            ).replace("0.040024", "1.000000")
            path.write_text(bulk_like, encoding="utf-8")
            bulk_values, _, _ = read_pyatb_polarization(path)
            self.assertEqual(bulk_values[0], 0.0)

    def test_calculate_molecular_ir_spectrum(self):
        with tempfile.TemporaryDirectory() as tmp:
            qpoints = Path(tmp) / "qpoints.yaml"
            write_qpoints(qpoints)
            modes = load_gamma_modes(qpoints)
            result = calculate_molecular_ir_spectrum(
                modes,
                [1, 2],
                np.asarray([[0.0, 0.0, 0.0], [1.0, 2.0, 2.0]]),
                points=101,
            )
            np.testing.assert_allclose(result.activities, [0.0, 9.0])
            np.testing.assert_allclose(result.normalized_activities, [0.0, 1.0])
            self.assertAlmostEqual(float(np.max(result.spectrum)), 1.0)
            summary = write_molecular_ir_outputs(
                Path(tmp) / "molecular_ir", result, plot=False
            )
            self.assertEqual(summary["dimensionality"], 0)
            self.assertTrue(
                (Path(tmp) / "molecular_ir" / "ir_modes.csv").is_file()
            )

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

    def test_spectrum_writers_export_raster_and_vector_plots(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            qpoints = root / "qpoints.yaml"
            write_qpoints(qpoints)
            modes = load_gamma_modes(qpoints)
            born = BornData(
                tensors=np.asarray([np.eye(3)]),
                electronic_dielectric=np.eye(3),
                source="test",
            )
            ir = calculate_ir_spectrum(
                modes, born, acoustic_cutoff_cm1=0.0, points=101
            )
            ir_summary = write_ir_outputs(root / "ir", ir)
            raman = calculate_raman_spectrum(
                modes, [1], np.asarray([np.eye(3)]), points=101
            )
            raman_summary = write_raman_outputs(root / "raman", raman)

            for summary in (ir_summary, raman_summary):
                for key in ("plot", "plot_pdf", "plot_svg"):
                    self.assertTrue(Path(summary["files"][key]).is_file())

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

            numbers, tensors, kind = collect_raman_tensors(
                root / "raman",
                dimensionality=0,
                cell_volume_angstrom3=500.0,
            )
            np.testing.assert_array_equal(numbers, [1])
            np.testing.assert_allclose(
                tensors[0], np.eye(3) * 1000.0 / (4.0 * np.pi)
            )
            self.assertIn("molecular polarizability derivative", kind)

            numbers, tensors, kind = collect_raman_tensors(
                root / "raman",
                dimensionality=1,
                cell_cross_section_angstrom2=25.0,
            )
            np.testing.assert_array_equal(numbers, [1])
            np.testing.assert_allclose(
                tensors[0], np.eye(3) * 50.0 / (4.0 * np.pi)
            )
            self.assertIn("1D line polarizability derivative", kind)


if __name__ == "__main__":
    unittest.main()
