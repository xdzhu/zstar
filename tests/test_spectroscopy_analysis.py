from pathlib import Path
import tempfile
import unittest

import numpy as np

from zstar.spectroscopy_analysis import (
    calculate_polarized_raman_spectrum,
    optical_constants_from_dielectric,
    read_dielectric_response,
    write_optical_constants,
)


class SpectroscopyAnalysisTests(unittest.TestCase):
    def test_polarized_raman_obeys_parallel_crossed_selection(self):
        frequencies = [500.0, 600.0]
        tensors = np.asarray(
            [
                [[2.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            ]
        )
        parallel = calculate_polarized_raman_spectrum(
            frequencies,
            tensors,
            incident_polarization=(1, 0, 0),
            scattered_polarization=(1, 0, 0),
            points=101,
        )
        crossed = calculate_polarized_raman_spectrum(
            frequencies,
            tensors,
            incident_polarization=(1, 0, 0),
            scattered_polarization=(0, 1, 0),
            points=101,
        )
        np.testing.assert_allclose(parallel.activities, [1.0, 0.0])
        np.testing.assert_allclose(crossed.activities, [0.0, 1.0])

    def test_optical_constants_match_scalar_identities(self):
        frequency = np.asarray([1000.0, 2000.0])
        dielectric = np.zeros((2, 3, 3), dtype=complex)
        dielectric[:] = np.eye(3) * 4.0
        result = optical_constants_from_dielectric(frequency, dielectric)
        np.testing.assert_allclose(result.refractive_index, 2.0)
        np.testing.assert_allclose(result.extinction_coefficient, 0.0)
        np.testing.assert_allclose(result.absorption_cm1, 0.0)
        np.testing.assert_allclose(result.normal_incidence_reflectivity, 1.0 / 9.0)
        np.testing.assert_allclose(result.energy_loss_function, 0.0)

    def test_response_file_roundtrip_and_writer(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            frequency = np.asarray([0.0, 100.0])
            real_tensor = np.tile(np.eye(3) * 2.0, (2, 1, 1))
            imag_tensor = np.tile(np.eye(3) * 0.2, (2, 1, 1))
            np.savetxt(root / "real.dat", np.column_stack([frequency, real_tensor.reshape(2, 9)]))
            np.savetxt(root / "imag.dat", np.column_stack([frequency, imag_tensor.reshape(2, 9)]))
            grid, tensor = read_dielectric_response(root / "real.dat", root / "imag.dat")
            result = optical_constants_from_dielectric(grid, tensor, polarization=(0, 0, 1))
            output = write_optical_constants(root / "optical.dat", result)
            self.assertTrue(output.is_file())
            self.assertEqual(np.loadtxt(output).shape, (2, 10))


if __name__ == "__main__":
    unittest.main()
