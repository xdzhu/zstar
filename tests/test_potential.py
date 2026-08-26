import unittest

import numpy as np

from zstar.potential import AxisProfile, build_arg_parser, estimate_vacuum_sides


class PotentialTests(unittest.TestCase):
    @staticmethod
    def _dipole_corrected_profile() -> tuple[AxisProfile, np.ndarray]:
        coord = np.arange(0.0, 30.0, 0.1)
        values = np.ones_like(coord)
        values[(coord >= 20.0) & (coord < 25.0)] = 2.0
        ramp = (coord >= 25.0) & (coord < 28.0)
        values[ramp] = 2.0 - (coord[ramp] - 25.0) / 3.0
        profile = AxisProfile(
            axis="z",
            index=np.arange(coord.size),
            coord_ang=coord,
            values_ev=values,
            cell_length_ang=30.0,
        )
        return profile, np.asarray([10.0, 20.0])

    def test_vacuum_sides_sample_local_surface_plateaus(self):
        profile, atoms = self._dipole_corrected_profile()
        result = estimate_vacuum_sides(
            profile,
            atoms,
            exclude_distance=2.0,
            plateau_width=0.75,
        )

        self.assertIsNotNone(result)
        self.assertAlmostEqual(result.lower_eV, 1.0, places=12)
        self.assertAlmostEqual(result.upper_eV, 2.0, places=12)
        self.assertAlmostEqual(
            result.delta_upper_minus_lower_eV,
            1.0,
            places=12,
        )
        self.assertLess(result.lower_std_eV, 1.0e-12)
        self.assertLess(result.upper_std_eV, 1.0e-12)
        self.assertGreaterEqual(result.lower_points, 7)
        self.assertGreaterEqual(result.upper_points, 7)

    def test_vacuum_sides_reject_nonpositive_window(self):
        profile, atoms = self._dipole_corrected_profile()
        with self.assertRaisesRegex(ValueError, "must be positive"):
            estimate_vacuum_sides(
                profile,
                atoms,
                exclude_distance=2.0,
                plateau_width=0.0,
            )

    def test_cli_accepts_vacuum_window(self):
        args = build_arg_parser().parse_args(["--vacuum-window", "0.5"])
        self.assertAlmostEqual(args.vacuum_window, 0.5)


if __name__ == "__main__":
    unittest.main()
