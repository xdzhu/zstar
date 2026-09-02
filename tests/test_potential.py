import unittest
from pathlib import Path
import tempfile

import numpy as np

from zstar.potential import (
    AxisProfile,
    DirectionProfile,
    analyze_mirror_symmetry,
    analyze_potential,
    build_arg_parser,
    estimate_vacuum_sides,
)


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

    def test_periodic_mirror_test_distinguishes_symmetric_profile(self):
        coordinate = (np.arange(160, dtype=float) + 0.5) / 160.0
        symmetric_values = (
            np.cos(2.0 * np.pi * coordinate)
            + 0.2 * np.cos(4.0 * np.pi * coordinate)
        )
        profile = DirectionProfile(
            label="a+b",
            safe_label="a_plus_b",
            method="linear",
            lattice_coeffs=np.asarray([1.0, 1.0, 0.0]),
            direction_cart_ang=np.asarray([1.0, 1.0, 0.0]),
            coord_ang=coordinate,
            values_ev=symmetric_values,
            counts=np.ones(160),
            sample_shape=(16, 16),
            smooth_sigma_ang=0.0,
        )
        symmetric = analyze_mirror_symmetry(profile)
        self.assertLess(symmetric.mirror_asymmetry, 1.0e-8)
        self.assertLess(symmetric.odd_rms_ev, 1.0e-8)

        profile.values_ev = symmetric_values + 0.17 * np.sin(
            6.0 * np.pi * coordinate + 0.37
        )
        asymmetric = analyze_mirror_symmetry(profile)
        self.assertGreater(asymmetric.mirror_asymmetry, 1.0e-3)
        self.assertGreater(asymmetric.odd_rms_ev, 1.0e-3)

    def test_mirror_flag_requires_direction(self):
        args = build_arg_parser().parse_args(["--mirror-test"])
        self.assertTrue(args.mirror_test)

    def test_cube_to_axis_plane_direction_and_mirror_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cube = root / "potential.cube"
            grid = np.indices((4, 4, 4), dtype=float)
            values = np.cos(2.0 * np.pi * (grid[0] + grid[1]) / 4.0).reshape(-1)
            lines = [
                "synthetic potential",
                "values in eV",
                "1 0 0 0",
                "4 1 0 0",
                "4 0 1 0",
                "4 0 0 1",
                "1 1.0 0.0 0.0 0.0",
            ]
            for start in range(0, len(values), 6):
                lines.append(" ".join(f"{value:.12g}" for value in values[start:start + 6]))
            cube.write_text("\n".join(lines) + "\n", encoding="utf-8")
            summary = analyze_potential(
                cube=cube,
                outdir=root / "out",
                axes=["z"],
                planes=["xy"],
                plane_average=True,
                directions=["a+b"],
                direction_bins=16,
                direction_samples=(8, 8),
                mirror_test=True,
                value_unit="ev",
                length_unit="angstrom",
                plot=False,
            )
            self.assertTrue(Path(summary["axis_profiles"]["z"]["dat"]).is_file())
            self.assertTrue(Path(summary["plane_maps"]["xy"]["dat"]).is_file())
            direction = summary["direction_profiles"]["a+b:linear"]
            self.assertTrue(Path(direction["dat"]).is_file())
            self.assertTrue(Path(direction["mirror"]["dat"]).is_file())
            self.assertLess(direction["mirror"]["mirror_asymmetry"], 1.0e-6)


if __name__ == "__main__":
    unittest.main()
