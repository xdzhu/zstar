import json
import math
import tempfile
from pathlib import Path
import unittest

import numpy as np

from zstar.polarization_2d import (
    BOHR_M,
    ELEMENTARY_CHARGE,
    calculate_hybrid_2d_born,
    compare_slab_charge_profiles,
    integrate_slab_dipole,
    write_slab_charge_difference,
)


def write_cube(
    path: Path,
    dipole_e_bohr: float,
    *,
    first_step=(1.0, 0.0, 0.0),
    second_step=(0.0, 1.0, 0.0),
    third_step=(0.0, 0.0, 1.0),
):
    nx, ny, nz = 2, 2, 8
    ionic_z = 4.0
    electron_center = ionic_z - dipole_e_bohr
    lower = math.floor(electron_center)
    upper = lower + 1
    upper_weight = electron_center - lower
    density = np.zeros((nx, ny, nz))
    density[0, 0, lower % nz] += 1.0 - upper_weight
    density[0, 0, upper % nz] += upper_weight
    lines = [
        "test cube",
        "density in e/bohr^3",
        "1 0.0 0.0 0.0",
        f"{nx} {' '.join(str(value) for value in first_step)}",
        f"{ny} {' '.join(str(value) for value in second_step)}",
        f"{nz} {' '.join(str(value) for value in third_step)}",
        f"1 1.0 0.0 0.0 {ionic_z}",
    ]
    flat = density.reshape(-1)
    for start in range(0, len(flat), 6):
        lines.append(" ".join(f"{value:.12e}" for value in flat[start : start + 6]))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_polarization(path: Path, values):
    path.parent.mkdir(parents=True, exist_ok=True)
    axes = "abc"
    path.write_text(
        "\n".join(
            f"The calculated polarization direction is in {axis}, "
            f"P = {value:.12e} (mod 1.000000000000e+02) C/m^2."
            for axis, value in zip(axes, values)
        )
        + "\n",
        encoding="utf-8",
    )


class Polarization2DTests(unittest.TestCase):
    def test_real_space_dipole_integral(self):
        with tempfile.TemporaryDirectory() as tmp:
            cube = Path(tmp) / "SPIN1_CHG.cube"
            write_cube(cube, 0.75)
            result = integrate_slab_dipole(cube)
            self.assertAlmostEqual(result.dipole_e_bohr, 0.75, places=9)
            expected = 0.75 / 4.0 * ELEMENTARY_CHARGE / BOHR_M
            self.assertAlmostEqual(result.polarization_C_per_m, expected)

    def test_slab_charge_profile_recovers_effective_charge(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            displacement_angstrom = 0.01
            displacement_bohr = displacement_angstrom * 1.0e-10 / BOHR_M
            reference = root / "reference.cube"
            displaced = root / "displaced.cube"
            write_cube(reference, 0.0)
            write_cube(displaced, 2.5 * displacement_bohr)

            result = compare_slab_charge_profiles(
                reference,
                displaced,
                displacement_angstrom=displacement_angstrom,
            )
            self.assertAlmostEqual(result.effective_charge_e, 2.5, places=8)
            self.assertLess(
                abs(result.diagnostics["profile_closure_error_e_angstrom"]),
                1.0e-10,
            )
            summary = write_slab_charge_difference(
                root / "profile", result, plot=False
            )
            self.assertTrue(
                (root / "profile" / summary["files"]["profile"]).is_file()
            )
            self.assertTrue(
                (root / "profile" / summary["files"]["summary"]).is_file()
            )
            with self.assertRaisesRegex(ValueError, "finite and non-zero"):
                compare_slab_charge_profiles(
                    reference,
                    displaced,
                    displacement_angstrom=0.0,
                )

    def test_hybrid_tensor_uses_berry_xy_and_cube_z(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference = root / "0.no-move"
            atom = root / "1.X"
            write_cube(reference / "OUT.TEST" / "SPIN1_CHG.cube", 0.0)
            write_polarization(
                reference / "pyatb" / "Out" / "Polarization" / "polarization.dat",
                [0.0, 0.0, 0.0],
            )
            lattice_ang = np.diag([2.0, 2.0, 4.0])
            input_json = reference / "pyatb" / "Out" / "input.json"
            input_json.parent.mkdir(parents=True, exist_ok=True)
            input_json.write_text(
                json.dumps({"LATTICE": {"lattice_vector": lattice_ang.tolist()}}),
                encoding="utf-8",
            )

            displacement_ang = 0.01
            displacement_m = displacement_ang * 1.0e-10
            volume_m3 = abs(np.linalg.det(lattice_ang)) * 1.0e-30
            targets = np.asarray(
                [
                    [2.0, 0.2, 0.0],
                    [0.3, 3.0, 0.0],
                    [0.4, 0.5, 0.0],
                ]
            )
            z_targets = [0.6, 0.7, 4.0]
            displacement_bohr = displacement_ang * 1.0e-10 / BOHR_M
            for beta, direction in enumerate("xyz"):
                delta_p = (
                    targets[beta, :2]
                    * (2.0 * displacement_m)
                    * ELEMENTARY_CHARGE
                    / volume_m3
                )
                plus_values = [delta_p[0] / 2.0, delta_p[1] / 2.0, 0.0]
                minus_values = [-delta_p[0] / 2.0, -delta_p[1] / 2.0, 0.0]
                plus = atom / f"{direction}+"
                minus = atom / f"{direction}-"
                write_polarization(
                    plus / "pyatb" / "Out" / "Polarization" / "polarization.dat",
                    plus_values,
                )
                write_polarization(
                    minus / "pyatb" / "Out" / "Polarization" / "polarization.dat",
                    minus_values,
                )
                half_delta = z_targets[beta] * displacement_bohr
                write_cube(plus / "OUT.TEST" / "SPIN1_CHG.cube", half_delta)
                write_cube(minus / "OUT.TEST" / "SPIN1_CHG.cube", -half_delta)

            result = calculate_hybrid_2d_born(
                atom,
                reference,
                method="central",
                displacement_angstrom=displacement_ang,
            )
            expected = np.asarray(
                [
                    [2.0, 0.3, 0.4],
                    [0.2, 3.0, 0.5],
                    [0.6, 0.7, 4.0],
                ]
            )
            np.testing.assert_allclose(result.tensor, expected, atol=1.0e-8)

            y_only = calculate_hybrid_2d_born(
                atom,
                reference,
                method="central",
                displacement_angstrom=displacement_ang,
                directions=("y",),
            )
            expected_y = np.zeros((3, 3))
            expected_y[:, 1] = expected[:, 1]
            np.testing.assert_allclose(y_only.tensor, expected_y, atol=1.0e-8)

    def test_hybrid_tensor_rejects_a_tilted_slab_normal(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference = root / "0.no-move"
            write_cube(
                reference / "OUT.TEST" / "SPIN1_CHG.cube",
                0.0,
                first_step=(1.0, 0.0, 1.0),
            )
            write_polarization(
                reference / "pyatb" / "Out" / "Polarization" / "polarization.dat",
                [0.0, 0.0, 0.0],
            )
            input_json = reference / "pyatb" / "Out" / "input.json"
            input_json.parent.mkdir(parents=True, exist_ok=True)
            input_json.write_text(
                json.dumps(
                    {
                        "LATTICE": {
                            "lattice_vector": [
                                [2.0, 0.0, 2.0],
                                [0.0, 2.0, 0.0],
                                [0.0, 0.0, 8.0],
                            ]
                        }
                    }
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "Cartesian z"):
                calculate_hybrid_2d_born(
                    root / "1.X",
                    reference,
                    method="forward",
                    directions=("x",),
                )


if __name__ == "__main__":
    unittest.main()
