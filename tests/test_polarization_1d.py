import json
import math
import tempfile
from pathlib import Path
import unittest

import numpy as np

from zstar.polarization_1d import (
    calculate_hybrid_1d_born,
    integrate_transverse_dipole,
    write_reference_1d_polarization,
)
from zstar.polarization_2d import (
    BOHR_M,
    ELEMENTARY_CHARGE,
    _periodic_weighted_center,
)


def write_wire_cube(path: Path, dipole_x_e_bohr: float, dipole_y_e_bohr: float):
    nx, ny, nz = 8, 8, 2
    ionic = np.asarray([4.0, 4.0, 1.0])
    electron = ionic - np.asarray([dipole_x_e_bohr, dipole_y_e_bohr, 0.0])
    x0 = math.floor(electron[0])
    y0 = math.floor(electron[1])
    wx = electron[0] - x0
    wy = electron[1] - y0
    density = np.zeros((nx, ny, nz))
    for dx, x_weight in ((0, 1.0 - wx), (1, wx)):
        for dy, y_weight in ((0, 1.0 - wy), (1, wy)):
            density[(x0 + dx) % nx, (y0 + dy) % ny, 1] += x_weight * y_weight
    lines = [
        "wire cube",
        "density in e/bohr^3",
        "1 0.0 0.0 0.0",
        f"{nx} 1.0 0.0 0.0",
        f"{ny} 0.0 1.0 0.0",
        f"{nz} 0.0 0.0 1.0",
        "1 1.0 4.0 4.0 1.0",
    ]
    flat = density.reshape(-1)
    for start in range(0, len(flat), 6):
        lines.append(" ".join(f"{value:.12e}" for value in flat[start : start + 6]))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_boundary_spanning_wire_cube(path: Path):
    density = np.zeros((8, 1, 1))
    density[0, 0, 0] = 2.0
    lines = [
        "boundary-spanning wire cube",
        "density in e/bohr^3",
        "2 0.0 0.0 0.0",
        "8 1.0 0.0 0.0",
        "1 0.0 1.0 0.0",
        "1 0.0 0.0 1.0",
        "1 1.0 0.2 0.0 0.0",
        "1 1.0 7.8 0.0 0.0",
    ]
    lines.extend(f"{value:.12e}" for value in density.reshape(-1))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_polarization(path: Path, values):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            f"The calculated polarization direction is in {axis}, "
            f"P = {value:.12e} (mod 1.000000000000e+02) C/m^2."
            for axis, value in zip("abc", values)
        )
        + "\n",
        encoding="utf-8",
    )


class Polarization1DTests(unittest.TestCase):
    def test_periodic_center_keeps_boundary_spanning_wire_contiguous(self):
        center = _periodic_weighted_center(
            np.asarray([0.2, 7.8]),
            np.asarray([1.0, 1.0]),
            origin=0.0,
            period=8.0,
        )
        distance_to_boundary = min(abs(center), abs(center - 8.0))
        self.assertLess(distance_to_boundary, 1.0e-12)

    def test_transverse_dipole_does_not_cut_boundary_spanning_wire(self):
        with tempfile.TemporaryDirectory() as tmp:
            cube = Path(tmp) / "SPIN1_CHG.cube"
            write_boundary_spanning_wire_cube(cube)
            result = integrate_transverse_dipole(cube, "x")
            self.assertAlmostEqual(result.dipole_e_bohr, 0.0, places=12)

    def test_transverse_cube_dipoles(self):
        with tempfile.TemporaryDirectory() as tmp:
            cube = Path(tmp) / "SPIN1_CHG.cube"
            write_wire_cube(cube, 0.75, -0.35)
            self.assertAlmostEqual(
                integrate_transverse_dipole(cube, "x").dipole_e_bohr,
                0.75,
                places=9,
            )
            self.assertAlmostEqual(
                integrate_transverse_dipole(cube, "y").dipole_e_bohr,
                -0.35,
                places=9,
            )
            overridden = integrate_transverse_dipole(
                cube,
                "x",
                unwrap_center_bohr=1.25,
            )
            self.assertAlmostEqual(overridden.unwrap_center_bohr, 1.25)

    def test_hybrid_tensor_uses_cube_xy_and_berry_z(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference = root / "0.no-move"
            atom = root / "1.X"
            write_wire_cube(reference / "OUT.TEST" / "SPIN1_CHG.cube", 0.0, 0.0)
            write_polarization(
                reference / "pyatb" / "Out" / "Polarization" / "polarization.dat",
                [0.0, 0.0, 0.0],
            )
            lattice_angstrom = np.diag([8.0, 8.0, 2.0])
            input_json = reference / "pyatb" / "Out" / "input.json"
            input_json.parent.mkdir(parents=True, exist_ok=True)
            input_json.write_text(
                json.dumps(
                    {"LATTICE": {"lattice_vector": lattice_angstrom.tolist()}}
                ),
                encoding="utf-8",
            )

            reference_summary = write_reference_1d_polarization(reference)
            reference_data = json.loads(reference_summary.read_text(encoding="utf-8"))
            self.assertEqual(reference_data["periodic_axis"], "z")
            self.assertAlmostEqual(reference_data["cross_section_angstrom2"], 64.0)

            target = np.asarray(
                [
                    [2.0, 0.2, 0.4],
                    [0.3, 3.0, 0.5],
                    [0.6, 0.7, 4.0],
                ]
            )
            displacement_angstrom = 0.01
            displacement_m = displacement_angstrom * 1.0e-10
            displacement_bohr = displacement_m / BOHR_M
            volume_m3 = abs(np.linalg.det(lattice_angstrom)) * 1.0e-30
            for beta, direction in enumerate("xyz"):
                plus = atom / f"{direction}+"
                minus = atom / f"{direction}-"
                write_wire_cube(
                    plus / "OUT.TEST" / "SPIN1_CHG.cube",
                    target[beta, 0] * displacement_bohr,
                    target[beta, 1] * displacement_bohr,
                )
                write_wire_cube(
                    minus / "OUT.TEST" / "SPIN1_CHG.cube",
                    -target[beta, 0] * displacement_bohr,
                    -target[beta, 1] * displacement_bohr,
                )
                pz = target[beta, 2] * displacement_m * ELEMENTARY_CHARGE / volume_m3
                write_polarization(
                    plus / "pyatb" / "Out" / "Polarization" / "polarization.dat",
                    [0.0, 0.0, pz],
                )
                write_polarization(
                    minus / "pyatb" / "Out" / "Polarization" / "polarization.dat",
                    [0.0, 0.0, -pz],
                )

            result = calculate_hybrid_1d_born(
                atom,
                reference,
                method="central",
                displacement_angstrom=displacement_angstrom,
            )
            np.testing.assert_allclose(result.tensor, target, atol=1.0e-8)
            self.assertEqual(result.periodic_axis, "z")


if __name__ == "__main__":
    unittest.main()
