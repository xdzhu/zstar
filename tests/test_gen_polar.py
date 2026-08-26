import os
import tempfile
from pathlib import Path
import unittest
from contextlib import redirect_stdout
from io import StringIO

from zstar.gen_polar import (
    _abacus_assets_from_stru,
    _copy_input_sets_to_here,
    gen_input_in_folder,
    print_modified_coordinates,
)
from zstar.deal_polar import _infer_displacement_angstrom
from zstar.phonopy_stru import write_phonopy_compatible_stru
from zstar.stru_analyzer import stru_analyzer


class GenPolarTests(unittest.TestCase):
    def test_dim1_generation_requests_high_precision_charge_cube(self):
        with tempfile.TemporaryDirectory() as tmp:
            previous = Path.cwd()
            try:
                os.chdir(tmp)
                gen_input_in_folder(
                    5,
                    dimension=1,
                    input_mode="pyatb",
                    xc="pbe",
                )
                text = Path("INPUT-scf").read_text(encoding="utf-8")
            finally:
                os.chdir(previous)
            self.assertIn("out_chg             1 10", text)

    def test_deal_infers_generated_half_displacement(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            audit = root / "1.X" / "disp_Angstrom.out"
            audit.parent.mkdir()
            audit.write_text(
                "0.050000 0 0\n-0.050000 0 0\n",
                encoding="utf-8",
            )
            previous = Path.cwd()
            try:
                os.chdir(root)
                displacement, source = _infer_displacement_angstrom()
            finally:
                os.chdir(previous)
            self.assertAlmostEqual(displacement, 0.05)
            self.assertTrue(source.endswith("disp_Angstrom.out"))

    def test_magnetic_initialization_is_preserved_in_displaced_structure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "STRU"
            source.write_text(
                """ATOMIC_SPECIES
Mn 54.938 Mn.upf

NUMERICAL_ORBITAL
Mn.orb

LATTICE_CONSTANT
1.0

LATTICE_VECTORS
5 0 0
0 5 0
0 0 5

ATOMIC_POSITIONS
Direct

Mn
5.0
2
0 0 0 m 1 1 1
0.5 0.5 0.5 m 1 1 1 mag -5
""",
                encoding="utf-8",
            )
            data = stru_analyzer(str(source))
            coordinates = data[5]
            movements = data[6]
            magnetisms = data[7]
            self.assertEqual(magnetisms["Mn"], [5.0, -5.0])

            output = StringIO()
            with redirect_stdout(output):
                print_modified_coordinates(
                    coordinates,
                    [0.01, 0.0, 0.0],
                    1,
                    "x+",
                    movements,
                    magnetisms,
                )
            text = output.getvalue()
            self.assertIn("m 1 1 1 mag 5", text)
            self.assertIn("m 1 1 1 mag -5", text)

    def test_relative_abacus_assets_are_discovered_from_stru(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "Ca.upf").write_text("pseudo", encoding="utf-8")
            (root / "Ca.orb").write_text("orbital", encoding="utf-8")
            (root / "STRU").write_text(
                """ATOMIC_SPECIES
Ca 40.078 Ca.upf upf201

NUMERICAL_ORBITAL
Ca.orb

LATTICE_CONSTANT
1.0
""",
                encoding="utf-8",
            )
            self.assertEqual(
                _abacus_assets_from_stru("STRU", root),
                [str((root / "Ca.upf").resolve()), str((root / "Ca.orb").resolve())],
            )

    def test_input_set_directory_with_spaces_is_copied(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            assets = root / "input assets"
            target = root / "work"
            assets.mkdir()
            target.mkdir()
            (assets / "A.upf").write_text("pseudo", encoding="utf-8")
            previous = Path.cwd()
            try:
                os.chdir(target)
                _copy_input_sets_to_here(str(assets))
            finally:
                os.chdir(previous)
            self.assertEqual((target / "A.upf").read_text(encoding="utf-8"), "pseudo")

    def test_phonopy_stru_adds_missing_movement_fields(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "STRU"
            source.write_text(
                """ATOMIC_SPECIES
In 114.82 In.upf
Se 78.971 Se.upf

NUMERICAL_ORBITAL
In.orb
Se.orb

LATTICE_CONSTANT
1.889726

LATTICE_VECTORS
4.0 0.0 0.0
-2.0 3.5 0.0
0.0 0.0 30.0

ATOMIC_POSITIONS
Direct

In
0
1
0.0 0.0 0.55

Se
0
1
0.3333333333 0.6666666667 0.45
""",
                encoding="utf-8",
            )
            destination = root / "STRU.phonopy"
            write_phonopy_compatible_stru(source, destination)
            text = destination.read_text(encoding="utf-8")
            self.assertIn(
                "0.000000000000 0.000000000000 0.550000000000 m 1 1 1",
                text,
            )
            self.assertIn(
                "0.333333333300 0.666666666700 0.450000000000 m 1 1 1",
                text,
            )

    def test_phonopy_stru_uses_rectangular_rows_for_mixed_magnetism(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "STRU"
            source.write_text(
                """ATOMIC_SPECIES
A 1 A.upf
B 2 B.upf

LATTICE_CONSTANT
1.0

LATTICE_VECTORS
5 0 0
0 5 0
0 0 5

ATOMIC_POSITIONS
Direct

A
0
1
0 0 0 m 1 1 1

B
1
1
0.5 0.5 0.5 m 1 1 1
""",
                encoding="utf-8",
            )
            destination = root / "STRU.phonopy"
            write_phonopy_compatible_stru(source, destination)
            coordinate_rows = [
                line.split()
                for line in destination.read_text(encoding="utf-8").splitlines()
                if line.count(" m ") == 1
            ]
            self.assertEqual([len(row) for row in coordinate_rows], [9, 9])
            self.assertEqual(coordinate_rows[0][-2:], ["mag", "0"])
            self.assertEqual(coordinate_rows[1][-2:], ["mag", "1"])


if __name__ == "__main__":
    unittest.main()
