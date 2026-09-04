import os
import tempfile
from pathlib import Path
import unittest
from contextlib import redirect_stdout
from io import StringIO
from unittest.mock import patch

from zstar.cli import zstar_cli
from zstar.gen_polar import (
    _abacus_assets_from_stru,
    _copy_input_sets_to_here,
    gen_input_in_folder,
    gen_polar,
    print_modified_coordinates,
)
from zstar.deal_polar import _infer_displacement_angstrom
from zstar.phonopy_stru import write_phonopy_compatible_stru
from zstar.stru_analyzer import stru_analyzer
from zstar.symmetry_reduction import reduce_abacus_atoms, write_reduction_report


class GenPolarTests(unittest.TestCase):
    def test_gen_polar_uses_spglib_representatives_for_task_directories(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            structure = root / "STRU"
            structure.write_text(
                """ATOMIC_SPECIES
Ba 137.327 Ba.upf
Ti 47.867 Ti.upf
O 15.999 O.upf

LATTICE_CONSTANT
1.0

LATTICE_VECTORS
4.0 0.0 0.0
0.0 4.0 0.0
0.0 0.0 4.0

ATOMIC_POSITIONS
Direct

Ba
0
1
0.0 0.0 0.0

Ti
0
1
0.5 0.5 0.5

O
0
3
0.5 0.5 0.0
0.5 0.0 0.5
0.0 0.5 0.5
""",
                encoding="utf-8",
            )
            previous = Path.cwd()
            try:
                os.chdir(root)
                with redirect_stdout(StringIO()):
                    gen_polar(
                        f_stru="STRU",
                        force_delete=True,
                        dimension=3,
                        input_mode="custom",
                        extract_starred_atoms_only=True,
                        method="forward",
                    )
            finally:
                os.chdir(previous)
            self.assertTrue((root / "0.no-move").is_dir())
            self.assertEqual(
                sorted(path.name for path in root.iterdir() if path.is_dir()),
                ["0.no-move", "1.Ba", "2.Ti", "3.O"],
            )
            self.assertIn(
                "representatives = 1 2 3",
                (root / "reduced_atom.out").read_text(encoding="utf-8"),
            )

    def test_spglib_reduction_returns_stable_representatives(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            structure = root / "STRU"
            structure.write_text(
                """ATOMIC_SPECIES
Ba 137.327 Ba.upf
Ti 47.867 Ti.upf
O 15.999 O.upf

LATTICE_CONSTANT
1.0

LATTICE_VECTORS
4.0 0.0 0.0
0.0 4.0 0.0
0.0 0.0 4.0

ATOMIC_POSITIONS
Direct

Ba
0
1
0.0 0.0 0.0

Ti
0
1
0.5 0.5 0.5

O
0
3
0.5 0.5 0.0
0.5 0.0 0.5
0.0 0.5 0.5
""",
                encoding="utf-8",
            )
            result = reduce_abacus_atoms(structure, dimensionality=3)
            self.assertEqual(result.representatives, (1, 2, 3))
            self.assertEqual(result.space_group, "Pm-3m")
            report = write_reduction_report(root / "reduced_atom.out", result)
            self.assertIn("representatives = 1 2 3", report.read_text(encoding="utf-8"))

    def test_molecular_reduction_does_not_infer_periodic_vacuum_symmetry(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            structure = root / "STRU"
            structure.write_text(
                """ATOMIC_SPECIES
O 15.999 O.upf
H 1.008 H.upf

LATTICE_CONSTANT
1.0

LATTICE_VECTORS
20.0 0.0 0.0
0.0 20.0 0.0
0.0 0.0 20.0

ATOMIC_POSITIONS
Cartesian

O
0
1
10.0 10.0 10.0

H
0
2
10.7 10.0 10.0
9.8 10.6 10.0
""",
                encoding="utf-8",
            )
            result = reduce_abacus_atoms(structure, dimensionality=0)
            self.assertEqual(result.engine, "none-molecular")
            self.assertEqual(result.representatives, (1, 2, 3))

    def test_symmetry_failure_explains_all_atom_fallback(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            structure = root / "STRU"
            structure.write_text(
                """ATOMIC_SPECIES
A 1.0 A.upf

LATTICE_CONSTANT
1.0

LATTICE_VECTORS
4 0 0
0 4 0
0 0 4

ATOMIC_POSITIONS
Direct

A
0
1
0 0 0
""",
                encoding="utf-8",
            )
            with patch(
                "zstar.symmetry_reduction.spglib.get_symmetry_dataset",
                return_value=None,
            ):
                with self.assertRaisesRegex(RuntimeError, "--all.*--symmprec"):
                    reduce_abacus_atoms(structure, dimensionality=3)

    def test_cli_preserves_input_functional_without_explicit_xc(self):
        with patch("zstar.gen_polar.gen_polar") as run_gen:
            zstar_cli(["gen", "--ensemble", "cartesian", "-i", "INPUT.seed", "--stru", "STRU"])
        self.assertIsNone(run_gen.call_args.kwargs["xc"])

    def test_cli_explicit_xc_overrides_input_functional(self):
        with patch("zstar.gen_polar.gen_polar") as run_gen:
            zstar_cli(
                [
                    "gen", "--ensemble", "cartesian", "-i", "INPUT.seed", "--stru", "STRU",
                    "--xc", "pbesol",
                ]
            )
        self.assertEqual(run_gen.call_args.kwargs["xc"], "pbesol")

    def test_cli_defaults_generated_input_to_pbe(self):
        with patch("zstar.gen_polar.gen_polar") as run_gen:
            zstar_cli(["gen", "--ensemble", "cartesian", "--stru", "STRU"])
        self.assertEqual(run_gen.call_args.kwargs["xc"], "pbe")

    def test_cli_resolves_abacus_assets_before_generation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "STRU"
            source.write_text(
                """ATOMIC_SPECIES
Si 28.0855 Si.upf

NUMERICAL_ORBITAL
Si.orb

LATTICE_CONSTANT
1.0

LATTICE_VECTORS
1 0 0
0 1 0
0 0 1

ATOMIC_POSITIONS
Direct

Si
0
1
0 0 0
""",
                encoding="utf-8",
            )
            pp = root / "PSEUDO"
            orb = root / "ORBITAL"
            pp.mkdir()
            orb.mkdir()
            (pp / "Si_PBE.upf").write_text("pseudo", encoding="utf-8")
            (orb / "Si_gga.orb").write_text("orbital", encoding="utf-8")
            previous = Path.cwd()
            try:
                os.chdir(root)
                with patch("zstar.gen_polar.gen_polar") as run_gen:
                    zstar_cli([
                        "gen", "--ensemble", "cartesian", "--stru", "STRU", "--pp", str(pp), "--orb", str(orb)
                    ])
            finally:
                os.chdir(previous)
            kwargs = run_gen.call_args.kwargs
            self.assertEqual(Path(kwargs["f_stru"]).name, "STRU.resolved")
            self.assertIn(str(pp / "Si_PBE.upf"), kwargs["input_sets"])
            self.assertIn(str(orb / "Si_gga.orb"), kwargs["input_sets"])

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
