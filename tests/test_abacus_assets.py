import tempfile
from pathlib import Path
import unittest

from zstar.abacus_assets import AbacusAssetError, prepare_stru_assets


STRU = """ATOMIC_SPECIES
Si 28.0855 Si.upf
O 15.999 O.upf

NUMERICAL_ORBITAL
Si.orb
O.orb

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

O
0
1
0.5 0.5 0.5
"""


class AbacusAssetTests(unittest.TestCase):
    def test_explicit_directories_prepare_a_copy_without_touching_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "case" / "STRU"
            source.parent.mkdir()
            source.write_text(STRU, encoding="utf-8")
            pp = root / "PSEUDO"
            orb = root / "ORBITAL"
            pp.mkdir()
            orb.mkdir()
            (pp / "Si_ONCV_PBE.upf").write_text("si", encoding="utf-8")
            (pp / "O_ONCV_PBE.upf").write_text("o", encoding="utf-8")
            (orb / "Si_gga.orb").write_text("si orb", encoding="utf-8")
            (orb / "O_gga.orb").write_text("o orb", encoding="utf-8")

            result = prepare_stru_assets(
                source, pp_dir=pp, orb_dir=orb, output_dir=root / ".zstar"
            )

            self.assertTrue(result.changed)
            self.assertEqual(result.path, (root / ".zstar" / "STRU.resolved").resolve())
            self.assertIn("Si_ONCV_PBE.upf", result.path.read_text(encoding="utf-8"))
            self.assertIn("Si_gga.orb", result.path.read_text(encoding="utf-8"))
            self.assertEqual(source.read_text(encoding="utf-8"), STRU)
            self.assertEqual(len(result.assets), 4)
            self.assertTrue((root / ".zstar" / "assets.json").is_file())

    def test_exact_filename_in_directory_disambiguates_versions(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "STRU"
            source.write_text(STRU, encoding="utf-8")
            pp = root / "PSEUDO"
            orb = root / "ORBITAL"
            pp.mkdir()
            orb.mkdir()
            (pp / "Si.upf").write_text("si", encoding="utf-8")
            (pp / "O.upf").write_text("o", encoding="utf-8")
            (orb / "Si.orb").write_text("si", encoding="utf-8")
            (orb / "O.orb").write_text("o", encoding="utf-8")
            result = prepare_stru_assets(source, pp_dir=pp, orb_dir=orb)
            self.assertFalse(result.changed)
            self.assertEqual(len(result.assets), 4)
            self.assertTrue((root / ".zstar" / "assets.json").is_file())

    def test_ambiguous_prefix_match_has_actionable_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "STRU"
            source.write_text(STRU, encoding="utf-8")
            pp = root / "PSEUDO"
            pp.mkdir()
            orb = root / "ORBITAL"
            orb.mkdir()
            (pp / "Si_PBE.upf").write_text("1", encoding="utf-8")
            (pp / "Si_PBEsol.upf").write_text("2", encoding="utf-8")
            (pp / "O.upf").write_text("o", encoding="utf-8")
            (orb / "Si.orb").write_text("si", encoding="utf-8")
            (orb / "O.orb").write_text("o", encoding="utf-8")
            with self.assertRaisesRegex(AbacusAssetError, "Multiple pseudopotential.*Si"):
                prepare_stru_assets(source, pp_dir=pp, orb_dir=orb)

    def test_no_directories_preserves_existing_route(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "STRU"
            source.write_text(STRU, encoding="utf-8")
            result = prepare_stru_assets(source)
            self.assertEqual(result.path, source.resolve())
            self.assertFalse(result.changed)


if __name__ == "__main__":
    unittest.main()
