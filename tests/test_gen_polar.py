import tempfile
from pathlib import Path
import unittest

from zstar.phonopy_stru import write_phonopy_compatible_stru


class GenPolarTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
