import os
from pathlib import Path
import tempfile
import unittest

import numpy as np

from zstar.deal_polar import _write_born_for_phonopy
from zstar.verify_born_symmetry import run_symcheck


class BornPrecisionTests(unittest.TestCase):
    def test_phonopy_born_writer_preserves_eight_decimals(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "BORN"
            _write_born_for_phonopy(
                np.eye(3) * 1.234567891,
                [np.eye(3) * 2.345678912],
                out_path=target,
            )

            text = target.read_text(encoding="utf-8")
            self.assertIn("1.23456789", text)
            self.assertIn("2.34567891", text)

    def test_symmetry_writer_does_not_truncate_reduced_tensors(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "STRU").write_text(
                "LATTICE_CONSTANT\n"
                "5.0\n\n"
                "LATTICE_VECTORS\n"
                "1 0 0\n"
                "0 1 0\n"
                "0 0 1\n\n"
                "ATOMIC_POSITIONS\n"
                "Direct\n\n"
                "Na\n"
                "0.0\n"
                "1\n"
                "0.0 0.0 0.0\n\n"
                "Cl\n"
                "0.0\n"
                "1\n"
                "0.5 0.5 0.5\n",
                encoding="utf-8",
            )
            (root / "Z-BORN-reduced.out").write_text(
                "No. Atom xx xy xz yx yy yz zx zy zz\n"
                "* 1 Na 1.123456789 0.012345678 0 0 1.123456789 0 0 0 1.123456789\n"
                "* 2 Cl -0.923456789 -0.002345678 0 0 -0.923456789 0 0 0 -0.923456789\n",
                encoding="utf-8",
            )

            previous = Path.cwd()
            try:
                os.chdir(root)
                run_symcheck(
                    stru="STRU",
                    reduced="Z-BORN-reduced.out",
                    all=None,
                    out="born_symmetry_report.txt",
                    json_path="born_symmetry_report.json",
                    csv_path=None,
                )
            finally:
                os.chdir(previous)

            text = (root / "Z-BORN-symm.out").read_text(encoding="utf-8")
            self.assertIn("1.02345679", text)
            self.assertIn("0.00734568", text)
            self.assertIn("-1.02345679", text)


if __name__ == "__main__":
    unittest.main()
