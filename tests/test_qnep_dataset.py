import csv
import json
from pathlib import Path
import tempfile
import unittest

from zstar.qnep_dataset import (
    augment_qnep_dataset,
    check_qnep_dataset,
    write_qnep_input,
)


FRAME = """2
Lattice="5 0 0 0 5 0 0 0 5" energy=-10 Properties=species:S:1:pos:R:3:force:R:3
Na 0 0 0 0.1 0.2 0.3
Cl 2.5 2.5 2.5 -0.1 -0.2 -0.3
"""


ZBORN = """No. Atom xx xy xz yx yy yz zx zy zz
 1 Na 1 2 3 4 5 6 7 8 9
 2 Cl -1 -2 -3 -4 -5 -6 -7 -8 -9
"""


class QnepDatasetTests(unittest.TestCase):
    def test_augment_transposes_zstar_tensor_for_qnep(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "train.xyz"
            source.write_text(FRAME, encoding="utf-8")
            bec = root / "Z-BORN-all.out"
            bec.write_text(ZBORN, encoding="utf-8")
            output = root / "train_qnep.xyz"
            summary = augment_qnep_dataset(source, output, bec=bec)
            lines = output.read_text().splitlines()
            self.assertIn(":bec:R:9", lines[1])
            self.assertEqual(
                lines[2].split()[-9:],
                ["1", "4", "7", "2", "5", "8", "3", "6", "9"],
            )
            self.assertEqual(summary["labeled_frames"], 1)
            self.assertTrue(Path(summary["audit_output"]).is_file())

    def test_partial_labels_from_frame_map_are_allowed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "train.xyz"
            source.write_text(FRAME + FRAME, encoding="utf-8")
            bec = root / "Z-BORN-all.out"
            bec.write_text(ZBORN, encoding="utf-8")
            mapping = root / "map.csv"
            with mapping.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["frame", "bec"])
                writer.writerow([1, bec.name])
            output = root / "qnep.xyz"
            summary = augment_qnep_dataset(source, output, bec_map=mapping)
            self.assertEqual(summary["labeled_frames"], 1)
            self.assertEqual(summary["unlabeled_frames"], 1)
            checked = check_qnep_dataset(output)
            self.assertEqual(checked["labeled_frames"], 1)

    def test_atom_order_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "train.xyz"
            source.write_text(FRAME, encoding="utf-8")
            data = {
                "tensor_convention": "rows=atomic displacement/force; columns=polarization/electric field",
                "atoms": [
                    {"label": "Cl", "tensor": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]},
                    {"label": "Na", "tensor": [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]},
                ],
            }
            bec = root / "bec.json"
            bec.write_text(json.dumps(data), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Atom order mismatch"):
                augment_qnep_dataset(source, root / "out.xyz", bec=bec)

    def test_minimal_qnep_input_uses_dataset_elements(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "train.xyz"
            source.write_text(FRAME, encoding="utf-8")
            bec = root / "Z-BORN-all.out"
            bec.write_text(ZBORN, encoding="utf-8")
            labeled = root / "qnep.xyz"
            augment_qnep_dataset(source, labeled, bec=bec)
            nep = write_qnep_input(labeled, root / "nep.in", charge_mode=2)
            text = nep.read_text()
            self.assertIn("type 2 Cl Na", text)
            self.assertIn("charge_mode 2", text)
            self.assertIn("lambda_z 0.5", text)


if __name__ == "__main__":
    unittest.main()
