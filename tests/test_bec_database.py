import csv
import json
from pathlib import Path
import tempfile
import unittest

from zstar.bec_database import collect_database, read_born, read_zborn, write_manifest_template


class BecDatabaseTests(unittest.TestCase):
    def test_read_born_and_collect_ranked_bulk(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case = root / "case"
            case.mkdir()
            (case / "BORN").write_text(
                "# dielectric and two atoms\n"
                "5 0 0 0 6 0 0 0 7\n"
                "2 0 0 0 2 0 0 0 2\n"
                "-2 0 0 0 -2 0 0 0 -2\n",
                encoding="utf-8",
            )
            (case / "0.no-move").mkdir()
            (case / "0.no-move" / "zstar_insulation.json").write_text(
                json.dumps({"gap_eV": 2.4, "insulating": True}), encoding="utf-8"
            )
            response = case / "dielectric_response"
            response.mkdir()
            (response / "ir_response_real.dat").write_text(
                "0.0 20 0 0 0 30 0 0 0 40\n1.0 1 0 0 0 1 0 0 0 1\n",
                encoding="utf-8",
            )
            manifest = root / "candidates.csv"
            with manifest.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["material_id", "formula", "dimensionality", "workspace"])
                writer.writerow(["test-1", "AB", 3, "case"])
            summary = collect_database(manifest, root / "database")
            self.assertEqual(summary["complete"], 1)
            self.assertEqual(summary["ranked_high_k_3d"], 1)
            with (root / "database" / "materials.csv").open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(float(rows[0]["k_static_mean"]), 30.0)
            self.assertEqual(float(rows[0]["acoustic_sum_max_abs_e"]), 0.0)
            _epsilon, tensors = read_born(case / "BORN")
            self.assertEqual(tensors.shape, (2, 3, 3))

    def test_template_and_missing_results_are_auditable(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = write_manifest_template(root / "candidates.csv")
            summary = collect_database(manifest, root / "database")
            self.assertEqual(summary["incomplete"], 1)
            row = json.loads((root / "database" / "materials.jsonl").read_text().splitlines()[0])
            self.assertIn("missing_bec_tensor", row["quality_flags"])

    def test_full_zborn_takes_precedence_over_representative_born(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case = root / "case"
            case.mkdir()
            (case / "STRU").write_text(
                "ATOMIC_POSITIONS\nDirect\n\nA\n0\n1\n0 0 0\n\nB\n0\n2\n0 0 0\n0 0 0\n",
                encoding="utf-8",
            )
            (case / "BORN").write_text(
                "5 0 0 0 5 0 0 0 5\n2 0 0 0 2 0 0 0 2\n-1 0 0 0 -1 0 0 0 -1\n",
                encoding="utf-8",
            )
            (case / "Z-BORN-symm.out").write_text(
                "No. Atom xx xy xz yx yy yz zx zy zz\n"
                "* 1 A 2 0 0 0 2 0 0 0 2\n"
                "* 2 B -1 0 0 0 -1 0 0 0 -1\n"
                "  3 B -1 0 0 0 -1 0 0 0 -1\n",
                encoding="utf-8",
            )
            manifest = root / "manifest.csv"
            manifest.write_text(
                "material_id,formula,dimensionality,workspace\nabc,AB2,3,case\n",
                encoding="utf-8",
            )
            summary = collect_database(manifest, root / "database")
            record = json.loads((root / "database" / "materials.jsonl").read_text().splitlines()[0])
            self.assertEqual(summary["atom_tensors"], 3)
            self.assertEqual(record["tensor_scope"], "full_cell")
            self.assertEqual(record["acoustic_sum_max_abs_e"], 0.0)
            self.assertEqual(read_zborn(case / "Z-BORN-symm.out").shape, (3, 3, 3))

    def test_molecular_spectra_are_complete_auxiliary(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case = root / "molecule" / "ir"
            case.mkdir(parents=True)
            (case / "ir_modes.csv").write_text("mode,frequency\n1,1000\n", encoding="utf-8")
            manifest = root / "manifest.csv"
            manifest.write_text(
                "material_id,formula,dimensionality,workspace\nch4,CH4,0,molecule\n",
                encoding="utf-8",
            )
            summary = collect_database(manifest, root / "database")
            record = json.loads((root / "database" / "materials.jsonl").read_text().splitlines()[0])
            self.assertEqual(summary["complete_auxiliary"], 1)
            self.assertEqual(record["status"], "complete_auxiliary")
            self.assertNotIn("missing_bec_tensor", record["quality_flags"])

    def test_one_dimensional_response_is_stored_but_not_ranked(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case = root / "wire"
            (case / "0.no-move").mkdir(parents=True)
            (case / "BORN").write_text(
                "1.1 0 0 0 1.2 0 0 0 1.3\n"
                "2 0 0 0 2 0 0 0 2\n"
                "-2 0 0 0 -2 0 0 0 -2\n",
                encoding="utf-8",
            )
            (case / "0.no-move" / "zstar_insulation.json").write_text(
                json.dumps({"gap_eV": 2.1, "insulating": True}),
                encoding="utf-8",
            )
            (case / "zstar_response.json").write_text(
                json.dumps(
                    {
                        "quantities": [
                            {
                                "name": "line_polarizability",
                                "unit": "angstrom^2",
                                "normalization": "isolated_object",
                                "convention": "gaussian",
                                "values": [[1, 0, 0], [0, 2, 0], [0, 0, 3]],
                            }
                        ]
                    }
                ),
                encoding="utf-8",
            )
            manifest = root / "manifest.csv"
            manifest.write_text(
                "material_id,formula,dimensionality,workspace\n"
                "wire-1,GaAs,1,wire\n",
                encoding="utf-8",
            )

            summary = collect_database(manifest, root / "database")
            record = json.loads(
                (root / "database" / "materials.jsonl").read_text().splitlines()[0]
            )

            self.assertEqual(summary["complete"], 1)
            self.assertEqual(summary["ranked_high_k_3d"], 0)
            self.assertEqual(record["response_kind"], "line_1d")
            self.assertEqual(record["intrinsic_response"]["unit"], "angstrom^2")
            self.assertEqual(record["intrinsic_response_mean_diagonal"], 2.0)
            self.assertEqual(record["high_k_rank_basis"], "not_applicable")
            self.assertNotIn("k_electronic_mean", record)


if __name__ == "__main__":
    unittest.main()
