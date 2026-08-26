import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from zstar.vasp_bec import (
    collect_vasp_bec,
    compare_vasp_bec,
    generate_vasp_backend_script,
    parse_vasp_outcar,
    prepare_vasp_bec,
    render_incar,
)


POSCAR = """BaTiO3
1.0
4.0 0.0 0.0
0.0 4.0 0.0
0.0 0.0 4.0
Ba Ti O
1 1 1
Direct
0.0 0.0 0.0
0.5 0.5 0.5
0.5 0.5 0.0
"""


def outcar_text():
    blocks = []
    matrices = [
        [[2.0, 1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]],
        [[-1.0, 0.0, 0.0], [0.0, -2.0, 0.0], [0.0, 0.0, -3.0]],
        [[-1.0, -1.0, -2.0], [-3.0, -2.0, -5.0], [-6.0, -7.0, -5.0]],
    ]
    for index, matrix in enumerate(matrices, start=1):
        blocks.append(f" ion {index}\n")
        for row, values in enumerate(matrix, start=1):
            blocks.append(f" {row} " + " ".join(str(value) for value in values) + "\n")
    return (
        " MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)\n"
        " ------------------------------------------------------\n"
        " 5.0 0.1 0.2\n 0.1 6.0 0.3\n 0.2 0.3 7.0\n"
        " BORN EFFECTIVE CHARGES (including local field effects) (in |e|)\n"
        " ------------------------------------------------------\n"
        + "".join(blocks)
        + " total drift: 0 0 0\n"
    )


class VaspBecTests(unittest.TestCase):
    def test_render_incar_removes_conflicting_response_tags(self):
        rendered = render_incar(
            "ENCUT = 500\nLEPSILON = .FALSE.\nEDIFF = 1E-10\n",
            updates={"LCALCEPS": ".TRUE."},
            remove=("LEPSILON",),
            defaults={"EDIFF": "1E-8"},
        )
        self.assertIn("ENCUT = 500", rendered)
        self.assertNotIn("LEPSILON", rendered)
        self.assertIn("LCALCEPS = .TRUE.", rendered)
        self.assertEqual(rendered.count("EDIFF"), 1)

    def test_prepare_reference_first_dfpt_workflow(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "input"
            source.mkdir()
            (source / "INCAR").write_text("ENCUT = 500\n", encoding="utf-8")
            (source / "POSCAR").write_text(POSCAR, encoding="utf-8")
            (source / "KPOINTS").write_text("Gamma\n0\nGamma\n2 2 2\n0 0 0\n")
            (source / "POTCAR").write_text("test", encoding="utf-8")
            root = prepare_vasp_bec(source, Path(tmp) / "work")
            reference = (root / "reference" / "INCAR").read_text()
            response = (root / "response" / "INCAR").read_text()
            self.assertIn("LCHARG = .TRUE.", reference)
            self.assertIn("LEPSILON = .TRUE.", response)
            self.assertIn("ICHARG = 1", response)
            manifest = json.loads((root / "vasp_bec_manifest.json").read_text())
            self.assertEqual(manifest["labels"], ["Ba", "Ti", "O"])

    def test_finite_field_replaces_tetrahedron_occupations(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "input"
            source.mkdir()
            (source / "INCAR").write_text(
                "ENCUT = 500\nISMEAR = -5\nSIGMA = 0.2\n", encoding="utf-8"
            )
            (source / "POSCAR").write_text(POSCAR, encoding="utf-8")
            (source / "KPOINTS").write_text("Gamma\n0\nGamma\n2 2 2\n0 0 0\n")
            (source / "POTCAR").write_text("test", encoding="utf-8")
            root = prepare_vasp_bec(
                source, Path(tmp) / "work", method="finite-field"
            )
            for stage in ("reference", "response"):
                rendered = (root / stage / "INCAR").read_text()
                self.assertIn("ISMEAR = 0", rendered)
                self.assertIn("SIGMA = 0.05", rendered)
                self.assertNotIn("ISMEAR = -5", rendered)
            manifest = json.loads((root / "vasp_bec_manifest.json").read_text())
            self.assertIn("ISMEAR=-5", manifest["occupation_override"])
            self.assertIn("EDIFF=unset", manifest["convergence_override"])
            self.assertIn("EDIFF = 1E-8", (root / "response" / "INCAR").read_text())

    def test_outcar_parser_transposes_vasp_to_zstar_convention(self):
        with tempfile.TemporaryDirectory() as tmp:
            outcar = Path(tmp) / "OUTCAR"
            outcar.write_text(outcar_text(), encoding="utf-8")
            epsilon, tensors = parse_vasp_outcar(outcar)
            self.assertTrue(np.allclose(np.diag(epsilon), [5.0, 6.0, 7.0]))
            self.assertTrue(
                np.allclose(
                    tensors[0],
                    [[2.0, 3.0, 6.0], [1.0, 4.0, 7.0], [2.0, 5.0, 8.0]],
                )
            )

    def test_parser_accepts_lepsilon_cumulative_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            outcar = Path(tmp) / "OUTCAR"
            text = outcar_text().replace(
                "BORN EFFECTIVE CHARGES (including local field effects) (in |e|)",
                "BORN EFFECTIVE CHARGES (in e, cumulative output)",
            )
            outcar.write_text(text, encoding="utf-8")
            _epsilon, tensors = parse_vasp_outcar(outcar)
            self.assertEqual(tensors.shape, (3, 3, 3))

    def test_parser_does_not_replace_epsilon_infinity_with_ionic_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            outcar = Path(tmp) / "OUTCAR"
            text = outcar_text() + (
                " MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION\n"
                " ------------------------------------------------------\n"
                " 20.0 0.0 0.0\n 0.0 30.0 0.0\n 0.0 0.0 40.0\n"
            )
            outcar.write_text(text, encoding="utf-8")
            epsilon, _tensors = parse_vasp_outcar(outcar)
            self.assertTrue(np.allclose(np.diag(epsilon), [5.0, 6.0, 7.0]))

    def test_collect_writes_zstar_json_and_phonopy_born(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "response").mkdir()
            (root / "response" / "OUTCAR").write_text(outcar_text(), encoding="utf-8")
            (root / "vasp_bec_manifest.json").write_text(
                json.dumps({"method": "dfpt", "labels": ["Ba", "Ti", "O"]}),
                encoding="utf-8",
            )
            result = collect_vasp_bec(root)
            self.assertTrue((root / "Z-BORN-all.out").is_file())
            self.assertTrue((root / "BORN").is_file())
            self.assertTrue((root / "vasp_bec.json").is_file())
            self.assertTrue((root / "zstar_response.json").is_file())
            response = json.loads((root / "zstar_response.json").read_text())
            self.assertEqual(response["schema"], "zstar-response")
            self.assertEqual(response["dimensionality"]["value"], 3)
            self.assertEqual(response["quantities"][0]["shape"], [3, 3, 3])
            self.assertAlmostEqual(
                np.max(np.abs(np.asarray(result["acoustic_sum_tensor"]))), 0.0
            )

    def test_slurm_script_is_one_serial_driver(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "vasp_bec_manifest.json").write_text(
                json.dumps({"stages": []}), encoding="utf-8"
            )
            script = generate_vasp_backend_script(
                root, backend="slurm", tasks=8, cpus_per_task=2
            )
            text = script.read_text()
            self.assertIn("#SBATCH --ntasks=8", text)
            self.assertIn("srun --ntasks=8 vasp_std", text)
            self.assertEqual(text.count("zstar vasp-bec run"), 1)

    def test_compare_normalized_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = {
                "epsilon_infinity": [[5, 0, 0], [0, 5, 0], [0, 0, 5]],
                "atoms": [{"label": "Si", "tensor": [[2, 0, 0], [0, 2, 0], [0, 0, 2]]}],
            }
            second = {
                "epsilon_infinity": [[5.1, 0, 0], [0, 5.1, 0], [0, 0, 5.1]],
                "atoms": [{"label": "Si", "tensor": [[2.2, 0, 0], [0, 2.2, 0], [0, 0, 2.2]]}],
            }
            first_path = root / "first.json"
            second_path = root / "second.json"
            first_path.write_text(json.dumps(first), encoding="utf-8")
            second_path.write_text(json.dumps(second), encoding="utf-8")
            result = compare_vasp_bec(first_path, second_path)
            self.assertAlmostEqual(result["bec_max_abs_e"], 0.2)
            self.assertAlmostEqual(result["epsilon_max_abs"], 0.1)


if __name__ == "__main__":
    unittest.main()
