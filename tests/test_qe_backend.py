import json
from pathlib import Path
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO

import numpy as np

from zstar.cli import zstar_cli
from zstar.qe_backend import (
    collect_qe_response,
    generate_qe_backend_script,
    parse_qe_dynmat_output,
    parse_qe_gap,
    parse_qe_ph_output,
    prepare_qe_response,
    qe_namelist_values,
    render_qe_namelist,
)


PW_INPUT = """&CONTROL
  calculation = 'relax',
  prefix = 'SiC',
  pseudo_dir = './pseudo_src',
  outdir = './tmp',
/
&SYSTEM
  ibrav = 1,
  celldm(1) = 8.2,
  nat = 2,
  ntyp = 2,
  ecutwfc = 40,
  nbnd = 8,
/
&ELECTRONS
  conv_thr = 1.0d-10,
/
ATOMIC_SPECIES
Si 28.085 Si.UPF
C 12.011 C.UPF
ATOMIC_POSITIONS crystal
Si 0 0 0
C 0.25 0.25 0.25
K_POINTS automatic
2 2 2 0 0 0
"""


PH_OUTPUT = """
     Dielectric constant in cartesian axis

     (  6.500000 0.100000 0.000000 )
     (  0.100000 6.600000 0.000000 )
     (  0.000000 0.000000 6.700000 )

     Polarizability (a.u.)^3                    Polarizability (A^3)
     20.0 0.0 0.0       2.963694 0.0 0.0
     0.0 21.0 0.0       0.0 3.111879 0.0
     0.0 0.0 22.0       0.0 0.0 3.260063

     Effective charges (d Force / dE) in cartesian axis

      atom 1 Si
 Ex ( 2.0 0.1 0.2 )
 Ey ( 0.3 2.1 0.4 )
 Ez ( 0.5 0.6 2.2 )
      atom 2 C
 Ex ( -2.0 -0.1 -0.2 )
 Ey ( -0.3 -2.1 -0.4 )
 Ez ( -0.5 -0.6 -2.2 )

     Diagonalizing the dynamical matrix
     freq ( 1) = 1.000000 [THz] = 33.356410 [cm-1]
     freq ( 2) = 10.000000 [THz] = 333.564100 [cm-1]
     JOB DONE.
"""


DYNMAT_OUTPUT = """
 # mode   [cm-1]    [THz]      IR        Raman   depol.fact
 1 33.3564 1.0000 0.0000 0.0000 0.7500
 2 333.5641 10.0000 2.5000 12.0000 0.1000

 JOB DONE.
"""


class QeBackendTests(unittest.TestCase):
    def test_namelist_update_preserves_cards(self):
        updated = render_qe_namelist(PW_INPUT, "control", {"calculation": "'scf'", "outdir": "'../scratch'"})
        values = qe_namelist_values(updated, "control")
        self.assertEqual(values["calculation"], "'scf'")
        self.assertEqual(values["outdir"], "'../scratch'")
        self.assertIn("ATOMIC_POSITIONS crystal", updated)

    def test_prepare_disables_gamma_only_restart_format(self):
        gamma_input = PW_INPUT.replace(
            "K_POINTS automatic\n2 2 2 0 0 0", "K_POINTS gamma"
        )
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source"
            (source / "pseudo_src").mkdir(parents=True)
            (source / "pw.in").write_text(gamma_input, encoding="utf-8")
            for name in ("Si.UPF", "C.UPF"):
                (source / "pseudo_src" / name).write_text("pseudo", encoding="utf-8")
            root = prepare_qe_response(source / "pw.in", Path(tmp) / "work")
            prepared = (root / "scf" / "pw.in").read_text()
            self.assertIn("K_POINTS automatic\n1 1 1 0 0 0", prepared)

    def test_prepare_copies_pseudos_and_builds_three_stages(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source"
            (source / "pseudo_src").mkdir(parents=True)
            (source / "pw.in").write_text(PW_INPUT, encoding="utf-8")
            for name in ("Si.UPF", "C.UPF"):
                (source / "pseudo_src" / name).write_text("pseudo", encoding="utf-8")
            root = prepare_qe_response(source / "pw.in", Path(tmp) / "work", raman=True)
            manifest = json.loads((root / "qe_response_manifest.json").read_text())
            self.assertEqual([stage["name"] for stage in manifest["stages"]], ["scf", "phonon", "dynmat"])
            self.assertTrue((root / "pseudo" / "Si.UPF").is_file())
            self.assertIn("recover = .true.", (root / "phonon" / "ph.in").read_text())

    def test_gap_and_response_parsers(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pw = root / "pw.out"
            pw.write_text("highest occupied, lowest unoccupied level (ev): 4.2 5.1\n", encoding="utf-8")
            self.assertAlmostEqual(parse_qe_gap(pw), 0.9)
            ph_path = root / "ph.out"
            ph_path.write_text(PH_OUTPUT, encoding="utf-8")
            parsed = parse_qe_ph_output(ph_path)
            self.assertEqual(parsed["born_tensors"].shape, (2, 3, 3))
            self.assertTrue(np.allclose(parsed["born_tensors"][0][0], [2.0, 0.3, 0.5]))
            dm_path = root / "dynmat.out"
            dm_path.write_text(DYNMAT_OUTPUT, encoding="utf-8")
            modes = parse_qe_dynmat_output(dm_path)
            self.assertEqual(modes["raman_activities"], [0.0, 12.0])

    def test_collect_writes_common_response_and_spectra(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "phonon").mkdir()
            (root / "dynmat").mkdir()
            (root / "phonon" / "ph.out").write_text(PH_OUTPUT, encoding="utf-8")
            (root / "dynmat" / "dynmat.out").write_text(DYNMAT_OUTPUT, encoding="utf-8")
            (root / "qe_response_manifest.json").write_text(
                json.dumps({
                    "prefix": "SiC", "raman_requested": True,
                    "dimensionality": {"value": 3, "periodic_axes": ["x", "y", "z"]},
                }),
                encoding="utf-8",
            )
            result = collect_qe_response(root, plot=False, points=101)
            self.assertTrue((root / "zstar_response.json").is_file())
            response = json.loads((root / "zstar_response.json").read_text())
            names = [quantity["name"] for quantity in response["quantities"]]
            self.assertIn("born_effective_charge", names)
            self.assertIn("raman_activity", names)
            self.assertEqual(result["raman_activities"], [0.0, 12.0])

    def test_molecular_qe_response_names_atomic_polar_tensor(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "phonon").mkdir()
            (root / "dynmat").mkdir()
            (root / "phonon" / "ph.out").write_text(PH_OUTPUT, encoding="utf-8")
            (root / "dynmat" / "dynmat.out").write_text(
                DYNMAT_OUTPUT, encoding="utf-8"
            )
            (root / "qe_response_manifest.json").write_text(
                json.dumps(
                    {
                        "prefix": "molecule",
                        "raman_requested": True,
                        "dimensionality": {"value": 0, "periodic_axes": []},
                    }
                ),
                encoding="utf-8",
            )
            collect_qe_response(root, plot=False, points=101)
            response = json.loads((root / "zstar_response.json").read_text())
            names = [quantity["name"] for quantity in response["quantities"]]
            self.assertIn("atomic_polar_tensor", names)
            self.assertNotIn("born_effective_charge", names)

    def test_script_supports_slurm(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "qe_response_manifest.json").write_text("{}", encoding="utf-8")
            script = generate_qe_backend_script(root, backend="slurm", tasks=20)
            text = script.read_text()
            self.assertIn("#SBATCH --ntasks=20", text)
            self.assertIn("srun --ntasks=20 pw.x", text)
            self.assertIn("zstar qe-bec run", text)
            self.assertIn("zstar qe-bec collect", text)

    def test_qe_bec_and_legacy_aliases_prepare_the_same_workflow(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source"
            (source / "pseudo_src").mkdir(parents=True)
            input_path = source / "pw.in"
            input_path.write_text(PW_INPUT, encoding="utf-8")
            for name in ("Si.UPF", "C.UPF"):
                (source / "pseudo_src" / name).write_text("pseudo", encoding="utf-8")

            canonical_root = Path(tmp) / "canonical"
            with redirect_stdout(StringIO()):
                zstar_cli([
                    "qe-bec", "prepare", "--input", str(input_path),
                    "--root", str(canonical_root),
                ])
            self.assertTrue((canonical_root / "qe_response_manifest.json").is_file())

            legacy_root = Path(tmp) / "legacy"
            stderr = StringIO()
            with redirect_stdout(StringIO()), redirect_stderr(stderr):
                zstar_cli([
                    "qe", "prepare", "--input", str(input_path),
                    "--root", str(legacy_root),
                ])
            self.assertTrue((legacy_root / "qe_response_manifest.json").is_file())
            self.assertIn("DEPRECATED", stderr.getvalue())

            backend_legacy_root = Path(tmp) / "backend-legacy"
            stderr = StringIO()
            with redirect_stdout(StringIO()), redirect_stderr(stderr):
                zstar_cli([
                    "backend", "qe", "prepare", "--input", str(input_path),
                    "--root", str(backend_legacy_root),
                ])
            self.assertTrue(
                (backend_legacy_root / "qe_response_manifest.json").is_file()
            )
            self.assertIn("DEPRECATED", stderr.getvalue())

    def test_help_exposes_canonical_bec_and_keeps_backend_query_only(self):
        stdout = StringIO()
        with self.assertRaisesRegex(SystemExit, "0"), redirect_stdout(stdout):
            zstar_cli(["--help"])
        help_text = stdout.getvalue()
        self.assertIn("backend", help_text)
        self.assertIn("bec", help_text)
        self.assertNotIn("qe-bec", help_text)
        self.assertNotIn("==SUPPRESS==", help_text)

        stdout = StringIO()
        with self.assertRaisesRegex(SystemExit, "0"), redirect_stdout(stdout):
            zstar_cli(["backend", "--help"])
        backend_help = stdout.getvalue()
        self.assertIn("{list}", backend_help)
        self.assertNotIn("{list,qe}", backend_help)


if __name__ == "__main__":
    unittest.main()
