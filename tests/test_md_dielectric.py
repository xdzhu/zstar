import unittest
import sys
import tempfile
from pathlib import Path

import numpy as np

from zstar.md_dielectric import (
    ExternalCommandBECProvider,
    MDDielectricResult,
    MDTrajectory,
    combine_dielectric_contributions,
    dipoles_from_bec,
    load_bec_provider,
    write_outputs,
)


class MDDielectricTests(unittest.TestCase):
    def test_bec_dipole_contraction_uses_displacement_rows(self):
        displacements = np.asarray([[[1.0, 2.0, 3.0]]])
        canonical = np.asarray(
            [[[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]]]
        )

        dipoles = dipoles_from_bec(displacements, canonical)

        np.testing.assert_allclose(dipoles, [[30.0, 36.0, 42.0]])

    @staticmethod
    def trajectory() -> MDTrajectory:
        return MDTrajectory(
            steps=np.asarray([0, 1]),
            positions=np.zeros((2, 1, 3)),
            cells=np.tile(np.eye(3) * 10.0, (2, 1, 1)),
            volumes=np.asarray([1000.0, 1000.0]),
            atom_ids=np.asarray([1]),
            elements=np.asarray(["X"]),
        )

    def test_total_combines_electronic_and_ionic_contributions(self):
        epsilon_ionic = np.diag([4.0, 3.0, 2.0])
        epsilon_electronic = np.diag([5.0, 4.0, 3.0])

        total, electronic, chi_ionic = combine_dielectric_contributions(
            epsilon_ionic, epsilon_electronic
        )

        np.testing.assert_allclose(chi_ionic, np.diag([3.0, 2.0, 1.0]))
        np.testing.assert_allclose(electronic, epsilon_electronic)
        np.testing.assert_allclose(total, np.diag([8.0, 6.0, 4.0]))

    def test_identity_is_used_when_epsilon_infinity_is_omitted(self):
        epsilon_ionic = np.diag([4.0, 3.0, 2.0])
        total, electronic, _ = combine_dielectric_contributions(epsilon_ionic)

        np.testing.assert_allclose(electronic, np.eye(3))
        np.testing.assert_allclose(total, epsilon_ionic)

    def test_external_command_provider_uses_batch_npz_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            script = Path(tmp) / "predict.py"
            script.write_text(
                "import os, numpy as np\n"
                "request=np.load(os.environ['ZSTAR_MD_REQUEST'])\n"
                "shape=(len(request['steps']),len(request['atom_ids']),3,3)\n"
                "out=np.tile(np.eye(3),(shape[0],shape[1],1,1))\n"
                "np.save(os.environ['ZSTAR_MD_OUTPUT'],out)\n",
                encoding="utf-8",
            )
            provider = ExternalCommandBECProvider(
                f'"{sys.executable}" "{script}"'
            )
            values = provider.provide(self.trajectory())
            self.assertEqual(values.shape, (2, 1, 3, 3))
            np.testing.assert_allclose(
                values[:, 0], np.tile(np.eye(3), (2, 1, 1))
            )

    def test_module_callable_provider_is_supported(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "demo_provider.py").write_text(
                "import numpy as np\n"
                "def provider(trajectory):\n"
                "    return np.zeros((len(trajectory.steps),len(trajectory.atom_ids),3,3))\n",
                encoding="utf-8",
            )
            sys.path.insert(0, str(root))
            try:
                provider = load_bec_provider("demo_provider:provider")
                self.assertEqual(provider.provide(self.trajectory()).shape, (2, 1, 3, 3))
            finally:
                sys.path.remove(str(root))
                sys.modules.pop("demo_provider", None)

    def test_outputs_include_common_response_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = MDDielectricResult(
                epsilon=np.eye(3) * 5.0,
                epsilon_ionic=np.eye(3) * 2.0,
                ionic_susceptibility=np.eye(3),
                electronic_dielectric=np.eye(3) * 4.0,
                electronic_source="test",
                covariance_eA2=np.eye(3),
                dipoles_eA=np.zeros((2, 3)),
                steps=np.asarray([0, 1]),
                selected=np.asarray([True, True]),
                volume_A3_avg=1000.0,
                temperature_K=300.0,
                reference_mode="mean",
                bec_mode="external-command",
            )
            root = Path(tmp)
            write_outputs(root, result, [None, None])
            self.assertTrue((root / "zstar_response.json").is_file())


if __name__ == "__main__":
    unittest.main()
