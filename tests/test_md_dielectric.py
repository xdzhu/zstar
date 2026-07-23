import unittest

import numpy as np

from zstar.md_dielectric import combine_dielectric_contributions


class MDDielectricTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
