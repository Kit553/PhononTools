import unittest

import numpy as np

from cat_tools.phonopy.projections import Projection, element_projections, mode_character


class ProjectionTests(unittest.TestCase):
    def setUp(self):
        self.eigenvectors = np.eye(6, dtype=complex)[None, :, :]
        self.masses = [23.0, 35.0]

    def test_element_fractions_and_participation_ratio(self):
        ratio, projections = mode_character(
            self.eigenvectors,
            self.masses,
            element_projections(["Na", "Cl"]),
        )
        np.testing.assert_allclose(ratio, 0.5)
        np.testing.assert_allclose(projections["Na"][0, :3], 1.0)
        np.testing.assert_allclose(projections["Na"][0, 3:], 0.0)
        np.testing.assert_allclose(projections["Na"] + projections["Cl"], 1.0)

    def test_cartesian_projection(self):
        _, projections = mode_character(
            self.eigenvectors,
            self.masses,
            [Projection("x", direction=(1.0, 0.0, 0.0))],
        )
        np.testing.assert_allclose(projections["x"], [[1, 0, 0, 1, 0, 0]])

    def test_displacement_weighting(self):
        mixed = np.zeros((1, 6, 6), dtype=complex)
        mixed[0, 0, 0] = mixed[0, 3, 0] = 1 / np.sqrt(2)
        _, projections = mode_character(
            mixed,
            self.masses,
            element_projections(["Na", "Cl"]),
            weighting="displacement",
        )
        expected = (1 / 23.0) / (1 / 23.0 + 1 / 35.0)
        self.assertAlmostEqual(projections["Na"][0, 0], expected)


if __name__ == "__main__":
    unittest.main()
