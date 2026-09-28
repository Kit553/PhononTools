import unittest

import numpy as np
import pandas as pd

from cat_tools.phono3py.cloud import Crystal, first_brillouin_zone, prepare_cloud


class Phono3pyCloudTests(unittest.TestCase):
    def setUp(self):
        self.crystal = Crystal(np.eye(3), np.zeros((1, 3)), np.ones(1, dtype=int))

    def test_cubic_first_bz(self):
        zone = first_brillouin_zone(self.crystal.reciprocal_lattice)
        self.assertEqual(len(zone.vertices), 8)
        self.assertEqual(len(zone.faces), 6)
        self.assertTrue(np.allclose(abs(zone.vertices).max(axis=0), np.pi))

    def test_irreducible_point_expands_to_full_star(self):
        modes = pd.DataFrame(
            {
                "dataset": "mesh",
                "temperature_K": 300,
                "ir_index": [0, 0],
                "band": [0, 1],
                "qx": [0.25, 0.25],
                "qy": [0, 0],
                "qz": [0, 0],
                "grid_weight": [6, 6],
                "frequency_THz": [1, 2],
                "gamma_THz": [1, 3],
            }
        )
        cloud = prepare_cloud(modes, self.crystal, value="gamma_THz", reduction="max")
        self.assertEqual(len(cloud), 6)
        self.assertTrue(np.allclose(cloud["value"], 3))
        self.assertTrue(np.allclose(cloud["normalized"], 1))


if __name__ == "__main__":
    unittest.main()
