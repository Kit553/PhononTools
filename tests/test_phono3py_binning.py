import unittest

import numpy as np
import pandas as pd

from cat_tools.artifacts import TableArtifact
from cat_tools.phono3py.binning import (
    accumulation_artifact,
    density_artifact,
    profile_artifact,
    spectrum_artifact,
)


class Phono3pyBinningTests(unittest.TestCase):
    def setUp(self):
        self.modes = TableArtifact(
            "mesh_modes",
            pd.DataFrame(
                {
                    "system": "test",
                    "dataset": "mesh",
                    "temperature_K": 300,
                    "ir_index": [0, 1, 2, 3],
                    "band": [0, 0, 0, 0],
                    "frequency_THz": [0.25, 0.75, 1.25, 1.75],
                    "grid_weight": [1, 3, 1, 1],
                    "gamma_THz": [1, 3, 2, 4],
                    "contribution": [1, -0.5, 2, -0.5],
                    "mean_free_path_A": [1, 10, 100, 1000],
                }
            ),
            ("ir_index", "band"),
            {"frequency_THz": "THz", "gamma_THz": "THz", "contribution": "W/(m K)"},
        )

    def test_profile_uses_grid_weights(self):
        profile = profile_artifact(self.modes, name="gamma", value="gamma_THz", bins=[0, 1, 2])
        profile.validate()
        self.assertAlmostEqual(profile.data["value_mean"].iloc[0], 2.5)
        self.assertEqual(profile.data["mode_count"].tolist(), [2, 2])

    def test_density_is_complete_and_normalized(self):
        density = density_artifact(
            self.modes,
            name="density",
            value="gamma_THz",
            frequency_bins=[0, 1, 2],
            value_bins=[0, 2, 4, 6],
        )
        density.validate()
        self.assertEqual(len(density.data), 6)
        self.assertAlmostEqual(density.data["weight_fraction"].sum(), 1.0)

    def test_additive_spectrum_preserves_sign(self):
        spectrum = spectrum_artifact(
            self.modes,
            name="spectrum",
            value="contribution",
            bins=[0, 1, 2],
        )
        spectrum.validate()
        np.testing.assert_allclose(spectrum.data["value_sum"], [0.5, 1.5])
        np.testing.assert_allclose(spectrum.data["value_cumulative"], [0.5, 2.0])

    def test_direct_accumulation_supports_log_mfp(self):
        accumulation = accumulation_artifact(
            self.modes,
            name="kappa_accum_mfp",
            value="contribution",
            coordinate="mean_free_path_A",
            bins={"start": 1, "stop": 1000, "count": 3, "scale": "log"},
            weight=None,
        )
        accumulation.validate()
        np.testing.assert_allclose(accumulation.data["value_cumulative"], [1, 0.5, 2])
        self.assertEqual(accumulation.metadata["coordinate"], "mean_free_path_A")


if __name__ == "__main__":
    unittest.main()
