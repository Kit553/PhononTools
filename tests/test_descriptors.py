import unittest

import numpy as np
import pandas as pd

from cat_tools.artifacts import TableArtifact
from cat_tools.phonopy.descriptors import acoustic_velocity_summary, averaged_descriptors


class DescriptorTests(unittest.TestCase):
    def test_averaged_descriptors_cover_full_dos_moment_set(self):
        dos = TableArtifact(
            "dos",
            pd.DataFrame(
                {
                    "system": "test",
                    "frequency_THz": [0.0, 1.0, 2.0],
                    "total_DOS": [1.0, 1.0, 1.0],
                    "Na_DOS": [0.5, 0.5, 0.5],
                }
            ),
            ("frequency_THz",),
        )
        result = averaged_descriptors(
            dos,
            min_frequency_THz=0.0,
            low_frequency_cutoff_THz=1.0,
            inverse_square_floor_THz=0.5,
        )
        row = result.data.iloc[0]
        self.assertEqual(result.name, "averaged_descriptors")
        self.assertEqual(result.keys, ("system", "source", "curve"))
        self.assertEqual(set(result.data["curve"]), {"total", "Na"})
        self.assertEqual(row["source"], "dos")
        self.assertAlmostEqual(row["integrated_weight"], 2.0)
        self.assertAlmostEqual(row["mean_frequency_THz"], 1.0)
        self.assertAlmostEqual(row["median_frequency_THz"], 1.0)
        self.assertAlmostEqual(row["std_frequency_THz"], np.sqrt(0.5))
        self.assertAlmostEqual(row["rms_frequency_THz"], np.sqrt(1.5))
        self.assertAlmostEqual(row["log_mean_frequency_THz"], 1.0)
        self.assertAlmostEqual(row["low_frequency_weight"], 1.0)
        self.assertAlmostEqual(row["low_frequency_fraction"], 0.5)
        self.assertAlmostEqual(row["low_frequency_average_THz"], 0.5)
        self.assertAlmostEqual(row["inverse_square_weight_per_THz2"], 3.125)
        self.assertAlmostEqual(row["inverse_square_moment_per_THz2"], 1.5625)
        self.assertAlmostEqual(row["inverse_square_effective_frequency_THz"], 0.8)
        self.assertAlmostEqual(row["inverse_square_average_frequency_THz"], 0.4)
        self.assertAlmostEqual(row["debye_equivalent_mean_THz"], 4 / 3)
        self.assertAlmostEqual(
            row["debye_equivalent_mean_K"],
            row["debye_equivalent_mean_THz"] * result.metadata["THz_to_K"],
        )
        result.validate()

    def test_exact_mesh_row_is_qpoint_weight_normalized(self):
        dos = TableArtifact(
            "dos",
            pd.DataFrame(
                {
                    "system": "test",
                    "frequency_THz": [1.0, 2.0],
                    "total_DOS": [1.0, 1.0],
                }
            ),
            ("frequency_THz",),
        )
        mesh = TableArtifact(
            "mesh_modes",
            pd.DataFrame(
                {
                    "system": "test",
                    "q_index": [0, 0, 1, 1],
                    "band": [0, 1, 0, 1],
                    "frequency_THz": [1.0, 2.0, 3.0, 4.0],
                    "weight": [1.0, 1.0, 3.0, 3.0],
                }
            ),
            ("q_index", "band"),
        )
        result = averaged_descriptors(
            dos,
            mesh,
            min_frequency_THz=0.0,
            low_frequency_cutoff_THz=2.0,
        )
        row = result.data.query("source == 'mesh_modes'").iloc[0]
        self.assertEqual(row["curve"], "all_modes")
        self.assertAlmostEqual(row["integrated_weight"], 2.0)
        self.assertAlmostEqual(row["mean_frequency_THz"], 3.0)
        self.assertEqual(result.metadata["mesh_qpoint_weight_normalization"], 4.0)
        result.validate()

    def test_descriptor_settings_are_checked(self):
        dos = TableArtifact(
            "dos",
            pd.DataFrame(
                {
                    "system": "test",
                    "frequency_THz": [0.0, 1.0],
                    "total_DOS": [1.0, 1.0],
                }
            ),
            ("frequency_THz",),
        )
        with self.assertRaisesRegex(ValueError, "inverse_square_floor"):
            averaged_descriptors(dos, inverse_square_floor_THz=0.0)

    def test_near_gamma_velocity_summary(self):
        rows = []
        for q_index, qx in enumerate((0.05, 0.10)):
            for band, speed in enumerate((10.0, 20.0, 30.0)):
                rows.append(
                    {
                        "system": "test",
                        "q_index": q_index,
                        "band": band,
                        "qx": qx,
                        "qy": 0.0,
                        "qz": 0.0,
                        "frequency_THz": 1.0,
                        "group_velocity": speed,
                    }
                )
        modes = TableArtifact("mesh_modes", pd.DataFrame(rows), ("q_index", "band"))
        result = acoustic_velocity_summary(modes, native_to_m_per_s=100.0)
        values = result.data.set_index("mode")["mean_speed_m_s"]
        self.assertEqual(values["transverse"], 1500.0)
        self.assertEqual(values["longitudinal"], 3000.0)


if __name__ == "__main__":
    unittest.main()
