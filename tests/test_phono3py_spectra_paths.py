import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import yaml

from cat_tools.artifacts import TableArtifact
from cat_tools.phono3py.lineshapes import lorentzian_path_artifact
from cat_tools.phono3py.paths import grid_path, path_artifact
from cat_tools.phono3py.spectra import ir_grid_weights, jdos_artifact, jdos_mesh_artifact


class Phono3pySpectraAndPathTests(unittest.TestCase):
    def test_jdos_channels_and_weighted_mesh(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            files = []
            for grid_point, factor in ((0, 1), (1, 3)):
                path = root / f"jdos-m211-g{grid_point}-s0p1-t300.dat"
                np.savetxt(path, [[0, factor, 2 * factor], [1, 2 * factor, 4 * factor]])
                files.append(path)
            address_file = root / "grid_address-m211.hdf5"
            with h5py.File(address_file, "w") as handle:
                handle["mesh"] = [2, 1, 1]
                handle["grid_address"] = [[0, 0, 0], [1, 0, 0]]
            weights_file = root / "ir_grid_points.yaml"
            weights_file.write_text(
                yaml.safe_dump(
                    {
                        "ir_grid_points": [
                            {"grid_point": 0, "weight": 1},
                            {"grid_point": 1, "weight": 1},
                        ]
                    }
                ),
                encoding="utf-8",
            )
            raw = jdos_artifact(files, system="test", grid_addresses=address_file)
            raw.validate()
            self.assertEqual(set(raw.data["channel"]), {"difference", "sum"})
            averaged = jdos_mesh_artifact(raw, ir_grid_weights(weights_file))
            averaged.validate()
            difference = averaged.data[averaged.data["channel"] == "difference"]
            np.testing.assert_allclose(difference["jdos_per_THz"], [2, 4])
            self.assertTrue(averaged.metadata["coverage_complete"])

    def test_exact_grid_path_and_mode_selection(self):
        path = grid_path([[[0, 0, 0], [0.5, 0, 0]]], [4, 4, 4], labels=["G", "X"])
        self.assertEqual(path["grid_address_a"].tolist(), [0, 1, 2])
        modes = TableArtifact(
            "modes",
            pd.DataFrame(
                {
                    "grid_address_a": np.repeat([0, 1, 2], 2),
                    "grid_address_b": 0,
                    "grid_address_c": 0,
                    "qx": np.repeat([0, 0.25, 0.5], 2),
                    "qy": 0,
                    "qz": 0,
                    "band": np.tile([0, 1], 3),
                    "frequency_THz": np.arange(6),
                }
            ),
            ("grid_address_a", "band"),
            {"frequency_THz": "THz"},
        )
        selected = path_artifact(modes, [[[0, 0, 0], [0.5, 0, 0]]], [4, 4, 4])
        selected.validate()
        self.assertEqual(len(selected.data), 6)
        self.assertEqual(selected.data["trace"].nunique(), 2)

    def test_noncommensurate_path_fails_explicitly(self):
        with self.assertRaisesRegex(ValueError, "not commensurate"):
            grid_path([[[0, 0, 0], [0.5, 0, 0]]], [21, 21, 21])

    def test_lorentzian_path_warns_about_overlapping_modes(self):
        source = TableArtifact(
            "grid_path_modes",
            pd.DataFrame(
                {
                    "system": "test",
                    "dataset": "mesh",
                    "temperature_K": 300,
                    "segment": 0,
                    "point_index": 0,
                    "distance": 0.0,
                    "band": [0, 1],
                    "frequency_THz": [1.0, 1.1],
                    "gamma_THz": [0.1, 0.1],
                }
            ),
            ("dataset", "temperature_K", "segment", "point_index", "band"),
            {"temperature_K": "K", "frequency_THz": "THz", "gamma_THz": "THz"},
        )
        with self.assertWarnsRegex(RuntimeWarning, "overlapping adjacent mode pairs"):
            broadened = lorentzian_path_artifact(
                source,
                frequency_points=5,
                frequency_range_THz=[0, 2],
            )
        broadened.validate()
        self.assertEqual(len(broadened.data), 5)
        self.assertEqual(broadened.data["overlapping_pair_count"].iloc[0], 1)
        self.assertGreater(broadened.data["spectral_weight_per_THz"].max(), 3)


if __name__ == "__main__":
    unittest.main()
