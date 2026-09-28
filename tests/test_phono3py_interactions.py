import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

from cat_tools.phono3py.interactions import interaction_artifact


class Phono3pyInteractionTests(unittest.TestCase):
    def test_explicit_band_selection_preserves_triplets(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "pp-m222-g17.hdf5"
            with h5py.File(path, "w") as handle:
                handle["pp"] = np.arange(16).reshape(2, 2, 2, 2)
                handle["triplet"] = [[17, 3, 4], [17, 5, 6]]
                handle["weight"] = [2, 4]
            artifact = interaction_artifact([path], system="test", bands={"i": [1]}, min_pp_eV2=8)
            artifact.validate()
            self.assertEqual(len(artifact.data), 4)
            self.assertEqual(set(artifact.data["band_i"]), {1})
            self.assertEqual(set(artifact.data["grid_point"]), {17})


if __name__ == "__main__":
    unittest.main()
