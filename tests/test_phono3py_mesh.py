import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

from cat_tools.phono3py.mesh import load_mesh, mesh_artifact
from cat_tools.phono3py.transport import pair_artifact, transport_artifact


def write_kappa(path: Path, *, pairs: bool = False) -> None:
    with h5py.File(path, "w") as handle:
        handle["mesh"] = [2, 1, 1]
        handle["qpoint"] = [[0, 0, 0], [0.5, 0, 0]]
        handle["frequency"] = [[1, 2], [3, 4]]
        handle["weight"] = [1, 1]
        handle["temperature"] = [100, 300]
        handle["boundary_mfp"] = 10.0
        handle["gamma"] = np.full((2, 2, 2), 0.25)
        handle["gamma_N"] = np.full((2, 2, 2), 0.10)
        handle["gamma_U"] = np.full((2, 2, 2), 0.15)
        handle["gamma_isotope"] = np.full((2, 2), 0.05)
        handle["heat_capacity"] = np.ones((2, 2, 2))
        handle["group_velocity"] = np.ones((2, 2, 3))
        handle["ave_pp"] = np.arange(4).reshape(2, 2)
        handle["mode_kappa"] = np.ones((2, 2, 2, 6))
        handle["kappa"] = np.arange(12).reshape(2, 6)
        if pairs:
            handle["mode_kappa_C"] = np.ones((2, 2, 2, 2, 6))


class Phono3pyMeshTests(unittest.TestCase):
    def test_modal_and_total_tables(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "kappa-m211-bmfp10.hdf5"
            write_kappa(path)
            source = load_mesh(path)
            modes = mesh_artifact([source], [300], system="test")
            modes.validate()
            self.assertEqual(len(modes.data), 4)
            self.assertTrue(np.allclose(modes.data["mode_kappa_iso_W_mK"], 0.5))
            self.assertTrue(np.allclose(modes.data["lifetime_ps"], 1 / np.pi))
            self.assertTrue(np.allclose(modes.data["mean_free_path_A"], np.sqrt(3) / np.pi))
            self.assertTrue(np.allclose(modes.data["gamma_isotope_THz"], 0.05))

            total = transport_artifact([source], system="test")
            total.validate()
            self.assertEqual(set(total.data["channel"]), {"total"})
            self.assertEqual(total.data["boundary_label"].iloc[0], "10 um")

    def test_pair_tensor_keeps_both_band_indices(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "kappa-m211.hdf5"
            write_kappa(path, pairs=True)
            pairs = pair_artifact(load_mesh(path), 300, system="test")
            pairs.validate()
            self.assertEqual(len(pairs.data), 8)
            self.assertIn("band_i", pairs.data)
            self.assertIn("band_j", pairs.data)
            self.assertTrue(np.allclose(pairs.data["mode_kappa_coherence_iso_W_mK"], 0.5))

    def test_generalized_grid_does_not_invent_addresses(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "kappa-m211-g4.hdf5"
            write_kappa(path)
            with h5py.File(path, "a") as handle:
                handle["grid_matrix"] = [[2, 1, 0], [0, 2, 0], [0, 0, 1]]
            source = load_mesh(path)
            modes = mesh_artifact([source], [300], system="test")
            self.assertNotIn("grid_address_a", modes.data)
            self.assertEqual(source.mode_normalization, 4)


if __name__ == "__main__":
    unittest.main()
