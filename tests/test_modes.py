import unittest

import numpy as np

from cat_tools.phonopy.modes import ModeData, mode_artifact


class ModeTableTests(unittest.TestCase):
    def test_mesh_schema(self):
        data = ModeData(
            source="mesh",
            qpoints=np.array([[0, 0, 0], [0.5, 0, 0]]),
            frequencies=np.array([[0, 1], [2, 3]], dtype=float),
            group_velocities=np.ones((2, 2, 3)),
            weights=np.array([1, 2]),
        )
        artifact = mode_artifact(
            data,
            system="test",
            symbols=["Na"],
            masses=[23],
            include_participation_ratio=False,
        )
        self.assertEqual(artifact.name, "mesh_modes")
        self.assertEqual(len(artifact.data), 4)
        self.assertEqual(list(artifact.data["band"]), [0, 1, 0, 1])
        np.testing.assert_allclose(artifact.data["group_velocity"], np.sqrt(3))
        artifact.validate()

    def test_band_trace_separates_segments(self):
        data = ModeData(
            source="band",
            qpoints=np.zeros((4, 3)),
            frequencies=np.ones((4, 2)),
            distances=np.arange(4),
            segments=np.array([0, 0, 1, 1]),
            point_indices=np.array([0, 1, 0, 1]),
        )
        artifact = mode_artifact(
            data,
            system="test",
            symbols=["Na"],
            masses=[23],
            include_participation_ratio=False,
        )
        self.assertEqual(sorted(artifact.data["trace"].unique()), [0, 1, 2, 3])
        artifact.validate()


if __name__ == "__main__":
    unittest.main()
