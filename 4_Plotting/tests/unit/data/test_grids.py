import unittest

import numpy as np
import pandas as pd

from catplot.data import interpolate_grid, rectilinear_grid
from catplot.errors import DataError


class GridTests(unittest.TestCase):
    def test_rectilinear_grid_has_y_x_orientation(self) -> None:
        frame = pd.DataFrame(
            {
                "x": [0, 1, 0, 1],
                "y": [0, 0, 1, 1],
                "z": [10, 20, 30, 40],
            }
        )
        grid = rectilinear_grid(frame, x="x", y="y", z="z")
        self.assertEqual(grid.z.shape, (2, 2))
        self.assertTrue(np.array_equal(grid.z, [[10, 20], [30, 40]]))
        self.assertEqual(grid.extent, (0.0, 1.0, 0.0, 1.0))

    def test_incomplete_grid_requires_explicit_interpolation(self) -> None:
        frame = pd.DataFrame({"x": [0, 1, 0], "y": [0, 0, 1], "z": [1, 2, 3]})
        with self.assertRaisesRegex(DataError, "interpolation"):
            rectilinear_grid(frame, x="x", y="y", z="z")

    def test_scattered_interpolation_shape_and_flag(self) -> None:
        frame = pd.DataFrame(
            {
                "x": [0, 1, 0, 1, 0.5],
                "y": [0, 0, 1, 1, 0.5],
                "z": [0, 1, 1, 2, 1],
            }
        )
        grid = interpolate_grid(
            frame,
            x="x",
            y="y",
            z="z",
            nx=11,
            ny=7,
            fill_nearest=True,
        )
        self.assertEqual(grid.z.shape, (7, 11))
        self.assertTrue(grid.interpolated)
        self.assertTrue(np.isfinite(grid.z).all())


if __name__ == "__main__":
    unittest.main()

