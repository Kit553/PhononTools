import unittest

import numpy as np
import pandas as pd

from catplot.data import apply_transforms
from catplot.errors import DataError


class TransformTests(unittest.TestCase):
    def setUp(self) -> None:
        self.frame = pd.DataFrame(
            {
                "x": np.linspace(0, 4, 5),
                "y": np.array([0.0, 1.0, 4.0, 9.0, 16.0]),
            }
        )

    def test_input_is_not_modified(self) -> None:
        transformed = apply_transforms(
            self.frame,
            [{"op": "scale", "column": "y", "output": "scaled", "factor": 2}],
        )
        self.assertNotIn("scaled", self.frame)
        self.assertEqual(transformed["scaled"].tolist(), [0, 2, 8, 18, 32])

    def test_area_normalization_uses_x(self) -> None:
        transformed = apply_transforms(
            self.frame,
            [{"op": "normalize", "column": "y", "output": "density", "mode": "area", "x": "x"}],
        )
        self.assertAlmostEqual(np.trapezoid(transformed["density"], transformed["x"]), 1.0)

    def test_derivative_and_sort_can_be_chained(self) -> None:
        reversed_frame = self.frame.iloc[::-1].reset_index(drop=True)
        transformed = apply_transforms(
            reversed_frame,
            [
                {"op": "sort", "by": "x"},
                {"op": "derivative", "column": "y", "output": "dy_dx", "x": "x"},
            ],
        )
        self.assertTrue(np.allclose(transformed["dy_dx"], 2 * transformed["x"]))

    def test_unknown_transform_is_rejected(self) -> None:
        with self.assertRaises(DataError):
            apply_transforms(self.frame, [{"op": "magic", "column": "y"}])


if __name__ == "__main__":
    unittest.main()

