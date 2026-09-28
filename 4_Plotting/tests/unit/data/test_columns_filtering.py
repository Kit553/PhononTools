import unittest

import numpy as np
import pandas as pd

from catplot.data import apply_filters, resolve_mapping, select_column
from catplot.errors import DataError


class ColumnTests(unittest.TestCase):
    def setUp(self) -> None:
        self.frame = pd.DataFrame(
            {
                "temperature": [250, 300, 350, 400],
                "conductivity": [0.1, 0.3, 0.7, 1.2],
                "compound": ["PSe4", "PSe4", "SbSe4", "SbSe4"],
            }
        )

    def test_missing_column_suggests_close_name(self) -> None:
        with self.assertRaisesRegex(DataError, "temperature"):
            select_column(self.frame, "temprature")

    def test_mapping_can_resolve_multiple_y_columns(self) -> None:
        frame = self.frame.assign(conductivity_2=self.frame["conductivity"] * 2)
        mapping = resolve_mapping(
            frame,
            {"x": "temperature", "y": ["conductivity", "conductivity_2"]},
            numeric_roles={"x", "y"},
        )
        self.assertEqual(len(mapping["y"]), 2)

    def test_finite_check_rejects_nan(self) -> None:
        frame = self.frame.copy()
        frame.loc[1, "conductivity"] = np.nan
        with self.assertRaisesRegex(DataError, "non-finite"):
            select_column(frame, "conductivity", finite=True)


class FilterTests(unittest.TestCase):
    def setUp(self) -> None:
        self.frame = pd.DataFrame(
            {
                "temperature": [250, 300, 350, 400],
                "compound": ["PSe4", "PSe4", "SbSe4", "SbSe4"],
            }
        )

    def test_structured_filters_combine_with_and(self) -> None:
        filtered = apply_filters(
            self.frame,
            [
                {"column": "temperature", "op": "between", "value": [300, 400]},
                {"column": "compound", "op": "eq", "value": "SbSe4"},
            ],
        )
        self.assertEqual(filtered["temperature"].tolist(), [350, 400])

    def test_contains_is_literal_not_regex(self) -> None:
        frame = pd.DataFrame({"label": ["a.b", "acb"]})
        filtered = apply_filters(
            frame,
            [{"column": "label", "op": "contains", "value": "."}],
        )
        self.assertEqual(filtered["label"].tolist(), ["a.b"])


if __name__ == "__main__":
    unittest.main()

