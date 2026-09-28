import unittest

import pandas as pd

from cat_tools.mode_tables import merge_mode_columns


class ModeJoinTests(unittest.TestCase):
    def test_periodic_qpoints_join_with_complete_report(self):
        modes = pd.DataFrame({"system": ["x"], "qx": [-0.5], "qy": [0.0], "qz": [0.0], "band": [1]})
        extra = pd.DataFrame(
            {
                "system": ["x"],
                "qx": [0.5],
                "qy": [0.0],
                "qz": [0.0],
                "band": [1],
                "direction_score": [0.8],
            }
        )
        joined, report = merge_mode_columns(modes, extra, columns=["direction_score"])
        self.assertEqual(joined.loc[0, "direction_score"], 0.8)
        self.assertEqual(report.coverage, 1.0)

    def test_partial_join_is_explicit(self):
        modes = pd.DataFrame(
            {
                "qx": [0.0, 0.25],
                "qy": [0.0, 0.0],
                "qz": [0.0, 0.0],
                "band": [0, 0],
            }
        )
        extra = pd.DataFrame({"qx": [0.0], "qy": [0.0], "qz": [0.0], "band": [0], "gamma": [0.1]})
        with self.assertRaisesRegex(ValueError, "matched 1/2"):
            merge_mode_columns(modes, extra, columns=["gamma"])
        joined, report = merge_mode_columns(
            modes,
            extra,
            columns=["gamma"],
            require_complete=False,
        )
        self.assertEqual(report.matched_rows, 1)
        self.assertTrue(pd.isna(joined.loc[1, "gamma"]))

    def test_path_indices_do_not_hide_qpoint_mismatch(self):
        common = {"segment": [0], "point_index": [0], "band": [0], "qy": [0.0], "qz": [0.0]}
        modes = pd.DataFrame({**common, "qx": [0.0]})
        extra = pd.DataFrame({**common, "qx": [0.25], "value": [1.0]})
        with self.assertRaisesRegex(ValueError, "matched 0/1"):
            merge_mode_columns(modes, extra, columns=["value"])


if __name__ == "__main__":
    unittest.main()
