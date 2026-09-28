from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from catplot.data import infer_format, load_many, load_table
from catplot.errors import DataError


class LoaderTests(unittest.TestCase):
    def test_csv_comments_and_header_whitespace(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "dos.csv"
            path.write_text("# generated\n frequency , NaDOS\n0,1\n1,2\n", encoding="utf-8")
            frame = load_table(path)
            self.assertEqual(list(frame.columns), ["frequency", "NaDOS"])
            self.assertEqual(frame.attrs["source_format"], "csv")
            self.assertEqual(frame["NaDOS"].tolist(), [1, 2])

    def test_whitespace_data(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "map.dat"
            path.write_text("x y z\n0 0 1\n1 0 2\n", encoding="utf-8")
            frame = load_table(path)
            self.assertEqual(frame.shape, (2, 3))

    def test_column_subset_is_checked(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "data.csv"
            path.write_text("frequency,NaDOS,total\n0,1,2\n", encoding="utf-8")
            frame = load_table(path, columns=["frequency", "NaDOS"])
            self.assertEqual(list(frame.columns), ["frequency", "NaDOS"])

    def test_many_sources_keep_explicit_ids(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "a.csv"
            second = root / "b.tsv"
            first.write_text("x,y\n0,1\n", encoding="utf-8")
            second.write_text("x\ty\n0\t2\n", encoding="utf-8")
            loaded = load_many({"first": first, "second": second})
            self.assertEqual(list(loaded), ["first", "second"])
            self.assertEqual(loaded["second"].loc[0, "y"], 2)

    def test_unknown_extension_requires_explicit_format(self) -> None:
        with self.assertRaises(DataError):
            infer_format("values.unknown")


if __name__ == "__main__":
    unittest.main()

