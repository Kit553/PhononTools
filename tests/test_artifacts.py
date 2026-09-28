import tempfile
import unittest
from pathlib import Path

import pandas as pd
import yaml

from cat_tools.artifacts import TableArtifact, write_artifacts


class ArtifactTests(unittest.TestCase):
    def test_writes_one_format_and_manifest(self):
        artifact = TableArtifact(
            "values",
            pd.DataFrame({"id": [0, 1], "value": [2.0, 3.0]}),
            ("id",),
            {"value": "arb."},
        )
        with tempfile.TemporaryDirectory() as temporary:
            written = write_artifacts([artifact], temporary)
            self.assertEqual({path.name for path in written}, {"values.csv", "manifest.yaml"})
            manifest = yaml.safe_load((Path(temporary) / "manifest.yaml").read_text())
            self.assertEqual(manifest["artifacts"][0]["keys"], ["id"])

    def test_rejects_duplicate_primary_key(self):
        artifact = TableArtifact("bad", pd.DataFrame({"id": [0, 0]}), ("id",))
        with self.assertRaises(ValueError):
            artifact.validate()

    def test_rejects_duplicate_artifact_names_before_writing(self):
        artifact = TableArtifact("values", pd.DataFrame({"id": [0]}), ("id",))
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(ValueError, "must be unique"):
                write_artifacts([artifact, artifact], temporary)
            self.assertFalse((Path(temporary) / "values.csv").exists())

    def test_rejects_unsafe_names_and_unknown_unit_columns(self):
        with self.assertRaises(ValueError):
            TableArtifact("../values", pd.DataFrame({"id": [0]}), ("id",)).validate()
        with self.assertRaisesRegex(ValueError, "units reference absent"):
            TableArtifact(
                "values", pd.DataFrame({"id": [0]}), ("id",), {"missing": "eV"}
            ).validate()

    def test_rerun_removes_only_files_owned_by_previous_manifest(self):
        first = TableArtifact("first", pd.DataFrame({"id": [0]}), ("id",))
        second = TableArtifact("second", pd.DataFrame({"id": [0]}), ("id",))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_artifacts([first, second], root)
            unrelated = root / "notes.txt"
            unrelated.write_text("keep me", encoding="utf-8")

            write_artifacts([first], root, table_format="parquet")

            self.assertFalse((root / "first.csv").exists())
            self.assertFalse((root / "second.csv").exists())
            self.assertTrue((root / "first.parquet").is_file())
            self.assertEqual(unrelated.read_text(encoding="utf-8"), "keep me")
            manifest = yaml.safe_load((root / "manifest.yaml").read_text())
            self.assertRegex(manifest["artifacts"][0]["sha256"], r"^[0-9a-f]{64}$")

    def test_invalid_later_artifact_does_not_touch_previous_run(self):
        good = TableArtifact("values", pd.DataFrame({"id": [0]}), ("id",))
        bad = TableArtifact("bad", pd.DataFrame({"id": [0, 0]}), ("id",))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_artifacts([good], root)
            before = (root / "values.csv").read_bytes()
            with self.assertRaises(ValueError):
                write_artifacts([good, bad], root)
            self.assertEqual((root / "values.csv").read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
