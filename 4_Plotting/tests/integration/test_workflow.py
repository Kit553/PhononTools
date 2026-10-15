from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import matplotlib.pyplot as plt
from pypdf import PdfReader

from catplot.cli.main import main
from catplot.recipe import load_recipe
from catplot.render import render_to_files


RECIPE = """\
schema: 1
style: {profile: paper, theme: default}
layout: {type: single, edge_mm: 85}
data:
  values: {path: data.csv}
panels:
  - id: main
    axes:
      x: {label: x}
      y: {label: y}
    legend: false
    layers:
      - {type: line, data: values, mapping: {x: x, y: y}}
output: {name: check, directory: output, formats: [pdf], preview_dpi: 100, final_dpi: 600}
"""


class WorkflowTests(unittest.TestCase):
    def test_draft_is_lightweight_and_final_is_portable(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "data.csv").write_text("x,y\n0,0\n1,1\n2,4\n", encoding="utf-8")
            recipe_path = root / "plot.yaml"
            recipe_path.write_text(RECIPE, encoding="utf-8")
            recipe = load_recipe(recipe_path)

            draft = render_to_files(recipe, recipe_path)
            plt.close(draft.rendered.figure)
            self.assertTrue((root / "output/draft/check.pdf").is_file())
            self.assertFalse(any((root / "output/draft").glob("*.lock.json")))

            final = render_to_files(recipe, recipe_path, final=True, command="catplot render plot.yaml --final")
            plt.close(final.rendered.figure)
            folder = root / "output/final/check"
            self.assertTrue((folder / "check.pdf").is_file())
            self.assertTrue((folder / "check.plot.yaml").is_file())
            self.assertTrue((folder / "check.plot.lock.json").is_file())

            page = PdfReader(folder / "check.pdf").pages[0]
            self.assertLess(abs(float(page.mediabox.width) - float(page.mediabox.height)), 0.01)

            # Render the copied YAML from its new directory to check relative paths.
            portable = load_recipe(folder / "check.plot.yaml")
            rerender = render_to_files(portable, folder / "check.plot.yaml")
            plt.close(rerender.rendered.figure)

    def test_quick_cli_creates_recipe_and_draft(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            data = root / "values.csv"
            data.write_text("x,y\n0,1\n1,2\n", encoding="utf-8")
            recipe = root / "recipes/quick.plot.yaml"
            result = main([
                "line", str(data), "--x", "x", "--y", "y",
                "--marker", "o", "--output", "quick", "--recipe", str(recipe),
            ])
            self.assertEqual(result, 0)
            self.assertTrue(recipe.is_file())
            self.assertTrue((root / "draft/quick.pdf").is_file())


if __name__ == "__main__":
    unittest.main()
