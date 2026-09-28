import importlib.util
import unittest
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
RECIPES = sorted((ROOT / "3_Analysis" / "recipes").rglob("*.plot.yaml"))


class RecipeTests(unittest.TestCase):
    def test_every_recipe_is_yaml_and_scatter_uses_canonical_size_option(self):
        self.assertTrue(RECIPES)
        for path in RECIPES:
            with self.subTest(path=path.relative_to(ROOT)):
                document = yaml.safe_load(path.read_text(encoding="utf-8"))
                self.assertEqual(document["schema"], 1)
                for panel in document["panels"]:
                    for layer in panel["layers"]:
                        if layer["type"] == "scatter":
                            self.assertNotIn("s", layer.get("options", {}))

    @unittest.skipUnless(importlib.util.find_spec("catplot"), "CATPlot is not installed")
    def test_every_recipe_validates_with_catplot(self):
        from catplot.recipe import load_recipe

        for path in RECIPES:
            with self.subTest(path=path.relative_to(ROOT)):
                load_recipe(path)


if __name__ == "__main__":
    unittest.main()
