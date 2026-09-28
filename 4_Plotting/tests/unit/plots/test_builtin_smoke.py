"""Render every built-in at least once to catch API drift in Matplotlib."""

from types import SimpleNamespace
import unittest

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from catplot.plots import REGISTRY, load_builtin_plots
from catplot.style import load_theme


class BuiltinSmokeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        load_builtin_plots()
        rows = []
        for y in (-1.0, 0.0, 1.0):
            for x in (-1.0, 0.0, 1.0):
                rows.append({"x": x, "y": y, "z": np.exp(-(x*x+y*y)), "u": -y,
                             "v": x, "color": x+y+2, "size": x*x+y*y+1})
        cls.grid = pd.DataFrame(rows)
        cls.series = pd.DataFrame({
            "x": np.linspace(0, 4, 9),
            "y": np.array([0, .4, .8, 1.0, .7, .3, .1, .5, .9]),
            "lower": np.array([0, .4, .8, 1.0, .7, .3, .1, .5, .9]) - .1,
            "upper": np.array([0, .4, .8, 1.0, .7, .3, .1, .5, .9]) + .1,
            "err": .08,
            "color": np.linspace(-1, 1, 9),
            "group": ["a"] * 5 + ["b"] * 4,
            "marker": ["first", "second", "first"] * 3,
            "size": np.arange(1, 10),
            "text": list("ABCDEFGHI"),
            "dx": .1,
            "dy": .1,
        })
        cls.categories = pd.DataFrame({
            "x": ["A"] * 4 + ["B"] * 4 + ["C"] * 4,
            "y": [1.0, 1.2, .9, 1.1, 1.4, 1.6, 1.5, 1.45, .7, .9, .8, .85],
        })

    def _draw(self, name, data, mapping, options=None, projection="2d"):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d" if projection == "3d" else None)
        context = SimpleNamespace(theme=load_theme("default"))
        result = REGISTRY.get(name).draw(ax=ax, data=data, mapping=mapping,
                                         options=options or {}, context=context)
        self.assertTrue(result.artists or result.mappable is not None)
        fig.canvas.draw()
        plt.close(fig)

    def test_cartesian_distribution_and_annotation_renderers(self):
        cases = [
            ("line", self.series, {"x": "x", "y": "y"}, {}),
            ("scatter", self.series, {"x": "x", "y": "y", "color": "color", "marker": "marker", "size": "size"}, {}),
            ("errorbar", self.series, {"x": "x", "y": "y", "yerr": "err"}, {}),
            ("band", self.series, {"x": "x", "lower": "lower", "upper": "upper"}, {}),
            ("area", self.series, {"x": "x", "y": "y"}, {}),
            ("step", self.series, {"x": "x", "y": "y"}, {}),
            ("stem", self.series, {"x": "x", "y": "y"}, {}),
            ("colored-line", self.series, {"x": "x", "y": "y", "color": "color"}, {}),
            ("histogram", self.series, {"x": "y"}, {}),
            ("density", self.series, {"x": "y"}, {}),
            ("ecdf", self.series, {"x": "y"}, {}),
            ("reference", pd.DataFrame(), {}, {"value": .5}),
            ("span", pd.DataFrame(), {}, {"minimum": .2, "maximum": .4}),
            ("text", self.series.iloc[:2], {"x": "x", "y": "y", "text": "text"}, {}),
            ("arrow", self.series.iloc[:2], {"x": "x", "y": "y", "dx": "dx", "dy": "dy"}, {}),
        ]
        for case in cases:
            with self.subTest(plot=case[0]):
                self._draw(*case)

    def test_categorical_renderers(self):
        for name in ("bar", "box", "violin", "strip", "swarm"):
            with self.subTest(plot=name):
                self._draw(name, self.categories, {"x": "x", "y": "y"})

    def test_field_renderers(self):
        cases = [
            ("heatmap", {"x": "x", "y": "y", "z": "z"}),
            ("contour", {"x": "x", "y": "y", "z": "z"}),
            ("quiver", {"x": "x", "y": "y", "u": "u", "v": "v", "color": "color"}),
            ("stream", {"x": "x", "y": "y", "u": "u", "v": "v", "color": "color"}),
        ]
        for name, mapping in cases:
            with self.subTest(plot=name):
                self._draw(name, self.grid, mapping)

    def test_three_dimensional_renderers(self):
        for name in ("line3d", "scatter3d", "surface", "wireframe"):
            mapping = {"x": "x", "y": "y", "z": "z"}
            data = self.series.assign(z=self.series["y"] * .5) if name == "line3d" else self.grid
            with self.subTest(plot=name):
                self._draw(name, data, mapping, projection="3d")


if __name__ == "__main__":
    unittest.main()
