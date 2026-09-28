"""Explicit built-in import list; importing a module registers its renderer."""

from __future__ import annotations

from importlib import import_module

_MODULES = (
    "catplot.plots.cartesian.line", "catplot.plots.cartesian.colored_line",
    "catplot.plots.cartesian.scatter", "catplot.plots.cartesian.errorbar",
    "catplot.plots.cartesian.band", "catplot.plots.cartesian.area",
    "catplot.plots.cartesian.step", "catplot.plots.cartesian.stem",
    "catplot.plots.categorical.bar", "catplot.plots.categorical.box",
    "catplot.plots.categorical.violin", "catplot.plots.categorical.strip",
    "catplot.plots.categorical.swarm", "catplot.plots.distributions.histogram",
    "catplot.plots.distributions.density", "catplot.plots.distributions.ecdf",
    "catplot.plots.fields2d.heatmap", "catplot.plots.fields2d.contour",
    "catplot.plots.fields2d.quiver", "catplot.plots.fields2d.stream",
    "catplot.plots.threed.line", "catplot.plots.threed.scatter",
    "catplot.plots.threed.surface", "catplot.plots.threed.wireframe",
    "catplot.plots.annotations.reference", "catplot.plots.annotations.span",
    "catplot.plots.annotations.text", "catplot.plots.annotations.arrow",
)

_loaded = False


def load_builtin_plots() -> None:
    global _loaded
    if _loaded:
        return
    for module in _MODULES:
        import_module(module)
    _loaded = True

