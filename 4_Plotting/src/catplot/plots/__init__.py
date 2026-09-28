"""Plot contracts and the global built-in registry."""

from .base import OptionSpec, PlotResult, PlotSpec, PlotType
from .builtins import load_builtin_plots
from .registry import PlotRegistry, REGISTRY, register_plot

__all__ = [
    "OptionSpec",
    "PlotRegistry",
    "PlotResult",
    "PlotSpec",
    "PlotType",
    "REGISTRY",
    "register_plot",
    "load_builtin_plots",
]
