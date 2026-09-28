"""Small public interface for CATPlot.

Most implementation modules stay importable for advanced use, but ordinary
scripts should begin with the names re-exported here.
"""

from .plots import (
    OptionSpec,
    PlotRegistry,
    PlotResult,
    PlotSpec,
    PlotType,
    REGISTRY,
    register_plot,
)
from .render import new_panel, save_panel
from .style import (
    DEFAULT_EDGE_MM,
    FINAL_DPI,
    apply_style,
    available_profiles,
    available_themes,
    load_theme,
    mm_to_inches,
)

__all__ = [
    "DEFAULT_EDGE_MM",
    "FINAL_DPI",
    "apply_style",
    "available_profiles",
    "available_themes",
    "load_theme",
    "mm_to_inches",
    "new_panel",
    "save_panel",
    "OptionSpec",
    "PlotRegistry",
    "PlotResult",
    "PlotSpec",
    "PlotType",
    "REGISTRY",
    "register_plot",
]

__version__ = "0.9.0"
