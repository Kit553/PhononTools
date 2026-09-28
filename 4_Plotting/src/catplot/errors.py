"""CATPlot's public exception hierarchy.

Callers can catch :class:`CATPlotError` for any expected user-facing failure,
or a narrower subclass when they can recover from a particular problem.
"""


class CATPlotError(Exception):
    """Base class for expected CATPlot failures."""


class ConfigurationError(CATPlotError):
    """A global configuration or style setting is invalid."""


class RecipeError(CATPlotError):
    """A plot recipe is malformed or internally inconsistent."""


class DataError(CATPlotError):
    """Input data cannot satisfy the requested mapping or transformation."""


class PlotSpecError(CATPlotError):
    """A plot type received invalid mappings, options, or renderer output."""


class UnknownPlotTypeError(PlotSpecError):
    """The requested plot type is not registered."""


class DuplicatePlotTypeError(PlotSpecError):
    """A plot name or alias is already owned by another renderer."""


class RenderError(CATPlotError):
    """A validated plot could not be rendered or exported."""

