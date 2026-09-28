"""Colour, marker, grouping, normalization, and label encodings."""

from .colors import ColorEncoding, encode_color
from .grouping import Group, iter_groups
from .labels import column_label, series_label
from .markers import DEFAULT_MARKERS, MarkerEncoding, encode_markers
from .normalization import make_norm

__all__ = [
    "ColorEncoding",
    "DEFAULT_MARKERS",
    "Group",
    "MarkerEncoding",
    "column_label",
    "encode_color",
    "encode_markers",
    "iter_groups",
    "make_norm",
    "series_label",
]

