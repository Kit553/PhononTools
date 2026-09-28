"""Scatter renderer with continuous/categorical colours and categorical markers."""

from __future__ import annotations

import numpy as np
import pandas as pd

from catplot.data.columns import select_column
from catplot.encoding import encode_color, encode_markers, iter_groups
from catplot.plots import OptionSpec, PlotResult, PlotSpec, register_plot

from ._common import color_at, context_theme, numeric

SPEC = PlotSpec(
    name="scatter",
    family="cartesian",
    description="Draw points with optional colour, marker, size, and group encodings.",
    required_mappings=("x", "y"),
    optional_mappings=("color", "marker", "size", "group"),
    aliases=("points",),
    options={
        "label": OptionSpec((str, type(None)), default=None),
        "color": OptionSpec((str, type(None)), default=None),
        "marker": OptionSpec(str, default="o"),
        "size": OptionSpec((int, float), default=18.0),
        "size_min": OptionSpec((int, float), default=10.0),
        "size_max": OptionSpec((int, float), default=60.0),
        "cmap": OptionSpec(str, default="viridis"),
        "norm": OptionSpec(str, default="linear"),
        "vmin": OptionSpec((int, float, type(None)), default=None),
        "vmax": OptionSpec((int, float, type(None)), default=None),
        "center": OptionSpec((int, float, type(None)), default=None),
        "alpha": OptionSpec((int, float), default=1.0),
        "edgecolor": OptionSpec((str, type(None)), default=None),
        "linewidth": OptionSpec((int, float), default=0.5),
        "zorder": OptionSpec((int, float), default=3.0),
        "colorbar": OptionSpec(bool, default=True),
        "colorbar_label": OptionSpec((str, type(None)), default=None),
    },
)


def _sizes(data, mapping, options):
    if "size" not in mapping:
        return np.full(len(data), float(options["size"]))
    values = numeric(data, mapping["size"])
    finite = np.isfinite(values)
    result = np.full(len(values), options["size_min"], dtype=float)
    if finite.any():
        lower, upper = np.nanmin(values), np.nanmax(values)
        scaled = np.zeros_like(values) if lower == upper else (values - lower) / (upper - lower)
        result = options["size_min"] + scaled * (options["size_max"] - options["size_min"])
    return result


@register_plot(SPEC)
def draw(*, ax, data, mapping, options, context):
    theme = context_theme(context)
    color_encoding = None
    if "color" in mapping:
        color_encoding = encode_color(
            select_column(data, mapping["color"]),
            theme,
            cmap=options["cmap"], norm=options["norm"],
            vmin=options["vmin"], vmax=options["vmax"], center=options["center"],
        )
    marker_encoding = encode_markers(select_column(data, mapping["marker"])) if "marker" in mapping else None
    size_values = pd.Series(_sizes(data, mapping, options), index=data.index)
    artists = []
    mappable = None
    for group in iter_groups(data, mapping.get("group")):
        subset = group.data
        marker_groups = [(options["marker"], subset)]
        if marker_encoding is not None:
            marker_series = pd.Series(marker_encoding.values, index=data.index).loc[subset.index]
            marker_groups = [(mark, subset.loc[marker_series == mark]) for mark in pd.unique(marker_series)]
        for marker, marked in marker_groups:
            if marked.empty:
                continue
            kwargs = {
                "s": size_values.loc[marked.index].to_numpy(),
                "marker": marker,
                "alpha": options["alpha"],
                "edgecolors": options["edgecolor"] if options["edgecolor"] is not None else "none",
                "linewidths": options["linewidth"],
                "zorder": options["zorder"],
                "label": options["label"] or group.label or None,
            }
            if color_encoding is None:
                kwargs["color"] = options["color"] or color_at(context, group.index)
            elif color_encoding.categorical:
                colors = pd.Series(color_encoding.values, index=data.index).loc[marked.index]
                kwargs["color"] = colors.to_numpy()
            else:
                values = select_column(marked, mapping["color"], numeric=True).to_numpy(dtype=float)
                kwargs.update(c=values, cmap=color_encoding.cmap, norm=color_encoding.norm)
            collection = ax.scatter(numeric(marked, mapping["x"]), numeric(marked, mapping["y"]), **kwargs)
            artists.append(collection)
            if color_encoding is not None and not color_encoding.categorical:
                mappable = collection
    return PlotResult(artists=tuple(artists), mappable=mappable,
                      metadata={"colorbar": options["colorbar"] and mappable is not None,
                                "colorbar_label": options["colorbar_label"],
                                "color_mapping": getattr(color_encoding, "mapping", {}),
                                "marker_mapping": getattr(marker_encoding, "mapping", {})})
