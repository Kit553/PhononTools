"""Coordinate data, layout, axes, plot modules, legends, and colorbars."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from catplot.axes import add_colorbar, configure_axes, configure_legend, get_axis_target
from catplot.data import apply_filters, apply_transforms, load_table
from catplot.layout import create_layout
from catplot.plots import REGISTRY, load_builtin_plots
from catplot.recipe.validate import validate_recipe
from catplot.style import apply_style

from .artists import LayerRecord
from .context import RenderContext


@dataclass
class RenderedFigure:
    figure: object
    axes: dict[str, object]
    sources: dict[str, pd.DataFrame]
    layers: list[LayerRecord]
    recipe: object


def prepare_sources(recipe, root: Path):
    prepared = {}
    for source_id, source in recipe.data.items():
        path = Path(source.path)
        if not path.is_absolute(): path = root / path
        frame = load_table(path, file_format=source.format, **source.load)
        frame = apply_filters(frame, source.filters)
        frame = apply_transforms(frame, source.transforms)
        prepared[source_id] = frame
    return prepared


def render_recipe(recipe, *, recipe_root="."):
    validate_recipe(recipe); load_builtin_plots()
    theme = apply_style(recipe.style.get("profile", "paper"), recipe.style.get("theme", "default"))
    root = Path(recipe_root).resolve(); sources = prepare_sources(recipe, root)
    projections = {panel.id: panel.projection for panel in recipe.panels}
    layout = create_layout([panel.id for panel in recipe.panels], projections, recipe.layout)
    records = []; twin_cache = {}
    for panel in recipe.panels:
        ax = layout.axes[panel.id]
        configure_axes(ax, panel.axes, projection=panel.projection)
        panel_twins = []
        for layer_index, layer in enumerate(panel.layers):
            if not layer.visible: continue
            target = get_axis_target(ax, layer.axis, twin_cache)
            if target is not ax and target not in panel_twins: panel_twins.append(target)
            if target is not ax and panel.axes.get("right_y"):
                right = panel.axes["right_y"]
                if right.get("label") is not None: target.set_ylabel(right["label"])
                if right.get("limits") is not None: target.set_ylim(*right["limits"])
                if right.get("scale") is not None: target.set_yscale(right["scale"])
            plot_type = REGISTRY.get(layer.type)
            data = sources[layer.data] if layer.data is not None else pd.DataFrame()
            context = RenderContext(layout.figure, target, panel.id, layer_index, theme, recipe, sources)
            result = plot_type.draw(ax=target, data=data, mapping=layer.mapping, options=layer.options, context=context)
            records.append(LayerRecord(panel.id, layer_index, plot_type.spec.name, layer.axis, result))
            if result.mappable is not None and result.metadata.get("colorbar", False):
                add_colorbar(layout.figure, target, result.mappable,
                             panel.colorbar if panel.colorbar is not None else True,
                             label=result.metadata.get("colorbar_label"))
        configure_legend(ax, panel.legend, other_axes=panel_twins)
        if panel.axes.get("panel_label"):
            label_text = getattr(ax, "text2D", ax.text)
            label_text(-0.12, 1.03, str(panel.axes["panel_label"]), transform=ax.transAxes,
                       ha="left", va="bottom", fontweight="bold", clip_on=False)
    return RenderedFigure(layout.figure, layout.axes, sources, records, recipe)
