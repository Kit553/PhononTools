"""Cross-reference validation beyond individual dataclass parsing."""

from __future__ import annotations

from catplot.errors import RecipeError
from catplot.plots import REGISTRY, load_builtin_plots


def validate_recipe(recipe):
    load_builtin_plots()
    if not recipe.panels: raise RecipeError("A recipe needs at least one panel")
    ids = [panel.id for panel in recipe.panels]
    if len(ids) != len(set(ids)): raise RecipeError("Panel IDs must be unique")
    if recipe.output.name in {"", ".", ".."}: raise RecipeError("Output name is invalid")
    for fmt in recipe.output.formats:
        if fmt.lower() not in {"pdf", "png", "svg"}: raise RecipeError(f"Unsupported output format {fmt!r}")
    if "pdf" not in {fmt.lower() for fmt in recipe.output.formats}:
        raise RecipeError("Every recipe must export PDF")
    for panel in recipe.panels:
        if panel.projection not in {"2d", "3d"}: raise RecipeError(f"Panel {panel.id!r} projection must be 2d or 3d")
        if not panel.layers: raise RecipeError(f"Panel {panel.id!r} has no layers")
        for layer in panel.layers:
            try: plot_type = REGISTRY.get(layer.type)
            except Exception as exc: raise RecipeError(str(exc)) from exc
            if panel.projection not in plot_type.spec.projections:
                raise RecipeError(f"Plot {layer.type!r} does not support projection {panel.projection!r}")
            if layer.data is not None and layer.data not in recipe.data:
                raise RecipeError(f"Layer references unknown data source {layer.data!r}")
            if layer.data is None and layer.mapping:
                raise RecipeError(f"Mapped layer {layer.type!r} needs a data source")
            if layer.axis not in {"left", "right"}: raise RecipeError("Layer axis must be left or right")
            plot_type.spec.validate_mapping(layer.mapping)
            plot_type.spec.resolve_options(layer.options)
    return recipe
