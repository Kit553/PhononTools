"""High-level draft/final render workflow."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import os
from pathlib import Path

from catplot.provenance import build_lock, write_lock
from catplot.recipe import dump_recipe

from .export import save_figure
from .pipeline import render_recipe


@dataclass
class RenderOutputs:
    rendered: object
    files: list[Path]
    directory: Path


def render_to_files(recipe, recipe_path, *, final=False, command=None):
    recipe_path = Path(recipe_path).resolve(); root = recipe_path.parent
    rendered = render_recipe(recipe, recipe_root=root)
    base = (root / recipe.output.directory).resolve()
    target = base / "final" / recipe.output.name if final else base / "draft"
    target.mkdir(parents=True, exist_ok=True)
    dpi = recipe.output.final_dpi if final else recipe.output.preview_dpi
    metadata = {str(k): str(v) for k, v in recipe.metadata.items() if v is not None}
    files = []
    for file_format in recipe.output.formats:
        output = target / f"{recipe.output.name}.{file_format.lower()}"
        files.append(save_figure(rendered.figure, output, dpi=dpi, metadata=metadata))
    if final:
        # Rebase relative data paths for the copied recipe in the final output directory.
        portable_recipe = deepcopy(recipe)
        for source in portable_recipe.data.values():
            source_path = Path(source.path)
            if not source_path.is_absolute():
                source.path = os.path.relpath((root / source_path).resolve(), target)
        portable_recipe.output.directory = "."
        final_recipe = dump_recipe(portable_recipe, target / f"{recipe.output.name}.plot.yaml")
        files.append(final_recipe)
        lock = build_lock(recipe, recipe_path, command=command)
        files.append(write_lock(lock, target / f"{recipe.output.name}.plot.lock.json"))
    return RenderOutputs(rendered, files, target)
