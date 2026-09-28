"""YAML recipe loading."""

from __future__ import annotations

from pathlib import Path
import yaml

from catplot.errors import RecipeError
from .model import Recipe
from .validate import validate_recipe


def load_recipe(path, *, validate=True):
    target = Path(path)
    if not target.is_file(): raise RecipeError(f"Recipe does not exist: {target}")
    try:
        with target.open("r", encoding="utf-8") as stream: raw = yaml.safe_load(stream)
    except (OSError, yaml.YAMLError) as exc:
        raise RecipeError(f"Could not read recipe {target}: {exc}") from exc
    recipe = Recipe.from_dict(raw)
    if validate: validate_recipe(recipe)
    return recipe

