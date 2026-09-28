"""Dotted-path overrides used by the CLI."""

from __future__ import annotations

from copy import deepcopy
import yaml

from catplot.errors import RecipeError
from .model import Recipe


def _set_path(root, path, value):
    parts = path.split("."); current = root
    for part in parts[:-1]:
        if isinstance(current, list):
            try: current = current[int(part)]
            except (ValueError, IndexError) as exc: raise RecipeError(f"Invalid override path {path!r}") from exc
        elif isinstance(current, dict):
            if part not in current: current[part] = {}
            current = current[part]
        else: raise RecipeError(f"Override path {path!r} crosses a scalar value")
    final = parts[-1]
    if isinstance(current, list): current[int(final)] = value
    elif isinstance(current, dict): current[final] = value
    else: raise RecipeError(f"Invalid override path {path!r}")


def apply_overrides(recipe, overrides):
    raw = deepcopy(recipe.to_dict())
    for override in overrides:
        if "=" not in override: raise RecipeError("Overrides must use path=value")
        path, text = override.split("=", 1); _set_path(raw, path, yaml.safe_load(text))
    return Recipe.from_dict(raw)

