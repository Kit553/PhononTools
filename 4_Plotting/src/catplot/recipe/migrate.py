"""Schema migration entry point reserved for future versions."""

from __future__ import annotations

from catplot.errors import RecipeError


def migrate_recipe(raw, target=1):
    current = int(raw.get("schema", 1))
    if current == target: return raw
    raise RecipeError(f"No migration path from schema {current} to {target}")

