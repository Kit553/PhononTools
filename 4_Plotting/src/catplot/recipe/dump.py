"""Stable human-readable YAML recipe writing."""

from __future__ import annotations

from pathlib import Path
import yaml


def dump_recipe(recipe, path):
    target = Path(path); target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="\n") as stream:
        yaml.safe_dump(recipe.to_dict(), stream, sort_keys=False, allow_unicode=True, width=100)
    return target

