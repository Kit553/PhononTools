"""Load and apply named colour themes."""

from __future__ import annotations

from importlib.resources import files
from typing import Any

import matplotlib as mpl
from cycler import cycler
import yaml


def _theme_dir():
    return files("catplot.style").joinpath("themes")


def available_themes() -> tuple[str, ...]:
    """Return installed colour themes in deterministic order."""
    return tuple(
        sorted(
            item.name.removesuffix(".yaml")
            for item in _theme_dir().iterdir()
            if item.name.endswith(".yaml")
        )
    )


def load_theme(name: str = "default") -> dict[str, Any]:
    """Load a theme and verify its minimum required structure."""
    target = _theme_dir().joinpath(f"{name}.yaml")
    if not target.is_file():
        choices = ", ".join(available_themes())
        raise ValueError(f"Unknown theme {name!r}. Available themes: {choices}")
    with target.open("r", encoding="utf-8") as stream:
        theme = yaml.safe_load(stream)
    if not isinstance(theme, dict) or not theme.get("series"):
        raise ValueError(f"Theme {name!r} must define a non-empty 'series' list")
    if not isinstance(theme.get("semantic"), dict):
        raise ValueError(f"Theme {name!r} must define a 'semantic' mapping")
    return theme


def apply_palette(theme: dict[str, Any]) -> None:
    """Install a validated theme's categorical cycle into Matplotlib."""
    mpl.rcParams["axes.prop_cycle"] = cycler(color=theme["series"])

