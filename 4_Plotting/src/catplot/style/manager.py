"""Stable style and square-panel helpers.

The rest of the plotting backend should call this module rather than mutating
Matplotlib rcParams ad hoc. Plot-specific options belong to recipes; global
graphical policy belongs here and in the packaged style/theme files.
"""

from __future__ import annotations

from importlib.resources import files
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt

from .palette import apply_palette, load_theme


def _style_dir():
    return files("catplot.style").joinpath("sheets")


def available_profiles() -> tuple[str, ...]:
    """Return installed named style profiles, excluding the internal base."""
    names = sorted(
        item.name.removesuffix(".mplstyle")
        for item in _style_dir().iterdir()
        if item.name.endswith(".mplstyle") and item.name != "base.mplstyle"
    )
    return tuple(names)


def apply_style(
    profile: str = "paper",
    theme: str = "default",
) -> dict[str, Any]:
    """Reset Matplotlib and apply one profile plus one semantic theme.

    Resetting first makes rendering independent of a user's matplotlibrc file or
    of style changes made by an earlier plot in the same Python process.
    """
    if profile not in available_profiles():
        choices = ", ".join(available_profiles())
        raise ValueError(f"Unknown profile {profile!r}. Available profiles: {choices}")

    mpl.rcdefaults()
    base = _style_dir().joinpath("base.mplstyle")
    selected = _style_dir().joinpath(f"{profile}.mplstyle")
    plt.style.use([str(base), str(selected)])

    selected_theme = load_theme(theme)
    apply_palette(selected_theme)
    return selected_theme
