"""Create Matplotlib panels with controlled physical geometry."""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt

from catplot.style.geometry import DEFAULT_EDGE_MM, mm_to_inches
from catplot.style.manager import apply_style


def new_panel(
    *,
    edge_mm: float = DEFAULT_EDGE_MM,
    profile: str = "paper",
    theme: str = "default",
    square_axes: bool = False,
    **subplot_kwargs: Any,
):
    """Return ``(figure, axes, theme)`` for one exact-square PDF artboard.

    ``square_axes`` concerns the inner data box. The outer artboard remains
    square regardless of this setting.
    """
    selected_theme = apply_style(profile=profile, theme=theme)
    edge_in = mm_to_inches(edge_mm)
    fig, ax = plt.subplots(
        figsize=(edge_in, edge_in),
        layout="constrained",
        **subplot_kwargs,
    )
    if square_axes:
        ax.set_box_aspect(1)
    return fig, ax, selected_theme

