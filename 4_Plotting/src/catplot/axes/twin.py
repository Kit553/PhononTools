"""Lazy right-y-axis creation."""

from __future__ import annotations


def get_axis_target(ax, side: str, cache: dict[int, object]):
    if side == "left":
        return ax
    if side != "right":
        raise ValueError("Layer axis must be 'left' or 'right'")
    key = id(ax)
    if key not in cache:
        cache[key] = ax.twinx()
    return cache[key]

