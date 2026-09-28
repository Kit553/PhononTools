"""Inset axes creation."""

from __future__ import annotations


def add_inset(ax, bounds, *, projection="2d"):
    kwargs = {"projection": "3d"} if projection == "3d" else {}
    return ax.inset_axes(bounds, **kwargs)

