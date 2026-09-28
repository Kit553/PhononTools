"""Selective rasterization helpers for dense artists."""

from __future__ import annotations


def rasterize_if_large(artist, count: int, *, threshold: int = 20_000) -> bool:
    if count >= threshold and hasattr(artist, "set_rasterized"):
        artist.set_rasterized(True)
        return True
    return False

