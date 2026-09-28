"""Colorbar placement."""

from __future__ import annotations


def add_colorbar(fig, ax, mappable, spec: dict | bool | None = None, *, label: str | None = None):
    if spec is False or mappable is None:
        return None
    options = {} if spec is None or spec is True else dict(spec)
    if label is not None and "label" not in options:
        options["label"] = label
    return fig.colorbar(mappable, ax=ax, **options)
