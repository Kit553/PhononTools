"""Context passed to every plot renderer."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass
class RenderContext:
    figure: Any
    axes: Any
    panel_id: str
    layer_index: int
    theme: dict[str, Any]
    recipe: Any
    sources: dict[str, Any]

