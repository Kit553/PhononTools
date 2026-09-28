"""Per-layer render record retained for inspection and tests."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass
class LayerRecord:
    panel_id: str
    layer_index: int
    plot_type: str
    axis: str
    result: Any

