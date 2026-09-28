"""Layout result shared by render pipelines."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass
class LayoutResult:
    figure: Any
    axes: dict[str, Any]
    rows: int
    columns: int

