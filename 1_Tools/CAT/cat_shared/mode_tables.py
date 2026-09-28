"""Strict joins for tables whose rows describe reciprocal-space modes."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class ModeJoinReport:
    """Coverage and columns contributed by one external mode table."""

    columns: tuple[str, ...]
    keys: tuple[str, ...]
    matched_rows: int
    total_rows: int

    @property
    def coverage(self) -> float:
        return self.matched_rows / self.total_rows if self.total_rows else 1.0


def read_table(path: str | Path) -> pd.DataFrame:
    """Read either supported canonical table representation."""
    source = Path(path)
    return pd.read_parquet(source) if source.suffix.lower() == ".parquet" else pd.read_csv(source)


def periodic_qkeys(frame: pd.DataFrame, tolerance: float = 1e-7) -> pd.DataFrame:
    """Integer keys for periodically equivalent reduced reciprocal coordinates."""
    if tolerance <= 0:
        raise ValueError("q-point tolerance must be positive")
    required = ["qx", "qy", "qz"]
    if not set(required).issubset(frame):
        raise KeyError("Mode joins require qx, qy, and qz")
    qpoints = np.mod(frame[required].to_numpy(float) + 0.5, 1.0) - 0.5
    return pd.DataFrame(
        np.rint(qpoints / tolerance).astype(np.int64),
        columns=["_qa", "_qb", "_qc"],
        index=frame.index,
    )


def _join_keys(
    left: pd.DataFrame, right: pd.DataFrame, tolerance: float
) -> tuple[list[str], pd.DataFrame, pd.DataFrame]:
    if "band" not in left or "band" not in right:
        raise KeyError("Mode joins require a zero-based band column")
    path_keys = ["segment", "point_index", "band"]
    if set(path_keys).issubset(left) and set(path_keys).issubset(right):
        left, right = left.copy(), right.copy()
        if {"qx", "qy", "qz"}.issubset(left) and {"qx", "qy", "qz"}.issubset(right):
            left[["_qa", "_qb", "_qc"]] = periodic_qkeys(left, tolerance)
            right[["_qa", "_qb", "_qc"]] = periodic_qkeys(right, tolerance)
            path_keys += ["_qa", "_qb", "_qc"]
        return path_keys, left, right
    left, right = left.copy(), right.copy()
    left[["_qa", "_qb", "_qc"]] = periodic_qkeys(left, tolerance)
    right[["_qa", "_qb", "_qc"]] = periodic_qkeys(right, tolerance)
    return ["_qa", "_qb", "_qc", "band"], left, right


def merge_mode_columns(
    modes: pd.DataFrame,
    extra: pd.DataFrame,
    *,
    columns: Sequence[str] | None = None,
    tolerance: float = 1e-7,
    require_complete: bool = True,
) -> tuple[pd.DataFrame, ModeJoinReport]:
    """Join external properties without silently accepting missing mode rows.

    Band paths use their exact ``segment, point_index, band`` identity when it is
    available in both tables. Other data use periodic q-point keys plus band.
    ``system`` and ``temperature_K`` are included when both sides provide them.
    """
    keys, left, right = _join_keys(modes, extra, tolerance)
    for optional in ("system", "temperature_K"):
        if optional in left and optional in right:
            keys.append(optional)

    identity = {
        "system",
        "dataset",
        "temperature_K",
        "q_index",
        "ir_index",
        "grid_point",
        "segment",
        "point_index",
        "trace",
        "band",
        "qx",
        "qy",
        "qz",
        "distance",
        "frequency_THz",
        "weight",
        "grid_weight",
        "_qa",
        "_qb",
        "_qc",
    }
    selected = list(columns) if columns is not None else sorted(set(right) - identity)
    missing_columns = set(selected) - set(right)
    if missing_columns:
        raise KeyError(f"External mode table has no columns {sorted(missing_columns)}")
    collisions = (set(selected) & set(left)) - set(keys)
    if collisions:
        raise ValueError(f"Mode join would overwrite columns {sorted(collisions)}")
    if right.duplicated(keys, keep=False).any():
        raise ValueError(f"External mode table has non-unique join keys {keys}")

    right = right[keys + selected].copy()
    right["_mode_joined"] = True
    left["_mode_row"] = np.arange(len(left))
    result = left.merge(right, on=keys, how="left", validate="many_to_one", sort=False)
    result.sort_values("_mode_row", inplace=True)
    matched = int(result["_mode_joined"].fillna(False).sum())
    report = ModeJoinReport(tuple(selected), tuple(keys), matched, len(result))
    if require_complete and matched != len(result):
        raise ValueError(
            f"External mode table matched {matched}/{len(result)} rows "
            f"({report.coverage:.1%}); set require_complete: false to allow a partial join"
        )
    if not matched and len(result):
        raise ValueError("External mode table did not match any mode rows")
    result.drop(
        columns=[
            column
            for column in ("_qa", "_qb", "_qc", "_mode_joined", "_mode_row")
            if column in result
        ],
        inplace=True,
    )
    result.reset_index(drop=True, inplace=True)
    return result, report
