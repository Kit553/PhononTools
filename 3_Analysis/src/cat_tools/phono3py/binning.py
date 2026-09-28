"""Small weighted reductions for dense reciprocal-mesh mode clouds."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from cat_shared.artifacts import TableArtifact


@dataclass(frozen=True)
class BinSpec:
    """Reproducible linear or logarithmic bin edges."""

    start: float | None = None
    stop: float | None = None
    width: float | None = None
    count: int | None = None
    scale: str = "linear"

    @classmethod
    def parse(cls, value: int | Sequence[float] | Mapping[str, Any]) -> BinSpec | np.ndarray:
        if isinstance(value, int):
            return cls(count=value)
        if isinstance(value, Mapping):
            return cls(**value)
        edges = np.asarray(value, dtype=float)
        if edges.ndim != 1 or len(edges) < 2 or np.any(np.diff(edges) <= 0):
            raise ValueError("Explicit bin edges must be a strictly increasing sequence")
        return edges

    def edges(self, values: np.ndarray) -> np.ndarray:
        finite = np.asarray(values, dtype=float)
        finite = finite[np.isfinite(finite)]
        if self.scale == "log":
            finite = finite[finite > 0]
        if not len(finite):
            raise ValueError("Cannot construct bins from empty data")
        start = float(finite.min()) if self.start is None else float(self.start)
        stop = float(finite.max()) if self.stop is None else float(self.stop)
        if stop <= start:
            raise ValueError("Bin stop must be larger than start")
        if self.scale == "log":
            if start <= 0:
                raise ValueError("Logarithmic bins require a positive start")
            count = self.count or 80
            return np.geomspace(start, stop, count + 1)
        if self.scale != "linear":
            raise ValueError("Bin scale must be 'linear' or 'log'")
        if self.width is not None:
            if self.width <= 0:
                raise ValueError("Bin width must be positive")
            count = max(1, int(np.ceil((stop - start) / self.width)))
            return start + np.arange(count + 1) * self.width
        return np.linspace(start, stop, (self.count or 120) + 1)


def bin_edges(values: np.ndarray, spec: int | Sequence[float] | Mapping[str, Any]) -> np.ndarray:
    parsed = BinSpec.parse(spec)
    return parsed if isinstance(parsed, np.ndarray) else parsed.edges(values)


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    order = np.argsort(values, kind="mergesort")
    values, weights = values[order], weights[order]
    cumulative = np.cumsum(weights) - 0.5 * weights
    return float(np.interp(quantile * weights.sum(), cumulative, values))


def _groups(data: pd.DataFrame, requested: Sequence[str] | None) -> list[str]:
    defaults = ("system", "dataset", "temperature_K", "boundary_mfp_um")
    return [column for column in (requested or defaults) if column in data]


def _units(source: TableArtifact, value: str) -> dict[str, str]:
    units = {
        "frequency_THz": source.units.get("frequency_THz", "THz"),
        "frequency_bin_left_THz": source.units.get("frequency_THz", "THz"),
        "frequency_bin_right_THz": source.units.get("frequency_THz", "THz"),
        "frequency_bin_width_THz": source.units.get("frequency_THz", "THz"),
    }
    if value in source.units:
        units.update(
            {
                name: source.units[value]
                for name in ("value_mean", "value_median", "value_q25", "value_q75", "value")
            }
        )
    return units


def profile_artifact(
    source: TableArtifact,
    *,
    name: str,
    value: str,
    bins: int | Sequence[float] | Mapping[str, Any],
    groups: Sequence[str] | None = None,
    weight: str = "grid_weight",
) -> TableArtifact:
    """Frequency-binned weighted mean, median, and interquartile range."""
    data = source.data
    if value not in data:
        raise KeyError(value)
    edges = bin_edges(data["frequency_THz"].to_numpy(float), bins)
    group_columns = _groups(data, groups)
    rows: list[dict[str, Any]] = []
    iterator = (
        data.groupby(group_columns, dropna=False, sort=False) if group_columns else [((), data)]
    )
    for group_key, group in iterator:
        keys = group_key if isinstance(group_key, tuple) else (group_key,)
        identity = dict(zip(group_columns, keys, strict=True))
        x = group["frequency_THz"].to_numpy(float)
        y = group[value].to_numpy(float)
        w = group[weight].to_numpy(float) if weight in group else np.ones(len(group))
        bin_index = np.searchsorted(edges, x, side="right") - 1
        bin_index[x == edges[-1]] = len(edges) - 2
        for index in range(len(edges) - 1):
            mask = (bin_index == index) & np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
            values, weights = y[mask], w[mask]
            row = {
                **identity,
                "frequency_bin": index,
                "frequency_THz": 0.5 * (edges[index] + edges[index + 1]),
                "frequency_bin_left_THz": edges[index],
                "frequency_bin_right_THz": edges[index + 1],
                "frequency_bin_width_THz": edges[index + 1] - edges[index],
                "mode_count": int(mask.sum()),
                "weight_sum": float(weights.sum()),
                "value_mean": np.nan,
                "value_median": np.nan,
                "value_q25": np.nan,
                "value_q75": np.nan,
            }
            if len(values):
                row.update(
                    value_mean=float(np.average(values, weights=weights)),
                    value_median=_weighted_quantile(values, weights, 0.5),
                    value_q25=_weighted_quantile(values, weights, 0.25),
                    value_q75=_weighted_quantile(values, weights, 0.75),
                )
            rows.append(row)
    data = pd.DataFrame(rows)
    keys = (*group_columns, "frequency_bin")
    return TableArtifact(
        name,
        data,
        keys,
        {column: unit for column, unit in _units(source, value).items() if column in data},
        {"property": value, "weight": weight},
    )


def density_artifact(
    source: TableArtifact,
    *,
    name: str,
    value: str,
    frequency_bins: int | Sequence[float] | Mapping[str, Any],
    value_bins: int | Sequence[float] | Mapping[str, Any],
    groups: Sequence[str] | None = None,
    weight: str = "grid_weight",
) -> TableArtifact:
    """Complete 2D frequency-property histogram for a CATPlot heatmap."""
    data = source.data
    xedges = bin_edges(data["frequency_THz"].to_numpy(float), frequency_bins)
    yedges = bin_edges(data[value].to_numpy(float), value_bins)
    group_columns = _groups(data, groups)
    rows: list[pd.DataFrame] = []
    iterator = (
        data.groupby(group_columns, dropna=False, sort=False) if group_columns else [((), data)]
    )
    for group_key, group in iterator:
        keys = group_key if isinstance(group_key, tuple) else (group_key,)
        identity = dict(zip(group_columns, keys, strict=True))
        x = group["frequency_THz"].to_numpy(float)
        y = group[value].to_numpy(float)
        w = group[weight].to_numpy(float) if weight in group else np.ones(len(group))
        finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
        weighted, _, _ = np.histogram2d(x[finite], y[finite], (xedges, yedges), weights=w[finite])
        counts, _, _ = np.histogram2d(x[finite], y[finite], (xedges, yedges))
        total = weighted.sum()
        area = np.diff(xedges)[:, None] * np.diff(yedges)[None, :]
        ix, iy = np.meshgrid(np.arange(len(xedges) - 1), np.arange(len(yedges) - 1), indexing="ij")
        frame = pd.DataFrame(
            {
                **{column: np.repeat(value_, ix.size) for column, value_ in identity.items()},
                "frequency_bin": ix.ravel(),
                "value_bin": iy.ravel(),
                "frequency_THz": ((xedges[:-1] + xedges[1:]) / 2)[ix].ravel(),
                "value": ((yedges[:-1] + yedges[1:]) / 2)[iy].ravel(),
                "mode_count": counts.ravel().astype(int),
                "weight_sum": weighted.ravel(),
                "weight_fraction": weighted.ravel() / total if total else np.nan,
                "probability_density": weighted.ravel() / (total * area).ravel()
                if total
                else np.nan,
            }
        )
        rows.append(frame)
    keys = (*group_columns, "frequency_bin", "value_bin")
    data = pd.concat(rows, ignore_index=True)
    return TableArtifact(
        name,
        data,
        keys,
        {column: unit for column, unit in _units(source, value).items() if column in data},
        {
            "property": value,
            "weight": weight,
            "frequency_edges": xedges.tolist(),
            "value_edges": yedges.tolist(),
        },
    )


def spectrum_artifact(
    source: TableArtifact,
    *,
    name: str,
    value: str,
    bins: int | Sequence[float] | Mapping[str, Any],
    groups: Sequence[str] | None = None,
    weighting: str = "none",
) -> TableArtifact:
    """Frequency spectrum and cumulative curve for an additive mode property."""
    data = source.data
    edges = bin_edges(data["frequency_THz"].to_numpy(float), bins)
    group_columns = _groups(data, groups)
    rows: list[pd.DataFrame] = []
    iterator = (
        data.groupby(group_columns, dropna=False, sort=False) if group_columns else [((), data)]
    )
    for group_key, group in iterator:
        keys = group_key if isinstance(group_key, tuple) else (group_key,)
        identity = dict(zip(group_columns, keys, strict=True))
        x = group["frequency_THz"].to_numpy(float)
        y = group[value].to_numpy(float)
        if weighting in {"bz", "multiplicity"}:
            if "grid_weight" not in group or "ir_index" not in group:
                raise KeyError("BZ/multiplicity weighting needs grid_weight and ir_index")
            y = y * group["grid_weight"].to_numpy(float)
            if weighting == "bz":
                denominator = group[["ir_index", "grid_weight"]].drop_duplicates()[
                    "grid_weight"
                ].sum()
                y = y / denominator
        elif weighting != "none":
            raise ValueError("Spectrum weighting must be 'none', 'multiplicity', or 'bz'")
        finite = np.isfinite(x) & np.isfinite(y)
        summed, _ = np.histogram(x[finite], edges, weights=y[finite])
        counts, _ = np.histogram(x[finite], edges)
        cumulative = np.cumsum(summed)
        total = summed.sum()
        frame = pd.DataFrame(
            {
                **{column: np.repeat(value_, len(summed)) for column, value_ in identity.items()},
                "frequency_bin": np.arange(len(summed)),
                "frequency_THz": (edges[:-1] + edges[1:]) / 2,
                "frequency_bin_left_THz": edges[:-1],
                "frequency_bin_right_THz": edges[1:],
                "frequency_bin_width_THz": np.diff(edges),
                "mode_count": counts,
                "value_sum": summed,
                "value_density": summed / np.diff(edges),
                "value_cumulative": cumulative,
                "value_cumulative_fraction": cumulative / total if total else np.nan,
            }
        )
        rows.append(frame)
    value_unit = source.units.get(value)
    units = _units(source, value)
    if value_unit:
        units.update({column: value_unit for column in ("value_sum", "value_cumulative")})
        units["value_density"] = f"{value_unit}/THz"
    data = pd.concat(rows, ignore_index=True)
    return TableArtifact(
        name,
        data,
        (*group_columns, "frequency_bin"),
        {column: unit for column, unit in units.items() if column in data},
        {"property": value, "weighting": weighting, "frequency_edges": edges.tolist()},
    )


def accumulation_artifact(
    source: TableArtifact,
    *,
    name: str,
    value: str,
    coordinate: str,
    bins: int | Sequence[float] | Mapping[str, Any],
    groups: Sequence[str] | None = None,
    weight: str | None = "grid_weight",
) -> TableArtifact:
    """Cumulative additive mode property along frequency, MFP, or another coordinate."""
    data = source.data
    for column in (value, coordinate):
        if column not in data:
            raise KeyError(column)
    edges = bin_edges(data[coordinate].to_numpy(float), bins)
    group_columns = _groups(data, groups)
    rows: list[pd.DataFrame] = []
    iterator = (
        data.groupby(group_columns, dropna=False, sort=False) if group_columns else [((), data)]
    )
    for group_key, group in iterator:
        keys = group_key if isinstance(group_key, tuple) else (group_key,)
        identity = dict(zip(group_columns, keys, strict=True))
        x = group[coordinate].to_numpy(float)
        y = group[value].to_numpy(float)
        if weight is not None:
            if weight not in group:
                raise KeyError(weight)
            y = y * group[weight].to_numpy(float)
        finite = np.isfinite(x) & np.isfinite(y)
        summed, _ = np.histogram(x[finite], edges, weights=y[finite])
        counts, _ = np.histogram(x[finite], edges)
        cumulative = np.cumsum(summed)
        total = summed.sum()
        rows.append(
            pd.DataFrame(
                {
                    **{
                        column: np.repeat(value_, len(summed))
                        for column, value_ in identity.items()
                    },
                    "coordinate_bin": np.arange(len(summed)),
                    "coordinate": edges[1:],
                    "coordinate_bin_left": edges[:-1],
                    "coordinate_bin_right": edges[1:],
                    "coordinate_bin_width": np.diff(edges),
                    "mode_count": counts,
                    "value_sum": summed,
                    "value_density": summed / np.diff(edges),
                    "value_cumulative": cumulative,
                    "value_cumulative_fraction": cumulative / total if total else np.nan,
                }
            )
        )
    coordinate_unit = source.units.get(coordinate)
    value_unit = source.units.get(value)
    units = {}
    if coordinate_unit:
        units.update(
            {
                column: coordinate_unit
                for column in (
                    "coordinate",
                    "coordinate_bin_left",
                    "coordinate_bin_right",
                    "coordinate_bin_width",
                )
            }
        )
    if value_unit:
        units.update(
            {
                "value_sum": value_unit,
                "value_cumulative": value_unit,
                "value_density": (
                    f"{value_unit}/({coordinate_unit})" if coordinate_unit else value_unit
                ),
            }
        )
    return TableArtifact(
        name,
        pd.concat(rows, ignore_index=True),
        (*group_columns, "coordinate_bin"),
        units,
        {
            "property": value,
            "coordinate": coordinate,
            "coordinate_edges": edges.tolist(),
            "weight": weight,
        },
    )
