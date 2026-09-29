#!/usr/bin/env python3
"""Plot FreqBall tables with Matplotlib; neither CAT nor CATPlot is required."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def plot_table(
    table: str | Path,
    output: str | Path,
    *,
    kind: str,
    equivalence_group: int | None = None,
    metric: str = "average_frequency_THz",
    ticks: list[float] | None = None,
    labels: list[str] | None = None,
) -> Path:
    """Render the same canonical CSV/Parquet outputs as the packaged workflow."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    source = Path(table).expanduser()
    frame = pd.read_parquet(source) if source.suffix == ".parquet" else pd.read_csv(source)
    if equivalence_group is not None:
        if "equivalence_group" not in frame:
            raise ValueError("This table does not contain equivalence groups")
        frame = frame[frame.equivalence_group == equivalence_group]
    if frame.empty:
        raise ValueError("No rows remain for the requested plot")
    required = {
        "fatband": {
            "distance",
            "frequency_THz",
            "mobile_participation",
            "direction_score",
            "trace",
        },
        "spectra": {"frequency_THz", "pdos_per_THz"},
        "mass-summary": {"mass_u", metric, "selection"},
    }
    if kind not in required:
        raise ValueError(f"Unknown plot kind: {kind}")
    missing = required[kind] - set(frame)
    if missing:
        raise ValueError(f"{source}: missing columns {sorted(missing)}")
    if (ticks is None) != (labels is None) or (ticks is not None and len(ticks) != len(labels)):
        raise ValueError("Supply equal-length --ticks and --labels together")

    figure, axes = plt.subplots(figsize=(6.0, 4.5), layout="constrained")
    try:
        if kind == "fatband":
            for _, branch in frame.groupby("trace", sort=False):
                branch = branch.sort_values("distance")
                axes.plot(branch.distance, branch.frequency_THz, color="0.75", linewidth=0.5)
            valid = np.isfinite(frame.direction_score) & np.isfinite(frame.mobile_participation)
            if not valid.any():
                raise ValueError("No classified modes: supply conduction and tetrahedral families")
            modes = frame[valid]
            points = axes.scatter(
                modes.distance,
                modes.frequency_THz,
                s=40 * np.clip(modes.mobile_participation.to_numpy(), 0, 1),
                c=modes.direction_score,
                cmap="coolwarm",
                vmin=-1,
                vmax=1,
                linewidths=0,
                alpha=0.85,
            )
            figure.colorbar(points, ax=axes, label="Conduction minus tetrahedral fraction")
            axes.axhline(0, color="0.4", linestyle="--", linewidth=0.6)
            axes.set(xlabel="Wave vector", ylabel="Frequency (THz)")
            if ticks is not None:
                axes.set_xticks(ticks, labels)
        else:
            if kind == "spectra":
                grouping = [
                    key
                    for key in ("mass_u", "equivalence_group", "category", "family", "selection")
                    if key in frame
                ]
                x, y = "frequency_THz", "pdos_per_THz"
                axes.set(xlabel="Frequency (THz)", ylabel="Projected DOS (THz⁻¹)")
            else:
                grouping = [key for key in ("equivalence_group", "selection") if key in frame]
                x, y = "mass_u", metric
                axes.set(xlabel="Assigned mass (u)", ylabel=metric.replace("_", " "))
            groups = frame.groupby(grouping, sort=False) if grouping else [("spectrum", frame)]
            for key, values in groups:
                values = values.sort_values(x)
                label = ", ".join(map(str, key)) if isinstance(key, tuple) else str(key)
                axes.plot(
                    values[x],
                    values[y],
                    label=label,
                    linewidth=1.2,
                    marker="o" if kind == "mass-summary" else None,
                    markersize=3,
                )
            axes.legend(fontsize="small", frameon=False)
        target = Path(output).expanduser().resolve()
        target.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(target, dpi=200)
    finally:
        plt.close(figure)
    return target


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=("fatband", "spectra", "mass-summary"))
    parser.add_argument("table", type=Path)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--equivalence-group", type=int)
    parser.add_argument("--metric", default="average_frequency_THz")
    parser.add_argument("--ticks", nargs="+", type=float)
    parser.add_argument("--labels", nargs="+")
    args = parser.parse_args()
    print(
        plot_table(
            args.table,
            args.output,
            kind=args.kind,
            equivalence_group=args.equivalence_group,
            metric=args.metric,
            ticks=args.ticks,
            labels=args.labels,
        )
    )


if __name__ == "__main__":
    main()
