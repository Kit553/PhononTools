# Exported from 3_Analysis/src/cat_tools/freqball/plot.py; update via sync_freqball_standalone.py.
# ruff: noqa -- generated copy; lint the shared implementation and exporter.
#!/usr/bin/env python3
"""Plot canonical FreqBall tables with the optional Matplotlib backend."""

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
    weight_column: str = "mobile_participation",
    color_column: str = "direction_score",
    color_range: tuple[float, float] | None = None,
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
            weight_column,
            color_column,
            "trace",
        },
        "spectra": {"frequency_THz", "pdos_per_THz"},
        "mass-summary": {"mass_u", metric},
    }
    if kind not in required:
        raise ValueError(f"Unknown plot kind: {kind}")
    missing = required[kind] - set(frame)
    if missing:
        raise ValueError(f"{source}: missing columns {sorted(missing)}")
    if (ticks is None) != (labels is None) or (ticks is not None and len(ticks) != len(labels)):
        raise ValueError("Supply equal-length --ticks and --labels together")
    limits = color_range or ((-1, 1) if color_column == "direction_score" else (0, 1))
    if len(limits) != 2 or not np.isfinite(limits).all() or limits[0] >= limits[1]:
        raise ValueError("Color range must contain two finite increasing values")

    figure, axes = plt.subplots(figsize=(6.0, 4.5), layout="constrained")
    try:
        if kind == "fatband":
            for _, branch in frame.groupby("trace", sort=False):
                branch = branch.sort_values("distance")
                axes.plot(branch.distance, branch.frequency_THz, color="0.75", linewidth=0.5)
            valid = np.isfinite(frame[color_column]) & np.isfinite(frame[weight_column])
            if not valid.any():
                raise ValueError(
                    "No finite fatband weights/colors: supply direction families or choose "
                    "named projection columns with --weight-column and --color-column"
                )
            modes = frame[valid]
            points = axes.scatter(
                modes.distance,
                modes.frequency_THz,
                s=40 * np.clip(modes[weight_column].to_numpy(), 0, 1),
                c=modes[color_column],
                cmap="coolwarm",
                vmin=limits[0],
                vmax=limits[1],
                linewidths=0,
                alpha=0.85,
            )
            color_label = (
                "Conduction minus tetrahedral fraction"
                if color_column == "direction_score"
                else color_column.replace("_", " ")
            )
            figure.colorbar(points, ax=axes, label=color_label)
            axes.axhline(0, color="0.4", linestyle="--", linewidth=0.6)
            axes.set(xlabel="Wave vector", ylabel="Frequency (THz)")
            if ticks is not None:
                axes.set_xticks(ticks, labels)
        else:
            if kind == "spectra":
                grouping = [
                    key
                    for key in (
                        "system",
                        "mass_u",
                        "equivalence_group",
                        "category",
                        "family",
                        "selection",
                        "projection",
                    )
                    if key in frame
                ]
                x, y = "frequency_THz", "pdos_per_THz"
                axes.set(xlabel="Frequency (THz)", ylabel="Projected DOS (THz⁻¹)")
            else:
                grouping = [
                    key
                    for key in ("system", "equivalence_group", "selection", "projection")
                    if key in frame
                ]
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


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=("fatband", "spectra", "mass-summary"))
    parser.add_argument("table", type=Path)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--equivalence-group", type=int)
    parser.add_argument("--metric", default="average_frequency_THz")
    parser.add_argument("--ticks", nargs="+", type=float)
    parser.add_argument("--labels", nargs="+")
    parser.add_argument("--weight-column", default="mobile_participation")
    parser.add_argument("--color-column", default="direction_score")
    parser.add_argument("--color-range", nargs=2, type=float)
    args = parser.parse_args(argv)
    print(
        plot_table(
            args.table,
            args.output,
            kind=args.kind,
            equivalence_group=args.equivalence_group,
            metric=args.metric,
            ticks=args.ticks,
            labels=args.labels,
            weight_column=args.weight_column,
            color_column=args.color_column,
            color_range=None if args.color_range is None else tuple(args.color_range),
        )
    )


if __name__ == "__main__":
    main()
