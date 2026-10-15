"""Clickable dispersion/DOS explorer with an embedded mode view."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np
from pandas.api.types import is_numeric_dtype

from cat_shared.mode_tables import merge_mode_columns, read_table
from cat_tools.phonopy.modes import calculate_band, calculate_mesh, mode_artifact
from cat_tools.phonopy.spectra import dos_artifact
from cat_tools.phonopy.workflow import projections_from_config

from .geometry import SceneOptions, build_scene
from .model import mode_from_band
from .render import draw_scene


@dataclass(frozen=True)
class ExplorerData:
    band: Any
    band_table: Any
    dos_table: Any
    color_by: str | None
    color_map: str = "viridis"
    color_range: tuple[float, float] | None = None
    color_label: str | None = None


def _join_properties(table, joins, base: Path | None):
    reports = []
    for item in joins:
        path = Path(item["path"]).expanduser()
        if not path.is_absolute():
            if base is None:
                raise ValueError("Relative dispersion joins require a configuration base path")
            path = (base / path).resolve()
        external = read_table(path)
        for column, value in item.get("where", {}).items():
            external = external[
                np.isclose(external[column], value)
                if isinstance(value, float)
                else external[column] == value
            ]
        table, report = merge_mode_columns(
            table,
            external,
            columns=item.get("columns"),
            tolerance=float(item.get("q_tolerance", 1e-7)),
            require_complete=bool(item.get("require_complete", True)),
        )
        reports.append(report)
    return table, tuple(reports)


def prepare_data(
    phonon,
    *,
    system: str,
    dispersion: dict,
    dos: dict,
    base: str | Path | None = None,
) -> ExplorerData:
    """Calculate the shared band eigenvectors, projections, and side DOS."""
    band_options = dict(dispersion)
    color_by = band_options.pop("color_by", None)
    color_map = str(band_options.pop("color_map", "viridis"))
    color_range_value = band_options.pop("color_range", None)
    color_range = (
        None if color_range_value is None else tuple(float(value) for value in color_range_value)
    )
    color_label = band_options.pop("color_label", None)
    joins = list(band_options.pop("joins", []))
    projection_config = band_options.pop("projections", None)
    weighting = band_options.pop("weighting", "eigenvector")
    velocity_directions = band_options.pop("velocity_directions", {})
    projections = projections_from_config(phonon, projection_config)
    band_options["eigenvectors"] = True
    band_options.setdefault(
        "group_velocities", bool(color_by and str(color_by).startswith("group_velocity"))
    )
    band = calculate_band(phonon, **band_options)
    band_table = mode_artifact(
        band,
        system=system,
        symbols=phonon.primitive.symbols,
        masses=phonon.primitive.masses,
        projections=projections,
        weighting=weighting,
        velocity_directions=velocity_directions,
    ).data
    band_table, _ = _join_properties(
        band_table,
        joins,
        None if base is None else Path(base),
    )
    if color_by is not None and color_by not in band_table:
        available = [
            column
            for column in band_table
            if column.endswith("_participation") or column.startswith("group_velocity")
        ]
        raise KeyError(f"Unknown dispersion color {color_by!r}; available: {available}")
    if color_by is not None and not is_numeric_dtype(band_table[color_by]):
        raise TypeError(f"Dispersion color {color_by!r} must be numeric")

    dos_options = dict(dos)
    mesh = dos_options.pop("mesh", [21, 21, 21])
    gamma_center = bool(dos_options.pop("gamma_center", True))
    dos_projections = projections_from_config(
        phonon, dos_options.pop("projections", projection_config)
    )
    calculate_mesh(
        phonon,
        mesh,
        eigenvectors=True,
        mesh_symmetry=False,
        gamma_center=gamma_center,
    )
    density = dos_artifact(
        phonon,
        system=system,
        projections=dos_projections,
        **dos_options,
    ).data
    return ExplorerData(
        band,
        band_table,
        density,
        color_by,
        color_map,
        color_range,
        None if color_label is None else str(color_label),
    )


def _ticks(data) -> tuple[list[float], list[str]]:
    labels = list(data.metadata.get("labels", []))
    connections = list(data.metadata.get("path_connections", []))
    positions: list[float] = []
    names: list[str] = []
    cursor = 0
    for segment in np.unique(data.segments):
        indices = np.flatnonzero(data.segments == segment)
        segment_positions = [float(data.distances[indices[0]]), float(data.distances[indices[-1]])]
        segment_labels = labels[cursor : cursor + 2]
        if len(segment_labels) < 2:
            segment_labels = [f"q{segment}", f"q{segment + 1}"]
        for position, label in zip(segment_positions, segment_labels, strict=True):
            if positions and abs(position - positions[-1]) < 1e-10:
                if label not in names[-1].split("|"):
                    names[-1] += f"|{label}"
            else:
                positions.append(position)
                names.append(str(label))
        cursor += 1 if connections[int(segment)] else 2
    return positions, names


def _dispersion(ax, data: ExplorerData):
    from matplotlib.cm import ScalarMappable
    from matplotlib.collections import LineCollection
    from matplotlib.colors import Normalize

    frequencies = np.asarray(data.band.frequencies)
    values = None
    normalizer = None
    if data.color_by:
        values = data.band_table[data.color_by].to_numpy(float).reshape(frequencies.shape)
        finite = values[np.isfinite(values)]
        if not len(finite):
            raise ValueError(f"Dispersion color {data.color_by!r} contains no finite values")
        limits = data.color_range or (float(finite.min()), float(finite.max()))
        normalizer = Normalize(*limits)
    for segment in np.unique(data.band.segments):
        indices = np.flatnonzero(data.band.segments == segment)
        distance = data.band.distances[indices]
        for band in range(frequencies.shape[1]):
            points = np.column_stack((distance, frequencies[indices, band]))
            if values is None:
                ax.plot(points[:, 0], points[:, 1], color="#242424", linewidth=0.8)
            else:
                line_segments = np.stack((points[:-1], points[1:]), axis=1)
                collection = LineCollection(line_segments, cmap=data.color_map, norm=normalizer)
                collection.set_array((values[indices[:-1], band] + values[indices[1:], band]) / 2)
                collection.set_linewidth(1.4)
                ax.add_collection(collection)
    positions, labels = _ticks(data.band)
    ax.set_xticks(positions, labels)
    for position in positions:
        ax.axvline(position, color="#999999", linewidth=0.6, alpha=0.45)
    ax.axhline(0, color="#777777", linewidth=0.7, linestyle="--")
    ax.set_xlim(float(np.min(data.band.distances)), float(np.max(data.band.distances)))
    ax.set_ylim(float(np.nanmin(frequencies)), float(np.nanmax(frequencies)))
    ax.set_xlabel("Wave vector")
    ax.set_ylabel("Frequency / THz")
    ax.set_title("Click a phonon branch")
    if values is None:
        return None
    return ScalarMappable(norm=normalizer, cmap=data.color_map)


def _dos(ax, table) -> None:
    frequency = table["frequency_THz"]
    columns = [column for column in table if column.endswith("_DOS")]
    for column in columns:
        label = column.removesuffix("_DOS").replace("total", "Total")
        width = 1.5 if column == "total_DOS" else 1.0
        ax.plot(table[column], frequency, linewidth=width, label=label)
    ax.axhline(0, color="#777777", linewidth=0.7, linestyle="--")
    ax.set_xlabel("DOS")
    ax.set_title("Density of states")
    ax.tick_params(labelleft=False)
    ax.legend(frameon=False, fontsize=8)


def _nearest(event, ax, data) -> tuple[int, int]:
    frequencies = np.asarray(data.frequencies)
    coordinates = np.column_stack(
        (
            np.repeat(np.asarray(data.distances), frequencies.shape[1]),
            frequencies.reshape(-1),
        )
    )
    pixels = ax.transData.transform(coordinates)
    closest = int(np.argmin(np.sum((pixels - [event.x, event.y]) ** 2, axis=1)))
    return np.unravel_index(closest, frequencies.shape)


def show_explorer(
    phonon,
    data: ExplorerData,
    options: SceneOptions,
    *,
    system: str,
    initial_q_index: int = 0,
    initial_band: int = 0,
    animation_interval_ms: int = 80,
) -> None:
    """Open the dispersion, DOS and animated-mode dashboard."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError(
            "Install PhonoVis with: pip install 'cat-tooling[visualization]'"
        ) from exc

    frequencies = np.asarray(data.band.frequencies)
    if initial_q_index not in range(len(frequencies)) or initial_band not in range(
        frequencies.shape[1]
    ):
        raise IndexError("Initial q-index or band is out of range")
    mode = mode_from_band(data.band, phonon, initial_q_index, initial_band, system=system)
    scene = build_scene(mode, options)
    state = {
        "q_index": initial_q_index,
        "band": initial_band,
        "phase": np.deg2rad(options.phase_deg),
        "playing": True,
        "options": options,
    }
    figure = plt.figure(figsize=(15, 6.2))
    grid = figure.add_gridspec(1, 3, width_ratios=(4.2, 1.25, 4.1), wspace=0.08)
    band_axis = figure.add_subplot(grid[0, 0])
    dos_axis = figure.add_subplot(grid[0, 1], sharey=band_axis)
    mode_axis = figure.add_subplot(grid[0, 2], projection="3d")
    mappable = _dispersion(band_axis, data)
    _dos(dos_axis, data.dos_table)
    if mappable is not None:
        figure.colorbar(
            mappable,
            ax=band_axis,
            label=data.color_label or data.color_by,
            pad=0.02,
        )
    marker = band_axis.plot(
        data.band.distances[initial_q_index],
        frequencies[initial_q_index, initial_band],
        marker="o",
        color="#D94F45",
        markersize=6,
        linestyle="none",
        zorder=20,
    )[0]

    def selected_mode():
        return mode_from_band(
            data.band,
            phonon,
            state["q_index"],
            state["band"],
            system=system,
        )

    def redraw() -> None:
        selected = selected_mode()
        draw_scene(
            mode_axis,
            selected,
            scene,
            state["options"],
            phase_rad=state["phase"],
            displaced_atoms=True,
        )
        marker.set_data(
            [data.band.distances[state["q_index"]]],
            [frequencies[state["q_index"], state["band"]]],
        )
        figure.canvas.draw_idle()

    redraw()
    figure.text(
        0.5,
        0.015,
        "click: select mode   space: pause   arrows: q/band   +/-: amplitude   L: labels   A: arrows",
        ha="center",
        fontsize=8,
        color="#555555",
    )

    def click(event) -> None:
        if event.inaxes is band_axis:
            state["q_index"], state["band"] = _nearest(event, band_axis, data.band)
            redraw()

    def key(event) -> None:
        if event.key in {" ", "space"}:
            state["playing"] = not state["playing"]
        elif event.key == "left":
            state["q_index"] = max(0, state["q_index"] - 1)
        elif event.key == "right":
            state["q_index"] = min(len(frequencies) - 1, state["q_index"] + 1)
        elif event.key == "down":
            state["band"] = max(0, state["band"] - 1)
        elif event.key == "up":
            state["band"] = min(frequencies.shape[1] - 1, state["band"] + 1)
        elif event.key in {"+", "="}:
            state["options"] = replace(
                state["options"], amplitude_A=1.2 * state["options"].amplitude_A
            )
        elif event.key == "-":
            state["options"] = replace(
                state["options"], amplitude_A=state["options"].amplitude_A / 1.2
            )
        elif event.key and event.key.lower() == "l":
            state["options"] = replace(
                state["options"], show_labels=not state["options"].show_labels
            )
        elif event.key and event.key.lower() == "a":
            state["options"] = replace(
                state["options"], show_arrows=not state["options"].show_arrows
            )
        else:
            return
        redraw()

    def animate() -> None:
        if state["playing"]:
            state["phase"] = (state["phase"] + 0.16) % (2 * np.pi)
            redraw()

    figure.canvas.mpl_connect("button_press_event", click)
    figure.canvas.mpl_connect("key_press_event", key)
    timer = figure.canvas.new_timer(interval=animation_interval_ms)
    timer.add_callback(animate)
    timer.start()
    plt.show()
