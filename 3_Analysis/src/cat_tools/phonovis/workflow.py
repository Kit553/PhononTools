"""Two-command PhonoVis workflow: interactive exploration or deterministic export."""

from __future__ import annotations

import argparse
import re
from dataclasses import asdict
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path
from typing import Any

import yaml

from cat_tools import __version__
from cat_tools.phonopy.input import load_phonon, resolve_files

from .explore import prepare_data, show_explorer
from .geometry import SceneOptions, build_scene, scene_options
from .model import Mode, commensurate_repeat, mode_from_band, mode_from_qpoint
from .render import render_figure, render_gif, show_figure
from .vmd import write_vmd


def _load(filename: str | Path):
    config_path = Path(filename).expanduser().resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    if not isinstance(config, dict):
        raise TypeError("The PhonoVis YAML must contain a mapping")
    system = dict(config.get("system", {}))
    root = Path(system.pop("path", ".")).expanduser()
    root = root.resolve() if root.is_absolute() else (config_path.parent / root).resolve()
    files = resolve_files(root, **system)
    return config_path, config, files, load_phonon(files, log_level=int(config.get("log_level", 0)))


def _mode(phonon, system: str, config: dict, qpoint=None, band: int | None = None) -> Mode:
    values = dict(config.get("mode", {}))
    selected_qpoint = qpoint if qpoint is not None else values.get("qpoint", [0, 0, 0])
    selected_band = int(values.get("band", 0) if band is None else band)
    return mode_from_qpoint(phonon, selected_qpoint, selected_band, system=system)


def _scene(
    mode: Mode,
    config: dict,
    *,
    amplitude_A: float | None = None,
    repeat: tuple[int, int, int] | None = None,
) -> SceneOptions:
    values = dict(config.get("scene", {}))
    from_qpoint = bool(values.pop("repeat_from_qpoint", False))
    denominator = int(values.pop("max_denominator", 24))
    if from_qpoint:
        values["repeat"] = commensurate_repeat(mode.qpoint, max_denominator=denominator)
    if repeat is not None:
        values["repeat"] = repeat
    if amplitude_A is not None:
        values["amplitude_A"] = amplitude_A
    return scene_options(values)


def _stem(mode: Mode) -> str:
    system = re.sub(r"[^A-Za-z0-9_.-]+", "-", mode.system).strip("-") or "mode"
    qpoint = "_".join(f"{value:+.4f}".replace("+", "p").replace("-", "m") for value in mode.qpoint)
    return f"{system}_q-{qpoint}_b{mode.band:03d}"


def _enabled(config: bool | dict | None) -> tuple[bool, dict[str, Any]]:
    if isinstance(config, dict):
        values = dict(config)
        return bool(values.pop("enabled", True)), values
    return bool(config), {}


def _write_manifest(directory, config_path, files, mode, options, outputs) -> Path:
    target = directory / "manifest.yaml"
    manifest = {
        "schema": 1,
        "created_utc": datetime.now(UTC).isoformat(),
        "run": {
            **files.metadata(),
            "configuration": config_path.name,
            "phonopy_version": version("phonopy"),
            "cat_tooling_version": __version__,
        },
        "artifacts": [
            {
                "name": "phonon_mode_visualization",
                "files": [path.name for path in outputs],
                "mode": {
                    "qpoint": mode.qpoint.tolist(),
                    "band": mode.band,
                    "frequency_THz": mode.frequency_THz,
                },
                "scene": asdict(options),
            }
        ],
    }
    target.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    return target


def run_export(
    filename: str | Path,
    *,
    qpoint=None,
    band: int | None = None,
    output_override: str | Path | None = None,
    amplitude_A: float | None = None,
    repeat: tuple[int, int, int] | None = None,
) -> list[Path]:
    """Render the configured exact mode without creating duplicate scientific tables."""
    config_path, config, files, phonon = _load(filename)
    mode = _mode(phonon, files.name, config, qpoint, band)
    options = _scene(mode, config, amplitude_A=amplitude_A, repeat=repeat)
    scene = build_scene(mode, options)
    export = dict(config.get("export", {}))
    directory_value = output_override or export.pop("directory", "phonovis-export")
    directory = Path(directory_value).expanduser()
    directory = (
        directory.resolve()
        if directory.is_absolute()
        else (config_path.parent / directory).resolve()
    )
    stem = str(export.pop("stem", _stem(mode)))
    formats = tuple(
        str(value).lower().lstrip(".") for value in export.pop("formats", ["pdf", "svg"])
    )
    unsupported = set(formats) - {"pdf", "svg", "png"}
    if unsupported:
        raise ValueError(f"Unsupported static formats: {sorted(unsupported)}")

    targets = [directory / f"{stem}.{suffix}" for suffix in formats]
    static_backend = str(export.pop("backend", "matplotlib"))
    transparent = bool(export.pop("transparent", True))
    width_in = float(export.pop("width_in", 6.0))
    dpi = int(export.pop("dpi", 300))
    window_size = tuple(int(value) for value in export.pop("window_size", [1600, 1600]))
    scale = int(export.pop("scale", 1))
    raster = bool(export.pop("raster", False))
    if static_backend == "matplotlib":
        written = render_figure(
            mode,
            scene,
            options,
            targets,
            width_in=width_in,
            dpi=dpi,
            transparent=transparent,
        )
    elif static_backend == "pyvista":
        from .pyvista_view import render_figure as render_pyvista_figure

        written = render_pyvista_figure(
            mode,
            scene,
            options,
            targets,
            window_size=window_size,
            transparent=transparent,
            scale=scale,
            raster=raster,
        )
    else:
        raise ValueError("Static backend must be 'matplotlib' or 'pyvista'")
    gif_enabled, gif = _enabled(export.pop("gif", False))
    if gif_enabled:
        backend = str(gif.pop("backend", "matplotlib"))
        target = directory / f"{stem}.gif"
        frames = int(gif.pop("frames", 36))
        fps = int(gif.pop("fps", 12))
        arrows = bool(gif.pop("arrows", False))
        gif_dpi = int(gif.pop("dpi", 160))
        if gif:
            raise TypeError(f"Unknown GIF options: {sorted(gif)}")
        if backend == "matplotlib":
            written.append(
                render_gif(
                    mode,
                    scene,
                    options,
                    target,
                    frames=frames,
                    fps=fps,
                    arrows=arrows,
                    dpi=gif_dpi,
                )
            )
        elif backend == "pyvista":
            from .pyvista_view import render_gif as render_pyvista_gif

            written.append(render_pyvista_gif(mode, scene, options, target, frames=frames, fps=fps))
        else:
            raise ValueError("GIF backend must be 'matplotlib' or 'pyvista'")

    vmd_enabled, vmd = _enabled(export.pop("vmd", False))
    if vmd_enabled:
        written.extend(write_vmd(mode, scene, options, directory / stem, **vmd))
    if export:
        raise TypeError(f"Unknown export options: {sorted(export)}")
    written.append(_write_manifest(directory, config_path, files, mode, options, written))
    return written


def run_explore(
    filename: str | Path,
    *,
    qpoint=None,
    band: int | None = None,
    viewer: str | None = None,
) -> None:
    """Open the clickable dashboard or one exact mode when a q-point is supplied."""
    config_path, config, files, phonon = _load(filename)
    explore = dict(config.get("explore", {}))
    configured_kind = str(explore.pop("kind", "dashboard"))
    configured_viewer = str(explore.pop("viewer", "pyvista"))
    if qpoint is not None or configured_kind == "mode":
        mode = _mode(phonon, files.name, config, qpoint, band)
        options = _scene(mode, config)
        scene = build_scene(mode, options)
        selected_viewer = viewer or configured_viewer
        if selected_viewer == "pyvista":
            from .pyvista_view import show_mode

            show_mode(mode, scene, options)
        elif selected_viewer == "matplotlib":
            show_figure(mode, scene, options)
        else:
            raise ValueError("Viewer must be 'matplotlib' or 'pyvista'")
        return
    if configured_kind != "dashboard":
        raise ValueError("Explore kind must be 'dashboard' or 'mode'")

    dispersion = dict(explore.pop("dispersion", {}))
    dos = dict(explore.pop("dos", {}))
    initial_q = int(explore.pop("initial_q_index", 0))
    initial_band = int(config.get("mode", {}).get("band", 0) if band is None else band)
    data = prepare_data(
        phonon,
        system=files.name,
        dispersion=dispersion,
        dos=dos,
        base=config_path.parent,
    )
    initial_mode = mode_from_band(data.band, phonon, initial_q, initial_band, system=files.name)
    options = _scene(initial_mode, config)
    show_explorer(
        phonon,
        data,
        options,
        system=files.name,
        initial_q_index=initial_q,
        initial_band=initial_band,
        **explore,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Explore phonon modes or recorded ionic trajectories")
    commands = parser.add_subparsers(dest="command", required=True)
    explore = commands.add_parser("explore", help="Clickable dispersion or explicit-q viewer")
    export = commands.add_parser("export", help="Static, GIF, and VMD outputs")
    for command in (explore, export):
        command.add_argument("config", type=Path, help="PhonoVis YAML workflow")
        command.add_argument("--q", nargs=3, type=float, metavar=("QX", "QY", "QZ"))
        command.add_argument("--band", type=int, help="Zero-based band index")
    explore.add_argument("--viewer", choices=("matplotlib", "pyvista"))
    export.add_argument("--out", type=Path, help="Override export directory")
    export.add_argument("--amplitude", type=float, help="Override displacement amplitude / A")
    export.add_argument("--repeat", nargs=3, type=int, metavar=("NX", "NY", "NZ"))
    trajectory = commands.add_parser("trajectory", help="View recorded ionic steps or save a GIF")
    trajectory.add_argument("input", type=Path, help="Extended XYZ trajectory")
    trajectory.add_argument("--gif", type=Path)
    trajectory.add_argument("--fps", type=int, default=6)
    args = parser.parse_args()
    if args.command == "trajectory":
        from .trajectory import run_trajectory
        if args.fps < 1:
            parser.error("--fps must be positive")
        run_trajectory(args.input, gif=args.gif, fps=args.fps)
        return
    if args.command == "explore":
        run_explore(args.config, qpoint=args.q, band=args.band, viewer=args.viewer)
    else:
        for path in run_export(
            args.config,
            qpoint=args.q,
            band=args.band,
            output_override=args.out,
            amplitude_A=args.amplitude,
            repeat=None if args.repeat is None else tuple(args.repeat),
        ):
            print(path)


if __name__ == "__main__":
    main()
