# Exported from 3_Analysis/src/cat_tools/freqball/mass_workflow.py; update via sync_freqball_standalone.py.
# ruff: noqa -- generated copy; lint the shared implementation and exporter.
"""CLI for FreqBall artificial-mass setup and cross-mass analysis."""

from __future__ import annotations

import argparse
import hashlib
import multiprocessing
import subprocess
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from .mass_setup import add_setup_arguments, run_setup
from .artifacts import TableArtifact, write_artifacts
from .scan_paths import configured_families, tensor_paths

from . import __version__

from .mass_analysis import MassAnalysisOptions, load_mass_points, mass_scan_artifacts
from .viewer import show_ball


def _scan_species(root: Path, explicit: str | None) -> str:
    if explicit:
        return explicit
    manifest = yaml.safe_load((root / "manifest.yaml").read_text(encoding="utf-8")) or {}
    species = manifest.get("run", {}).get("mass_species")
    if not species:
        raise ValueError("Mass species is absent from the setup manifest; pass --species")
    return str(species)


def run_jobs(scan: str | Path, *, workers: int = 1) -> list[Path]:
    """Execute generated jobs locally with bounded concurrency and per-point logs."""
    if workers < 1:
        raise ValueError("workers must be positive")
    root = Path(scan).expanduser().resolve()
    jobs = pd.read_csv(root / "freqball_mass_jobs.csv")
    scripts = [(root / str(value)).resolve() for value in jobs["job_file"]]
    for script in scripts:
        if not script.is_relative_to(root) or not script.is_file():
            raise ValueError(f"Missing or external mass job: {script}")

    def execute(script):
        log = script.parent / "run.log"
        with log.open("w", encoding="utf-8") as handle:
            try:
                process = subprocess.run(
                    ["bash", str(script)],
                    cwd=script.parent,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
                code = process.returncode
            except OSError as exc:
                handle.write(str(exc))
                code = 127
        return {
            "job_file": str(script.relative_to(root)),
            "returncode": code,
            "log_file": str(log.relative_to(root)),
        }

    status_path = root / "freqball_mass_status.csv"
    # Mark every job unfinished before launching: an interrupted rerun must not
    # leave a previous success table authorizing collection of stale tensors.
    pd.DataFrame({"job_file": jobs["job_file"], "returncode": -1}).to_csv(status_path, index=False)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        statuses = list(pool.map(execute, scripts))
    pd.DataFrame(statuses).to_csv(status_path, index=False)
    failed = [row for row in statuses if row["returncode"]]
    if failed:
        raise RuntimeError(f"{len(failed)} mass jobs failed; see {status_path}")
    return [status_path, *(root / row["log_file"] for row in statuses)]


def _point_artifacts(point, *, species, options):
    return mass_scan_artifacts([point], species, options)


def run_collect(
    scan: str | Path,
    output: str | Path,
    *,
    species: str | None = None,
    families: str | Path | None = None,
    require_complete: bool = True,
    table_format: str = "csv",
    options: MassAnalysisOptions | None = None,
    workers: int = 1,
) -> list[Path]:
    if workers < 1:
        raise ValueError("workers must be positive")
    root = Path(scan).expanduser().resolve()
    status = root / "freqball_mass_status.csv"
    if status.is_file() and (pd.read_csv(status)["returncode"] != 0).any():
        raise ValueError("Scan has failed jobs; inspect run.log and rerun before collecting")
    species = _scan_species(root, species)
    family_path = families or configured_families(root)
    paths = tensor_paths(root, require_complete=require_complete)
    points = load_mass_points(paths, species, family_path)
    job_table = root / "freqball_mass_jobs.csv"
    if job_table.is_file():
        jobs = pd.read_csv(job_table)
        expected = {
            (root / row.tensor_file).resolve(): float(row.mass_u) for row in jobs.itertuples()
        }
        configs = {
            (root / row.tensor_file).resolve(): root / row.config_file for row in jobs.itertuples()
        }
        for point in points:
            saved_hash = point.result.metadata.get("configuration_sha256")
            if (
                saved_hash
                and saved_hash != hashlib.sha256(configs[point.path].read_bytes()).hexdigest()
            ):
                raise ValueError(f"{point.path}: point configuration changed since computation")
            if not np.isclose(point.mass_u, expected[point.path], rtol=0, atol=1e-8):
                raise ValueError(f"{point.path}: stored mass differs from the scan job table")
    if workers == 1:
        artifacts = mass_scan_artifacts(points, species, options)
    else:
        with ProcessPoolExecutor(
            max_workers=workers, mp_context=multiprocessing.get_context("spawn")
        ) as pool:
            groups = list(
                pool.map(partial(_point_artifacts, species=species, options=options), points)
            )
        artifacts = [
            TableArtifact(
                first.name,
                pd.concat([group[i].data for group in groups], ignore_index=True),
                first.keys,
                first.units,
                first.metadata,
            )
            for i, first in enumerate(groups[0])
        ]
    return write_artifacts(
        artifacts,
        output,
        table_format=table_format,
        run_metadata={
            "mass_species": species,
            "mass_points": len(points),
            "scan": str(root),
            "freqball_bundle_version": __version__,
        },
    )


def run_render(
    scan: str | Path,
    output: str | Path,
    *,
    species: str | None = None,
    families: str | Path | None = None,
    masses: list[float] | None = None,
    metric: str = "average_frequency_THz",
    directions: int = 1000,
    equivalence_group: int | None = None,
    low_cutoff: float = 2.0,
    inverse_floor: float = 0.1,
    deformation: float = 0.65,
    value_range: tuple[float, float] | None = None,
    cmap: str = "viridis",
) -> list[Path]:
    """Render selected scan points through PhonoVis's sole 3D backend."""
    root = Path(scan).expanduser().resolve()
    species = _scan_species(root, species)
    family_path = families or configured_families(root)
    points = load_mass_points(tensor_paths(root), species, family_path)
    if masses:
        requested = np.asarray(masses, dtype=float)
        points = tuple(
            point
            for point in points
            if np.isclose(point.mass_u, requested, rtol=0, atol=1e-8).any()
        )
        missing = [
            mass
            for mass in requested
            if not any(np.isclose(mass, point.mass_u, rtol=0, atol=1e-8) for point in points)
        ]
        if missing:
            raise ValueError(f"Requested masses are absent from the scan: {missing}")
    target = Path(output).expanduser().resolve()
    target.mkdir(parents=True, exist_ok=True)
    written = []
    for point in points:
        image = target / f"freqball_{species}_{point.mass_u:g}u_{metric}.png"
        show_ball(
            point.result,
            point.families,
            metric=metric,
            count=directions,
            equivalence_group=equivalence_group,
            low_cutoff=low_cutoff,
            inverse_floor=inverse_floor,
            deformation=deformation,
            value_range=value_range,
            cmap=cmap,
            screenshot=image,
            no_show=True,
        )
        written.append(image)
    return written


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Set up and analyze FreqBall mass scans")
    commands = parser.add_subparsers(dest="action", required=True)
    setup = commands.add_parser("setup", help="Generate one independent workflow per mass")
    add_setup_arguments(setup)

    run = commands.add_parser("run", help="Execute generated jobs locally with bounded concurrency")
    run.add_argument("scan", type=Path)
    run.add_argument("--workers", type=int, default=1)

    collect = commands.add_parser("collect", help="Create canonical cross-mass tables")
    collect.add_argument("scan", type=Path)
    collect.add_argument("-o", "--output", type=Path)
    collect.add_argument("--species")
    collect.add_argument("--families", type=Path)
    collect.add_argument("--allow-missing", action="store_true")
    collect.add_argument("--format", choices=("csv", "parquet"), default="csv")
    collect.add_argument("--workers", type=int, default=1)
    collect.add_argument("--directions", type=int, default=1000)
    collect.add_argument("--conduction-threshold", type=float, default=0.95)
    collect.add_argument("--tetrahedral-threshold", type=float, default=0.95)
    collect.add_argument("--low-cutoff", type=float, default=2.0)
    collect.add_argument("--inverse-floor", type=float, default=0.1)

    render = commands.add_parser(
        "render", help="Render comparable 3D balls with the PhonoVis backend"
    )
    render.add_argument("scan", type=Path)
    render.add_argument("-o", "--output", type=Path, required=True)
    render.add_argument("--species")
    render.add_argument("--families", type=Path)
    render.add_argument("--masses", nargs="+", type=float)
    render.add_argument("--metric", default="average_frequency_THz")
    render.add_argument("--directions", type=int, default=1000)
    render.add_argument("--equivalence-group", type=int)
    render.add_argument("--low-cutoff", type=float, default=2.0)
    render.add_argument("--inverse-floor", type=float, default=0.1)
    render.add_argument("--deformation", type=float, default=0.65)
    render.add_argument("--range", nargs=2, type=float, dest="value_range")
    render.add_argument("--cmap", default="viridis")
    return parser


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    from .parallel import context

    if context()[2] > 1:
        raise ValueError("Run mass-scan commands once; use --workers for local concurrency")
    if args.action == "setup":
        written = run_setup(args)
    elif args.action == "run":
        written = run_jobs(args.scan, workers=args.workers)
    elif args.action == "collect":
        options = MassAnalysisOptions(
            directions=args.directions,
            conduction_threshold=args.conduction_threshold,
            tetrahedral_threshold=args.tetrahedral_threshold,
            low_frequency_cutoff_THz=args.low_cutoff,
            inverse_frequency_floor_THz=args.inverse_floor,
        )
        written = run_collect(
            args.scan,
            args.output or args.scan / "analysis",
            species=args.species,
            families=args.families,
            require_complete=not args.allow_missing,
            table_format=args.format,
            options=options,
            workers=args.workers,
        )
    else:
        written = run_render(
            args.scan,
            args.output,
            species=args.species,
            families=args.families,
            masses=args.masses,
            metric=args.metric,
            directions=args.directions,
            equivalence_group=args.equivalence_group,
            low_cutoff=args.low_cutoff,
            inverse_floor=args.inverse_floor,
            deformation=args.deformation,
            value_range=None if args.value_range is None else tuple(args.value_range),
            cmap=args.cmap,
        )
    for path in written:
        print(path)


if __name__ == "__main__":
    main()
