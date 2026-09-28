"""Generate independent, reproducible FreqBall mass-scan points."""

from __future__ import annotations

import hashlib
import shlex
from collections.abc import Sequence
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from cat_shared.artifacts import TableArtifact, write_artifacts
from cat_shared.mass_scan import _load


@dataclass(frozen=True)
class MassJob:
    index: int
    species: str
    mass_u: float
    directory: Path
    config: Path
    tensor: Path
    script: Path


def mass_values(
    *,
    values: Sequence[float] | None = None,
    interval: Sequence[float] | None = None,
    points: int = 10,
) -> tuple[float, ...]:
    """Return an explicit, positive and duplicate-free mass grid."""
    if (values is None) == (interval is None):
        raise ValueError("Select exactly one of explicit values or an interval")
    if interval is not None:
        if len(interval) != 2 or points < 2:
            raise ValueError("A mass interval needs two bounds and at least two points")
        start, stop = map(float, interval)
        selected = [start + index * (stop - start) / (points - 1) for index in range(points)]
    else:
        selected = [float(value) for value in values or ()]
    if not selected or any(value <= 0 for value in selected):
        raise ValueError("Masses must be non-empty and positive")
    if len(set(selected)) != len(selected):
        raise ValueError("Mass scan values must be unique")
    return tuple(selected)


def _absolute(base: Path, value: str | Path) -> str:
    path = Path(value).expanduser()
    return str(path.resolve() if path.is_absolute() else (base / path).resolve())


def _portable_config(base_config: Path, species: str, mass: float) -> dict[str, Any]:
    """Copy a workflow while anchoring its external inputs to the source YAML."""
    document = deepcopy(_load(base_config))
    base = base_config.parent
    system = document["system"]
    system["path"] = _absolute(base, system.get("path", "."))
    existing = system.get("masses", {})
    if existing and (not isinstance(existing, dict) or "primitive" in existing):
        raise ValueError("Mass scans require element-based system.masses overrides")
    system["masses"] = {**dict(existing), species: float(mass)}
    families = document.get("families")
    if isinstance(families, dict) and families.get("file"):
        families["file"] = _absolute(base, families["file"])
    modes = document.get("modes") or []
    modes = [modes] if isinstance(modes, dict) else modes
    for mode in modes:
        if isinstance(mode, dict) and mode.get("from_table"):
            mode["from_table"] = _absolute(base, mode["from_table"])
    document["output"] = {
        **dict(document.get("output", {})),
        "tensor": "freqball.h5",
        "directory": "analysis",
    }
    return document


def _token(value: float) -> str:
    return f"{value:.6f}".rstrip("0").rstrip(".").replace("-", "m").replace(".", "p")


def _render_script(template: str | None, job: MassJob, command: str) -> str:
    invocation = shlex.join([*shlex.split(command), str(job.config), "--action", "run"])
    if template is None:
        return f"#!/usr/bin/env bash\nset -euo pipefail\n{invocation}\n"
    replacements = {
        "FREQBALL_COMMAND": invocation,
        "CONFIG": shlex.quote(str(job.config)),
        "WORKDIR": shlex.quote(str(job.directory)),
        "POINT_DIR": shlex.quote(str(job.directory)),
        "OUTPUT": shlex.quote(str(job.tensor)),
        "SPECIES": shlex.quote(job.species),
        "MASS": f"{job.mass_u:.12g}",
        "INDEX": str(job.index),
        "JOB_NAME": f"FreqBall_{job.species}_{job.index:03d}",
    }
    result = template
    for name, value in replacements.items():
        result = result.replace("{{" + name + "}}", value)
    if "{{" in result or "}}" in result:
        raise ValueError("The job template contains an unknown {{TOKEN}}")
    if "{{FREQBALL_COMMAND}}" not in template:
        raise ValueError("The job template must contain {{FREQBALL_COMMAND}}")
    return result


def setup_mass_scan(
    base_config: str | Path,
    output: str | Path,
    *,
    species: str,
    masses: Sequence[float],
    job_template: str | Path | None = None,
    command: str = "cat-freqball",
    submit_command: str | None = None,
    force: bool = False,
) -> list[Path]:
    """Write one self-contained configuration and executable script per mass."""
    source = Path(base_config).expanduser().resolve()
    root = Path(output).expanduser().resolve()
    if not source.is_file() or not species or any(character.isspace() for character in species):
        raise ValueError("Supply an existing base configuration and one species token")
    selected = mass_values(values=masses)
    template = None
    if job_template is not None:
        template = Path(job_template).expanduser().resolve().read_text(encoding="utf-8")
    root.mkdir(parents=True, exist_ok=True)
    jobs: list[MassJob] = []
    for index, mass in enumerate(selected):
        folder = root / f"{index:03d}_{species}_{_token(mass)}u"
        config = folder / "freqball.yaml"
        script = folder / "run.sh"
        if not force and (config.exists() or script.exists()):
            raise FileExistsError(f"Mass point already exists: {folder}")
        folder.mkdir(parents=True, exist_ok=True)
        job = MassJob(index, species, mass, folder, config, folder / "freqball.h5", script)
        config.write_text(
            yaml.safe_dump(_portable_config(source, species, mass), sort_keys=False),
            encoding="utf-8",
        )
        script.write_text(_render_script(template, job, command), encoding="utf-8")
        script.chmod(script.stat().st_mode | 0o111)
        jobs.append(job)

    submit_script = None
    if submit_command:
        submit_script = root / "submit_all.sh"
        launcher = shlex.split(submit_command)
        if not launcher:
            raise ValueError("submit_command must contain an executable")
        lines = ["#!/usr/bin/env bash", "set -euo pipefail"]
        lines.extend(shlex.join([*launcher, str(job.script)]) for job in jobs)
        submit_script.write_text("\n".join(lines) + "\n", encoding="utf-8")
        submit_script.chmod(submit_script.stat().st_mode | 0o111)
    elif force and (root / "submit_all.sh").is_file():
        (root / "submit_all.sh").unlink()

    rows = pd.DataFrame(
        {
            "point_index": [job.index for job in jobs],
            "mass_species": species,
            "mass_u": [job.mass_u for job in jobs],
            "point_directory": [str(job.directory.relative_to(root)) for job in jobs],
            "config_file": [str(job.config.relative_to(root)) for job in jobs],
            "tensor_file": [str(job.tensor.relative_to(root)) for job in jobs],
            "job_file": [str(job.script.relative_to(root)) for job in jobs],
        }
    )
    artifact = TableArtifact(
        "freqball_mass_jobs",
        rows,
        ("mass_species", "point_index"),
        {"mass_u": "u"},
    )
    return [
        *(job.config for job in jobs),
        *(job.script for job in jobs),
        *((submit_script,) if submit_script is not None else ()),
        *write_artifacts(
            [artifact],
            root,
            run_metadata={
                "mass_species": species,
                "base_configuration": str(source),
                "base_configuration_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                "command": command,
                "submit_command": submit_command,
            },
        ),
    ]


def add_setup_arguments(setup) -> None:
    setup.add_argument("config", type=Path)
    setup.add_argument("-o", "--output", type=Path, required=True)
    setup.add_argument("--species", required=True)
    selection = setup.add_mutually_exclusive_group(required=True)
    selection.add_argument("--values", nargs="+", type=float)
    selection.add_argument("--range", nargs=2, type=float, dest="interval")
    setup.add_argument("--points", type=int, default=10)
    setup.add_argument("--job-template", type=Path)
    setup.add_argument("--command", default="cat-freqball")
    setup.add_argument(
        "--submit-command",
        help="Write submit_all.sh using this launcher, e.g. sbatch; never execute it",
    )
    setup.add_argument("--force", action="store_true")


def run_setup(args) -> list[Path]:
    selected = mass_values(values=args.values, interval=args.interval, points=args.points)
    return setup_mass_scan(
        args.config,
        args.output,
        species=args.species,
        masses=selected,
        job_template=args.job_template,
        command=args.command,
        submit_command=args.submit_command,
        force=args.force,
    )
