"""Export complete OUTCAR position/force blocks for relaxation visualization."""

from __future__ import annotations

import argparse
import csv
import re
import warnings
from itertools import groupby
from pathlib import Path

import numpy as np


def read_relaxation(outcar: Path, poscar: Path):
    """Keep recorded Cartesian coordinates and the cell current at each force block."""
    from ase import Atoms
    from ase.calculators.singlepoint import SinglePointCalculator
    from ase.io import read

    initial = read(poscar, format="vasp")
    symbols = initial.get_chemical_symbols()
    cell = initial.cell.array.copy()
    frames = []
    with outcar.open(errors="replace") as handle:
        lines = iter(handle)
        for line in lines:
            if "ions per type" in line:
                counts = [int(x) for x in line.split("=", 1)[1].split()]
                expected = [len(list(g)) for _, g in groupby(symbols)]
                if counts != expected:
                    raise ValueError("OUTCAR species counts do not match the supplied POSCAR")
            if "NIONS" in line:
                match = re.search(r"NIONS\s*=\s*(\d+)", line)
                if match and int(match[1]) != len(symbols):
                    raise ValueError("OUTCAR atom count does not match POSCAR")
            if "direct lattice vectors" in line:
                try:
                    cell = np.array([[float(x) for x in next(lines).split()[:3]] for _ in range(3)])
                    if cell.shape != (3, 3) or not np.isfinite(cell).all():
                        raise ValueError("Invalid lattice vectors")
                except StopIteration:
                    warnings.warn("Ignoring unfinished final lattice block", stacklevel=2)
                    break
            if "POSITION" not in line or "TOTAL-FORCE" not in line:
                continue
            rows = []
            try:
                next(lines)  # separator
                for _ in symbols:
                    raw = next(lines)
                    try:
                        row = [float(x.replace("D", "E")) for x in raw.split()]
                    except ValueError:
                        if not raw.endswith("\n"):
                            raise StopIteration from None
                        raise
                    if len(row) != 6 and not raw.endswith("\n"):
                        raise StopIteration
                    if len(row) != 6 or not np.isfinite(row).all():
                        raise ValueError("Malformed POSITION/TOTAL-FORCE block")
                    rows.append(row)
            except StopIteration:
                warnings.warn("Ignoring unfinished final position/force block", stacklevel=2)
                break
            values = np.asarray(rows)
            atoms = Atoms(symbols, positions=values[:, :3], cell=cell.copy(), pbc=True)
            atoms.calc = SinglePointCalculator(atoms, forces=values[:, 3:])
            atoms.info.update(ionic_step=len(frames) + 1, source=str(outcar))
            frames.append(atoms)
    if not frames:
        raise ValueError("No complete POSITION/TOTAL-FORCE blocks found")
    return frames


def discover_outcars(directory: Path) -> list[Path]:
    """Prefer a local archived/current chain; otherwise inspect numbered runs."""
    archives = sorted(
        (
            p
            for p in directory.glob("OUTCAR.s*")
            if p.is_file() and re.fullmatch(r"OUTCAR\.s\d+", p.name)
        ),
        key=lambda p: (int(p.name.split(".s")[-1]), p.name),
    )
    current = directory / "OUTCAR"
    if current.is_file() or archives:
        return archives + ([current] if current.is_file() else [])
    runs = sorted(
        (p for p in directory.iterdir() if p.is_dir() and re.fullmatch(r"relax_\d+", p.name)),
        key=lambda p: (int(p.name.split("_")[-1]), p.name),
    )
    files = []
    for run in runs:
        # Include archived segments within each numbered run too.
        if (run / "OUTCAR").is_file() or any(run.glob("OUTCAR.s[0-9]*")):
            files.extend(discover_outcars(run))
    if not files:
        raise ValueError(f"No OUTCAR, OUTCAR.sN or relax_N/OUTCAR found in {directory}")
    return files


def resolve_outcars(inputs) -> list[Path]:
    if isinstance(inputs, (str, Path)):
        inputs = [Path(inputs)]
    files = []
    for item in inputs:
        path = Path(item)
        files.extend(discover_outcars(path) if path.is_dir() else [path])

    def key(path):
        run = re.fullmatch(r"relax_(\d+)", path.parent.name)
        archive = re.fullmatch(r"OUTCAR\.s(\d+)", path.name)
        return (
            str(path.parent.parent if run else path.parent),
            int(run[1]) if run else -1,
            int(archive[1]) if archive else float("inf"),
            str(path),
        )

    files = sorted(dict.fromkeys(p.resolve() for p in files), key=key)
    if not files:
        raise ValueError("No OUTCAR inputs supplied")
    return files


def read_segments(outcars, poscar=None):
    frames = []
    symbols = None
    for segment, path in enumerate(outcars, 1):
        reference = poscar or path.parent / "POSCAR"
        try:
            current = read_relaxation(path, reference)
        except ValueError as exc:
            raise ValueError(f"{path}: {exc}") from exc
        current_symbols = current[0].get_chemical_symbols()
        if symbols is not None and current_symbols != symbols:
            raise ValueError(f"Species/order differs between relaxation segments: {path}")
        symbols = current_symbols
        for frame in current:
            frame.info.update(
                local_step=frame.info["ionic_step"],
                ionic_step=len(frames) + 1,
                segment=segment,
                segment_label=f"{path.parent.name}/{path.name}",
            )
            frames.append(frame)
    return frames


def export_relaxation(outcar: Path, poscar: Path, directory: Path, *, stride=1):
    from ase.io import write

    if stride < 1:
        raise ValueError("stride must be positive")
    outcars = resolve_outcars(outcar)
    all_frames = read_segments(outcars, poscar)
    indices = list(range(0, len(all_frames), stride))
    if indices[-1] != len(all_frames) - 1:
        indices.append(len(all_frames) - 1)
    frames = [all_frames[i] for i in indices]
    directory.mkdir(parents=True, exist_ok=True)
    write(directory / "relaxation.extxyz", frames, format="extxyz")
    write(directory / "relaxation.xyz", frames, format="xyz")
    # XYZ cannot carry lattice vectors. Supply per-frame lengths/angles for VMD.
    script = [
        "set script_dir [file dirname [info script]]",
        "mol new [file join $script_dir relaxation.xyz] type xyz waitfor all",
        "set molecule [molinfo top]",
        "set cells {",
    ]
    script += ["  {" + " ".join(f"{x:.12g}" for x in a.cell.cellpar()) + "}" for a in frames]
    script += [
        "}",
        "for {set i 0} {$i < [llength $cells]} {incr i} {",
        "  molinfo $molecule set frame $i",
        "  molinfo $molecule set {a b c alpha beta gamma} [lindex $cells $i]",
        "}",
        "mol delrep 0 top",
        "mol representation VDW 0.35 20",
        "mol color Element",
        "mol addrep top",
        "display projection Orthographic",
        "color Display Background white",
        "axes location Off",
        "animate style Once",
        "animate goto 0",
    ]
    (directory / "relaxation.vmd").write_text("\n".join(script) + "\n")
    with (directory / "steps.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["ionic_step", "max_force_eV_A", "volume_A3", "segment", "local_step", "source"]
        )
        for a in frames:
            writer.writerow(
                [
                    a.info["ionic_step"],
                    f"{np.linalg.norm(a.get_forces(), axis=1).max():.12g}",
                    f"{a.get_volume():.12g}",
                    a.info["segment"],
                    a.info["local_step"],
                    a.info["source"],
                ]
            )
    return frames


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "outcar",
        type=Path,
        nargs="*",
        help="Files or directories; default: discover in current directory",
    )
    parser.add_argument(
        "--poscar",
        type=Path,
        help="Shared species/order reference (default: POSCAR beside each OUTCAR)",
    )
    parser.add_argument("--out", type=Path, default=Path("relaxation-view"))
    parser.add_argument("--stride", type=int, default=1, help="Keep every Nth block and the last")
    parser.add_argument("--gif", action="store_true", help="Also render through PhonoVis")
    parser.add_argument("--fps", type=int, default=6)
    args = parser.parse_args(argv)
    if args.stride < 1 or args.fps < 1:
        parser.error("--stride and --fps must be positive")
    try:
        frames = export_relaxation(
            args.outcar or [Path.cwd()], args.poscar, args.out, stride=args.stride
        )
        if args.gif:
            from cat_tools.phonovis.trajectory import visualize

            visualize(frames, gif=args.out / "relaxation.gif", fps=args.fps)
    except (OSError, ValueError, ImportError) as exc:
        parser.exit(1, f"Error: {exc}\n")
    print(f"Exported {len(frames)} recorded geometries to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
