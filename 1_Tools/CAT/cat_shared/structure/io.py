"""Open structure I/O through Pymatgen rather than format-specific parsers."""

from __future__ import annotations

from pathlib import Path


def _structure_class():
    try:
        from pymatgen.core import Structure
    except ImportError as exc:
        raise ImportError(
            "Install StructureTools with: pip install 'cat-tooling[structure]'"
        ) from exc
    return Structure


def read_structure(path: str | Path, *, fmt: str | None = None):
    """Read any periodic structure understood by Pymatgen."""
    source = Path(path).expanduser().resolve()
    Structure = _structure_class()
    if fmt is None:
        return Structure.from_file(source)
    return Structure.from_str(source.read_text(encoding="utf-8"), fmt=fmt)


def write_structure(structure, path: str | Path, *, fmt: str | None = None) -> Path:
    """Write a Pymatgen Structure without constraining the output to VASP."""
    target = Path(path).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    if fmt is None:
        structure.to(filename=target)
    else:
        structure.to(filename=target, fmt=fmt)
    return target


def standardize_structure(
    structure,
    mode: str,
    *,
    symprec: float = 1e-3,
    angle_tolerance: float = 5.0,
):
    """Return an explicitly requested symmetry-standardized structure."""
    if mode == "none":
        return structure.copy()
    try:
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    except ImportError as exc:
        raise ImportError(
            "Install StructureTools with: pip install 'cat-tooling[structure]'"
        ) from exc
    analyzer = SpacegroupAnalyzer(
        structure,
        symprec=symprec,
        angle_tolerance=angle_tolerance,
    )
    methods = {
        "primitive": analyzer.get_primitive_standard_structure,
        "conventional": analyzer.get_conventional_standard_structure,
        "refined": analyzer.get_refined_structure,
    }
    if mode not in methods:
        raise ValueError("standardization must be none, primitive, conventional, or refined")
    return methods[mode]()
