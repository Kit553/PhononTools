"""Publishable crystallographic LaTeX tables from Pymatgen structures."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class LatexTableOptions:
    caption: str | None = None
    label: str | None = None
    decimals: int = 6
    unique_sites: bool = True
    symprec: float = 1e-3
    angle_tolerance: float = 5.0


def _escape(value: str) -> str:
    return str(value).replace("\\", r"\textbackslash{}").replace("_", r"\_")


def _formula(value: str) -> str:
    pieces = re.findall(r"([A-Z][a-z]?)([0-9.]*)", value)
    return "".join(
        element + (rf"$_{{{count}}}$" if count and count != "1" else "")
        for element, count in pieces
    )


def latex_table(
    structure, *, system: str | None = None, options: LatexTableOptions | None = None
) -> str:
    """Render a booktabs-compatible table of symmetry-unique atomic coordinates."""
    values = options or LatexTableOptions()
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    analyzer = SpacegroupAnalyzer(
        structure,
        symprec=values.symprec,
        angle_tolerance=values.angle_tolerance,
    )
    symmetrized = analyzer.get_symmetrized_structure()
    groups = list(zip(symmetrized.equivalent_indices, symmetrized.wyckoff_symbols, strict=True))
    if not values.unique_sites:
        wyckoff_by_site = {site: wyckoff for indices, wyckoff in groups for site in indices}
        groups = [([site], wyckoff_by_site[site]) for site in range(len(structure))]
    formula = structure.composition.reduced_formula
    caption = values.caption or f"Crystallographic parameters of {_formula(formula)}."
    label = (
        values.label
        or f"tab:structure-{re.sub(r'[^A-Za-z0-9]+', '-', system or formula).strip('-').lower()}"
    )
    a, b, c = structure.lattice.abc
    alpha, beta, gamma = structure.lattice.angles
    lines = [
        r"\begin{table}",
        r"  \centering",
        f"  \\caption{{{caption}}}",
        f"  \\label{{{_escape(label)}}}",
        r"  \begin{tabular}{llrrrrr}",
        r"    \toprule",
        r"    Atom & Wyckoff & $x$ & $y$ & $z$ & Occupancy & Multiplicity \\",
        r"    \midrule",
    ]
    for indices, wyckoff in groups:
        site = structure[int(indices[0])]
        species = "/".join(
            str(item) if abs(float(occupancy) - 1) < 1e-12 else f"{item}:{occupancy:g}"
            for item, occupancy in site.species.items()
        )
        occupancy = sum(float(item) for item in site.species.values())
        coordinates = site.frac_coords
        lines.append(
            f"    {_escape(site.label or species)} & {_escape(wyckoff)} & "
            f"{coordinates[0]:.{values.decimals}f} & {coordinates[1]:.{values.decimals}f} & "
            f"{coordinates[2]:.{values.decimals}f} & {occupancy:.6g} & {len(indices)} \\\\"
        )
    lines.extend(
        [
            r"    \bottomrule",
            r"  \end{tabular}",
            r"  \par\smallskip",
            (
                f"  \\footnotesize Space group: {_escape(analyzer.get_space_group_symbol())} "
                f"(No. {analyzer.get_space_group_number()}); "
                f"$a={a:.6g}$ \\AA, $b={b:.6g}$ \\AA, $c={c:.6g}$ \\AA; "
                f"$\\alpha={alpha:.6g}^\\circ$, $\\beta={beta:.6g}^\\circ$, "
                f"$\\gamma={gamma:.6g}^\\circ$."
            ),
            r"\end{table}",
            "",
        ]
    )
    return "\n".join(lines)


def write_latex_table(
    structure,
    path: str | Path,
    *,
    system: str | None = None,
    options: LatexTableOptions | None = None,
) -> Path:
    target = Path(path).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(latex_table(structure, system=system, options=options), encoding="utf-8")
    return target
