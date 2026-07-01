#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from pymatgen.core import Structure


# =========================
# Editable defaults
# =========================
DEFAULT_STRUCTURE = Path("POSCAR")
DEFAULT_OUTPUT = Path("cohpBetween_auto.txt")

DEFAULT_RMIN = 0.0
DEFAULT_RMAX = 3.5

# Leave empty for all element pairs.
# Examples: ["Na-S"], ["P-S", "Sb-Se"], ["S-S"]
DEFAULT_PAIRS: list[str] = []

DEFAULT_ORBITAL_WISE = False
# =========================


@dataclass(frozen=True)
class Contact:
    i: int                  # 0-based atom index in home cell
    j: int                  # 0-based atom index of partner atom
    cell: tuple[int, int, int]
    distance: float
    el_i: str
    el_j: str


def parse_pair(text: str) -> tuple[str, str]:
    """
    Parse element pair strings like:
      Na-S
      Na:S
      Na,S
    """
    for sep in ("-", ":", ","):
        if sep in text:
            a, b = text.split(sep, 1)
            return a.strip(), b.strip()

    raise ValueError(
        f"Could not parse pair '{text}'. Use e.g. Na-S, P-S, S-S."
    )


def pair_matches(el_i: str, el_j: str, allowed_pairs: set[frozenset[str]]) -> bool:
    if not allowed_pairs:
        return True

    return frozenset((el_i, el_j)) in allowed_pairs


def canonical_contact_key(
    i: int,
    j: int,
    cell: tuple[int, int, int],
) -> tuple[int, int, tuple[int, int, int]]:
    """
    Avoid double-counting.

    Pair i -> j in cell (a,b,c) is equivalent to
    pair j -> i in cell (-a,-b,-c).
    """
    reverse_cell = (-cell[0], -cell[1], -cell[2])

    forward = (i, j, cell)
    reverse = (j, i, reverse_cell)

    return min(forward, reverse)


def format_cohpbetween_line(contact: Contact, orbital_wise: bool) -> str:
    i_lobster = contact.i + 1
    j_lobster = contact.j + 1
    a, b, c = contact.cell

    line = f"cohpBetween atom {i_lobster} atom {j_lobster} cell {a} {b} {c}"

    if orbital_wise:
        line += " orbitalWise"

    return line


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate explicit LOBSTER cohpBetween lines from a POSCAR/CIF/etc. "
            "using a distance cutoff."
        )
    )

    parser.add_argument(
        "structure",
        nargs="?",
        type=Path,
        default=DEFAULT_STRUCTURE,
        help="POSCAR / CONTCAR / CIF / etc. Default: POSCAR",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Output file with cohpBetween lines.",
    )
    parser.add_argument(
        "--pairs",
        nargs="*",
        default=DEFAULT_PAIRS,
        help="Element pairs to include, e.g. --pairs Na-S P-S S-S. Default: all pairs.",
    )
    parser.add_argument(
        "--rmin",
        type=float,
        default=DEFAULT_RMIN,
        help="Minimum distance in Å.",
    )
    parser.add_argument(
        "--rmax",
        type=float,
        default=DEFAULT_RMAX,
        help="Maximum distance in Å.",
    )
    parser.add_argument(
        "--orbital-wise",
        action="store_true",
        default=DEFAULT_ORBITAL_WISE,
        help="Append orbitalWise to every cohpBetween line.",
    )
    parser.add_argument(
        "--include-self-images",
        action="store_true",
        help=(
            "Allow interactions of an atom with its own periodic images. "
            "Usually not needed."
        ),
    )
    parser.add_argument(
        "--no-comments",
        action="store_true",
        help="Write only raw cohpBetween lines without comment lines.",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-8,
        help="Numerical tolerance for distance and image handling.",
    )

    args = parser.parse_args()

    structure = Structure.from_file(str(args.structure))

    allowed_pairs = {frozenset(parse_pair(p)) for p in args.pairs}

    center_indices, point_indices, offset_vectors, distances = structure.get_neighbor_list(
        r=args.rmax,
        numerical_tol=args.tol,
        exclude_self=True,
    )

    contacts_by_key: dict[tuple[int, int, tuple[int, int, int]], Contact] = {}

    for i_raw, j_raw, offset_raw, dist_raw in zip(
        center_indices,
        point_indices,
        offset_vectors,
        distances,
    ):
        i = int(i_raw)
        j = int(j_raw)
        dist = float(dist_raw)

        if dist < args.rmin - args.tol:
            continue

        cell = tuple(int(x) for x in np.rint(offset_raw))

        # Skip direct self-pairs and, unless requested, self-image contacts.
        if i == j:
            if cell == (0, 0, 0):
                continue
            if not args.include_self_images:
                continue

        el_i = structure[i].specie.symbol
        el_j = structure[j].specie.symbol

        if not pair_matches(el_i, el_j, allowed_pairs):
            continue

        key = canonical_contact_key(i, j, cell)
        ci, cj, ccell = key

        contact = Contact(
            i=ci,
            j=cj,
            cell=ccell,
            distance=dist,
            el_i=structure[ci].specie.symbol,
            el_j=structure[cj].specie.symbol,
        )

        # Keep shortest distance if numerical duplicates appear.
        if key not in contacts_by_key or dist < contacts_by_key[key].distance:
            contacts_by_key[key] = contact

    contacts = sorted(
        contacts_by_key.values(),
        key=lambda x: (
            min(x.el_i, x.el_j),
            max(x.el_i, x.el_j),
            x.distance,
            x.i,
            x.j,
            x.cell,
        ),
    )

    lines: list[str] = []

    if not args.no_comments:
        pair_text = ", ".join(args.pairs) if args.pairs else "all"
        lines.extend(
            [
                "! Automatically generated cohpBetween lines",
                f"! Structure : {args.structure}",
                f"! Formula   : {structure.composition.reduced_formula}",
                f"! Distance  : {args.rmin:.4f} <= d <= {args.rmax:.4f} Å",
                f"! Pairs     : {pair_text}",
                f"! Count     : {len(contacts)}",
                "! Atom indices are 1-based, matching POSCAR / LOBSTER convention.",
                "",
            ]
        )

    for contact in contacts:
        if not args.no_comments:
            lines.append(
                f"! {contact.el_i}{contact.i + 1} - "
                f"{contact.el_j}{contact.j + 1} "
                f"cell {contact.cell[0]} {contact.cell[1]} {contact.cell[2]} "
                f"d = {contact.distance:.6f} Å"
            )

        lines.append(format_cohpbetween_line(contact, args.orbital_wise))

    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print()
    print(f"Structure : {args.structure}")
    print(f"Formula   : {structure.composition.reduced_formula}")
    print(f"Pairs     : {', '.join(args.pairs) if args.pairs else 'all'}")
    print(f"Distance  : {args.rmin:.4f} <= d <= {args.rmax:.4f} Å")
    print(f"Contacts  : {len(contacts)}")
    print(f"Wrote     : {args.output}")
    print()


if __name__ == "__main__":
    main()
