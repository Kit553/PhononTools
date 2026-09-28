"""Compact structure summaries and symmetry-equivalent site tables."""

from __future__ import annotations
from collections.abc import Sequence

import numpy as np
import pandas as pd

from cat_shared.artifacts import TableArtifact
from cat_shared.structure.model import StructureFrame



def _analyzer(structure, symprec, angle_tolerance):
    try:
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    except ImportError as exc:
        raise ImportError(
            "Install StructureTools with: pip install 'cat-tooling[structure]'"
        ) from exc
    return SpacegroupAnalyzer(structure, symprec=symprec, angle_tolerance=angle_tolerance)


def summary_artifact(
    frames: list[StructureFrame], *, symprec: float = 1e-3, angle_tolerance: float = 5.0
) -> TableArtifact:
    rows = []
    for item in frames:
        structure = item.structure
        analyzer = _analyzer(structure, symprec, angle_tolerance)
        a, b, c = structure.lattice.abc
        alpha, beta, gamma = structure.lattice.angles
        rows.append(
            {
                "system": item.system,
                "frame": item.frame,
                "formula": structure.composition.formula,
                "reduced_formula": structure.composition.reduced_formula,
                "site_count": len(structure),
                "atom_count": float(structure.composition.num_atoms),
                "a_A": a,
                "b_A": b,
                "c_A": c,
                "alpha_deg": alpha,
                "beta_deg": beta,
                "gamma_deg": gamma,
                "volume_A3": structure.volume,
                "density_g_cm3": float(structure.density),
                "atom_density_m3": float(structure.composition.num_atoms / structure.volume * 1e30),
                "space_group": analyzer.get_space_group_symbol(),
                "space_group_number": analyzer.get_space_group_number(),
                "crystal_system": analyzer.get_crystal_system(),
            }
        )
    return TableArtifact(
        "structure_summary",
        pd.DataFrame(rows),
        ("system", "frame"),
        {
            "a_A": "Angstrom",
            "b_A": "Angstrom",
            "c_A": "Angstrom",
            "alpha_deg": "degree",
            "beta_deg": "degree",
            "gamma_deg": "degree",
            "volume_A3": "Angstrom^3",
            "density_g_cm3": "g/cm^3",
            "atom_density_m3": "atoms/m^3",
        },
        {"symprec": symprec, "angle_tolerance": angle_tolerance},
    )

def symmetry_artifact(
    frames: list[StructureFrame], *, symprec: float = 1e-3, angle_tolerance: float = 5.0
) -> TableArtifact:
    rows = []
    for item in frames:
        symmetrized = _analyzer(
            item.structure, symprec, angle_tolerance
        ).get_symmetrized_structure()
        for group, (indices, wyckoff) in enumerate(
            zip(symmetrized.equivalent_indices, symmetrized.wyckoff_symbols, strict=True)
        ):
            for site in indices:
                original = item.structure[site]
                rows.append(
                    {
                        "system": item.system,
                        "frame": item.frame,
                        "site": int(site),
                        "label": original.label,
                        "species": str(original.species),
                        "equivalent_group": group,
                        "wyckoff": wyckoff,
                        "multiplicity": len(indices),
                    }
                )
    return TableArtifact(
        "symmetry_sites",
        pd.DataFrame(rows),
        ("system", "frame", "site"),
        {},
        {
            "site_indexing": "zero_based",
            "symprec": symprec,
            "angle_tolerance": angle_tolerance,
        },
    )

def symmetry_scan_artifact(
    frames: list[StructureFrame],
    symprec_values: Sequence[float],
    *,
    angle_tolerance: float = 5.0,
) -> TableArtifact:
    """Detect symmetry across a distance-tolerance ladder."""

    tolerances = sorted({float(value) for value in symprec_values})
    if not tolerances or tolerances[0] <= 0:
        raise ValueError("Symmetry tolerances must be positive")

    rows = []

    for item in frames:
        for symprec in tolerances:
            analyzer = _analyzer(
                item.structure,
                symprec=symprec,
                angle_tolerance=angle_tolerance,
            )
            dataset = analyzer.get_symmetry_dataset()
            primitive = analyzer.find_primitive()
            refined = analyzer.get_refined_structure()

            rows.append(
                {
                    "system": item.system,
                    "frame": item.frame,
                    "symprec_A": symprec,
                    "angle_tolerance_deg": angle_tolerance,
                    "space_group": analyzer.get_space_group_symbol(),
                    "space_group_number": analyzer.get_space_group_number(),
                    "crystal_system": analyzer.get_crystal_system(),
                    "symmetry_operation_count": len(
                        analyzer.get_symmetry_operations()
                    ),
                    "equivalent_site_groups": len(
                        np.unique(dataset.equivalent_atoms)
                    ),
                    "primitive_site_count": (
                        len(primitive)
                        if primitive is not None
                        else len(item.structure)
                    ),
                    "refined_site_count": len(refined),
                }
            )

    table = pd.DataFrame(rows).sort_values(
        ["system", "frame", "symprec_A"]
    )

    table["first_for_space_group"] = ~table.duplicated(
        ["system", "frame", "space_group_number"]
    )

    return TableArtifact(
        "symmetry_scan",
        table,
        ("system", "frame", "symprec_A"),
        {
            "symprec_A": "Angstrom",
            "angle_tolerance_deg": "degree",
        },
        {
            "interpretation": (
                "The first occurrence of a space group gives its minimum "
                "tested detection tolerance; symprec is not an atomic "
                "displacement measurement."
            )
        },
    )