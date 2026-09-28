"""Order-independent comparison of relaxed periodic structures."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from cat_shared.artifacts import TableArtifact


@dataclass(frozen=True)
class ComparisonSettings:
    length_tolerance: float = 0.2
    site_tolerance: float = 0.3
    angle_tolerance: float = 5.0
    displacement_warning_A: float = 0.1
    lattice_warning_percent: float = 1.0
    angle_warning_deg: float = 1.0
    volume_warning_percent: float = 3.0


@dataclass(frozen=True)
class StructureComparison:
    name: str
    reference: object
    aligned: object | None
    source_indices: np.ndarray
    displacements_A: np.ndarray
    summary: TableArtifact
    sites: TableArtifact
    major_differences: tuple[str, ...]

    @property
    def matched(self) -> bool:
        return self.aligned is not None


def _space_group(structure, tolerance: float) -> tuple[str, int]:
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    analyzer = SpacegroupAnalyzer(structure, symprec=tolerance)
    return analyzer.get_space_group_symbol(), analyzer.get_space_group_number()


def _empty_sites() -> pd.DataFrame:
    return pd.DataFrame(
        columns=(
            "system",
            "comparison",
            "reference_site",
            "other_site",
            "reference_label",
            "other_label",
            "species",
            "dx_A",
            "dy_A",
            "dz_A",
            "displacement_A",
        )
    )


def compare_structures(
    reference,
    other,
    *,
    system: str,
    name: str = "reference-vs-other",
    settings: ComparisonSettings | None = None,
) -> StructureComparison:
    """Separate cell changes from mapped internal atomic displacements."""
    values = settings or ComparisonSettings()
    from pymatgen.analysis.structure_matcher import StructureMatcher

    differences = []
    same_composition = (
        reference.composition.reduced_composition == other.composition.reduced_composition
    )
    if not same_composition:
        differences.append("composition changed")
    if len(reference) != len(other):
        differences.append(f"site count changed ({len(reference)} -> {len(other)})")

    tagged = other.copy()
    property_name = "_cat_source_index"
    if property_name in tagged.site_properties:
        tagged.remove_site_property(property_name)
    tagged.add_site_property(property_name, list(range(len(tagged))))
    matcher = StructureMatcher(
        ltol=values.length_tolerance,
        stol=values.site_tolerance,
        angle_tol=values.angle_tolerance,
        primitive_cell=False,
        scale=False,
        attempt_supercell=False,
    )
    aligned = matcher.get_s2_like_s1(reference, tagged) if len(reference) == len(other) else None
    if aligned is None:
        differences.append("atomic sites could not be mapped")

    comparison_lattice = other.lattice if aligned is None else aligned.lattice
    ref_lengths = np.asarray(reference.lattice.abc)
    other_lengths = np.asarray(comparison_lattice.abc)
    length_change = 100 * (other_lengths / ref_lengths - 1)
    angle_change = np.asarray(comparison_lattice.angles) - np.asarray(reference.lattice.angles)
    volume_change = 100 * (other.volume / reference.volume - 1)
    if np.max(np.abs(length_change)) > values.lattice_warning_percent:
        differences.append(f"maximum lattice-length change is {np.max(np.abs(length_change)):.3g}%")
    if np.max(np.abs(angle_change)) > values.angle_warning_deg:
        differences.append(
            f"maximum lattice-angle change is {np.max(np.abs(angle_change)):.3g} degree"
        )
    if abs(volume_change) > values.volume_warning_percent:
        differences.append(f"volume change is {volume_change:.3g}%")

    reference_sg = _space_group(reference, 1e-3)
    other_sg = _space_group(other, 1e-3)
    if reference_sg[1] != other_sg[1]:
        differences.append(f"space group changed ({reference_sg[0]} -> {other_sg[0]})")

    source_indices = np.array([], dtype=int)
    displacement = np.empty((0, 3))
    site_rows = _empty_sites()
    if aligned is not None:
        source_indices = np.asarray(aligned.site_properties[property_name], dtype=int)
        delta_fractional = aligned.frac_coords - reference.frac_coords
        delta_fractional -= np.round(delta_fractional)
        displacement = delta_fractional @ reference.lattice.matrix
        norms = np.linalg.norm(displacement, axis=1)
        if float(norms.max(initial=0)) > values.displacement_warning_A:
            differences.append(f"maximum internal displacement is {norms.max():.4g} A")
        site_rows = pd.DataFrame(
            {
                "system": system,
                "comparison": name,
                "reference_site": np.arange(len(reference)),
                "other_site": source_indices,
                "reference_label": reference.labels,
                "other_label": aligned.labels,
                "species": [str(site.species) for site in reference],
                "dx_A": displacement[:, 0],
                "dy_A": displacement[:, 1],
                "dz_A": displacement[:, 2],
                "displacement_A": norms,
            }
        )
    else:
        norms = np.array([])

    row = {
        "system": system,
        "comparison": name,
        "status": "different" if differences else "similar",
        "mapped": aligned is not None,
        "reference_formula": reference.composition.reduced_formula,
        "other_formula": other.composition.reduced_formula,
        "reference_sites": len(reference),
        "other_sites": len(other),
        "reference_space_group": reference_sg[0],
        "other_space_group": other_sg[0],
        "a_change_percent": length_change[0],
        "b_change_percent": length_change[1],
        "c_change_percent": length_change[2],
        "alpha_change_deg": angle_change[0],
        "beta_change_deg": angle_change[1],
        "gamma_change_deg": angle_change[2],
        "volume_change_percent": volume_change,
        "rms_displacement_A": np.sqrt(np.mean(norms**2)) if len(norms) else np.nan,
        "mean_displacement_A": np.mean(norms) if len(norms) else np.nan,
        "max_displacement_A": np.max(norms) if len(norms) else np.nan,
        "major_differences": " | ".join(differences),
    }
    summary = TableArtifact(
        "structure_comparison",
        pd.DataFrame([row]),
        ("system", "comparison"),
        {
            "a_change_percent": "%",
            "b_change_percent": "%",
            "c_change_percent": "%",
            "alpha_change_deg": "degree",
            "beta_change_deg": "degree",
            "gamma_change_deg": "degree",
            "volume_change_percent": "%",
            "rms_displacement_A": "Angstrom",
            "mean_displacement_A": "Angstrom",
            "max_displacement_A": "Angstrom",
        },
        {
            "cell_and_internal_motion_are_reported_separately": True,
            "length_tolerance": values.length_tolerance,
            "site_tolerance": values.site_tolerance,
            "angle_tolerance": values.angle_tolerance,
        },
    )
    sites = TableArtifact(
        "site_displacements",
        site_rows,
        ("system", "comparison", "reference_site"),
        {
            "dx_A": "Angstrom",
            "dy_A": "Angstrom",
            "dz_A": "Angstrom",
            "displacement_A": "Angstrom",
        },
        {
            "site_indexing": "zero_based",
            "rigid_periodic_translation_removed_by_structure_matching": True,
        },
    )
    return StructureComparison(
        name,
        reference,
        aligned,
        source_indices,
        displacement,
        summary,
        sites,
        tuple(differences),
    )


def phonovis_displacement(result: StructureComparison):
    """Create PhonoVis's static-motion model; StructureTools never renders it."""
    if not result.matched:
        raise ValueError("An unmapped comparison cannot be visualized")
    from cat_tools.phonovis.model import StaticDisplacement

    motion = StaticDisplacement(
        f"{result.summary.data.iloc[0]['system']}: {result.name}",
        np.asarray(result.reference.lattice.matrix),
        np.asarray(result.reference.frac_coords),
        tuple(next(iter(site.species)).symbol for site in result.reference),
        np.asarray(result.displacements_A),
        label="mapped internal displacement",
    )
    motion.validate()
    return motion
