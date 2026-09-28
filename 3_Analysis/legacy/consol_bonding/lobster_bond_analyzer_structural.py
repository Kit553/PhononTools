#!/usr/bin/env python3
"""
System-agnostic LOBSTER COHP/COBI/COOP analysis tool.

Core public methods:
    1. parse_files()
    2. run_diagnosis()
    3. run_integrated_analysis()
    4. run_energy_resolved_analysis()
    5. run_lobsterpy_analysis()

Main design choice:
    The analyzer does not know what "Na-S", "P-S", "S-S", etc. mean.
    All chemical meaning is injected through LobsterAnalysisConfig.

Minimal dependencies:
    pip install numpy pandas

Optional:
    pip install lobsterpy
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Literal

import numpy as np
import pandas as pd

from lobster_topology import (
    TopologyConfig,
    TopologyResult,
    build_polyhedra,
    classify_contact,
    load_manifest,
    read_poscar,
    topology_config_from_manifest,
    validate_manifest_structure,
)


Population = Literal["COHP", "COBI", "COOP"]


# =============================================================================
# Configuration layer: edit this, not the analyzer
# =============================================================================


@dataclass
class EnergyConfig:
    """Energy and sign conventions for output curves.

    LOBSTER CAR files are normally shifted so that E_F = 0 eV.
    `plot_energy_min/max` only filter the output CSVs; they do not affect parsing.
    """

    plot_energy_min: float | None = None
    plot_energy_max: float | None = None
    populations: tuple[Population, ...] = ("COHP", "COBI", "COOP")

    # Plot convention:
    # COHP is usually plotted as -COHP so bonding is positive.
    plot_sign: dict[str, float] = field(
        default_factory=lambda: {"COHP": -1.0, "COBI": 1.0, "COOP": 1.0}
    )


@dataclass
class InteractionRule:
    """Chemical interaction classifier.

    A rule can match either exact element pairs or pairs of element groups.

    Examples
    --------
    Exact elements:
        InteractionRule(
            name="Li-O",
            element_pairs=[("Li", "O")]
        )

    Element groups:
        InteractionRule(
            name="Pn-X",
            group_pairs=[("framework_cation", "anion")]
        )

    Split anion-anion pairs into intra/inter based on atom index groups:
        InteractionRule(
            name="anion_anion",
            group_pairs=[("anion", "anion")],
            split_by_atom_index_groups=True,
            same_group_name="anion_anion_intra_polyhedron",
            different_group_name="anion_anion_inter_polyhedron",
            unknown_group_name="anion_anion_unclassified",
        )
    """

    name: str
    element_pairs: list[tuple[str, str]] = field(default_factory=list)
    group_pairs: list[tuple[str, str]] = field(default_factory=list)
    split_by_atom_index_groups: bool = False
    same_group_name: str | None = None
    different_group_name: str | None = None
    unknown_group_name: str | None = None


@dataclass
class LobsterAnalysisConfig:
    """All system-specific choices live here."""

    name: str = "generic_lobster_analysis"

    # Example: {"mobile_cation": ["Na"], "framework_cation": ["P", "Sb"], "anion": ["S", "Se"]}
    element_groups: dict[str, list[str]] = field(default_factory=dict)

    # Deprecated compatibility field. Structural classification must not depend
    # on hard-coded atom ordering; topology is reconstructed from each structure.
    atom_index_groups: dict[str, list[int]] = field(default_factory=dict)

    # Structure-aware topology settings. These may be supplied by the optional
    # generator sidecar manifest and are always applied to the current structure.
    topology: dict[str, Any] = field(default_factory=dict)

    interaction_rules: list[InteractionRule] = field(default_factory=list)
    energy: EnergyConfig = field(default_factory=EnergyConfig)

    # The analyzer can run on any system, but these labels are useful in outputs.
    run_label_from_folder: bool = True

    @property
    def element_to_groups(self) -> dict[str, set[str]]:
        mapping: dict[str, set[str]] = {}
        for group, elements in self.element_groups.items():
            for element in elements:
                mapping.setdefault(element, set()).add(group)
        return mapping

    @property
    def atom_to_index_group(self) -> dict[int, str]:
        out: dict[int, str] = {}
        for group, indices in self.atom_index_groups.items():
            for idx in indices:
                out[int(idx)] = group
        return out

    @classmethod
    def from_json(cls, path: str | Path) -> "LobsterAnalysisConfig":
        raw = json.loads(Path(path).read_text())
        energy_raw = raw.get("energy", {})
        energy = EnergyConfig(**energy_raw)

        rules = []
        for r in raw.get("interaction_rules", []):
            rules.append(
                InteractionRule(
                    name=r["name"],
                    element_pairs=[tuple(x) for x in r.get("element_pairs", [])],
                    group_pairs=[tuple(x) for x in r.get("group_pairs", [])],
                    split_by_atom_index_groups=r.get("split_by_atom_index_groups", False),
                    same_group_name=r.get("same_group_name"),
                    different_group_name=r.get("different_group_name"),
                    unknown_group_name=r.get("unknown_group_name"),
                )
            )

        return cls(
            name=raw.get("name", "generic_lobster_analysis"),
            element_groups=raw.get("element_groups", {}),
            atom_index_groups=raw.get("atom_index_groups", {}),
            topology=raw.get("topology", {}),
            interaction_rules=rules,
            energy=energy,
            run_label_from_folder=raw.get("run_label_from_folder", True),
        )

    def to_json(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(asdict(self), indent=2))


def default_na3ps4_like_config() -> LobsterAnalysisConfig:
    """Current Na3(P/Sb)(S/Se)4-style example.

    To switch systems, change this config or use a JSON config file.
    The analyzer itself should not need changes.
    """
    return LobsterAnalysisConfig(
        name="Na3MCh4",
        element_groups={
            "mobile_cation": ["Na"],
            "framework_cation": ["P", "As", "Sb"],
            "anion": ["S", "Se"],
        },
        atom_index_groups={},
        topology={
            "framework_cations": ["P", "As", "Sb"],
            "anions": ["S", "Se"],
            "mobile_cations": ["Na"],
            "coordination": 4,
            "center_ligand_search_A": 3.4,
            "image_range": 2,
            "distance_tolerance_A": 0.03,
            "framework_class": "Pn-X",
            "mobile_class": "Na-X",
            "intra_class": "X-X_intra",
            "inter_class": "X-X_inter",
            "unclassified_class": "X-X_unclassified",
        },
        interaction_rules=[
            InteractionRule(
                name="Pn-X",
                group_pairs=[("framework_cation", "anion")],
            ),
            InteractionRule(
                name="Na-X",
                group_pairs=[("mobile_cation", "anion")],
            ),
            InteractionRule(
                name="X-X",
                group_pairs=[("anion", "anion")],
                split_by_atom_index_groups=False,
                same_group_name="X-X_intra",
                different_group_name="X-X_inter",
                unknown_group_name="X-X_unclassified",
            ),
        ],
        energy=EnergyConfig(
            plot_energy_min=-20.0,
            plot_energy_max=15.0,
            populations=("COHP", "COBI", "COOP"),
        ),
    )


# =============================================================================
# Internal data containers
# =============================================================================


@dataclass
class LobsterRun:
    run_id: str
    folder: Path
    files: dict[str, Path | None]
    structure: dict[str, Any] = field(default_factory=dict)
    lobsterin: dict[str, Any] = field(default_factory=dict)
    lobsterout: dict[str, Any] = field(default_factory=dict)
    total_interactions: pd.DataFrame = field(default_factory=pd.DataFrame)
    orbital_interactions: pd.DataFrame = field(default_factory=pd.DataFrame)
    energy_labels: pd.DataFrame = field(default_factory=pd.DataFrame)
    energy_curves: dict[str, pd.DataFrame] = field(default_factory=dict)
    band_overlap_quality: pd.DataFrame = field(default_factory=pd.DataFrame)
    topology: TopologyResult | None = None
    manifest: dict[str, Any] = field(default_factory=dict)
    manifest_validation: dict[str, Any] = field(default_factory=dict)


# =============================================================================
# Analyzer
# =============================================================================


class LobsterBondAnalyzer:
    def __init__(
        self,
        root: str | Path,
        out_dir: str | Path,
        config: LobsterAnalysisConfig | None = None,
        recursive: bool = True,
        manifest_mode: str = "auto",
    ) -> None:
        self.root = Path(root).resolve()
        self.out_dir = Path(out_dir).resolve()
        self.config = config or default_na3ps4_like_config()
        self.recursive = recursive
        if manifest_mode not in {"auto", "ignore", "require"}:
            raise ValueError("manifest_mode must be auto, ignore, or require")
        self.manifest_mode = manifest_mode
        self.runs: list[LobsterRun] = []
        self._parsed = False

        self.out_dir.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------------------------------
    # 1. Parsing
    # -------------------------------------------------------------------------

    def parse_files(self, write_raw: bool = True) -> list[LobsterRun]:
        """Parse all detected LOBSTER runs below root."""
        self.runs = []

        raw_dir = self.out_dir / "00_parsed_raw"
        if write_raw:
            raw_dir.mkdir(parents=True, exist_ok=True)

        for folder in self._find_run_folders():
            files = self._detect_files(folder)
            run = LobsterRun(
                run_id=self._make_run_id(folder),
                folder=folder,
                files=files,
            )

            # Structure-aware topology is reconstructed for every run. The optional
            # sidecar manifest supplies intended chemistry/settings but never replaces
            # the current structure geometry.
            structure_path = files.get("structure")
            if structure_path is None:
                run.structure = {"file": "", "formula": "", "n_sites": math.nan, "parse_error": "missing structure"}
            else:
                try:
                    structure_model = read_poscar(structure_path)
                    run.structure = structure_model.summary()
                    manifest = {}
                    if self.manifest_mode != "ignore" and files.get("manifest"):
                        manifest = load_manifest(files["manifest"])
                    run.manifest = manifest
                    run.manifest_validation = validate_manifest_structure(manifest, structure_model)
                    if self.manifest_mode == "require":
                        if not manifest:
                            raise RuntimeError(f"{run.run_id}: manifest required but missing")
                        if run.manifest_validation.get("structure_hash_matches") is not True:
                            raise RuntimeError(f"{run.run_id}: required manifest does not match current structure")

                    fallback_cfg = self._topology_config_from_analysis_config()
                    if manifest and run.manifest_validation.get("structure_hash_matches") is True:
                        top_cfg = topology_config_from_manifest(manifest, fallback=fallback_cfg)
                    else:
                        top_cfg = fallback_cfg
                    run.topology = build_polyhedra(structure_model, top_cfg)
                except Exception as exc:
                    if self.manifest_mode == "require":
                        raise
                    run.structure = {"file": str(structure_path), "formula": "", "n_sites": math.nan, "parse_error": repr(exc)}
                    run.topology = None
                    run.manifest_validation = {"manifest_present": bool(files.get("manifest")), "structure_hash_matches": None, "message": repr(exc)}

            run.lobsterin = self._parse_lobsterin(files.get("lobsterin"))
            run.lobsterout = self._parse_lobsterout(
                files.get("lobsterout") or files.get("lobster_log")
            )
            run.band_overlap_quality = self._parse_band_overlaps(files.get("bandOverlaps"))

            total, orbital = self._parse_and_merge_integrated(run)
            run.total_interactions = total
            run.orbital_interactions = orbital

            labels = []
            for population in self.config.energy.populations:
                car_file = files.get(f"{population}CAR")
                parsed = self._parse_energy_file(car_file, population)
                if parsed is None:
                    continue
                label_df, curve_df = parsed
                label_df.insert(0, "run_id", run.run_id)
                run.energy_curves[population] = curve_df
                labels.append(label_df)

            run.energy_labels = pd.concat(labels, ignore_index=True) if labels else pd.DataFrame()
            self.runs.append(run)

            if write_raw:
                prefix = raw_dir / self._safe(run.run_id)
                if not run.total_interactions.empty:
                    run.total_interactions.to_csv(f"{prefix}__total_interactions.csv", index=False)
                if not run.orbital_interactions.empty:
                    run.orbital_interactions.to_csv(f"{prefix}__orbital_interactions.csv", index=False)
                if not run.energy_labels.empty:
                    run.energy_labels.to_csv(f"{prefix}__energy_labels.csv", index=False)
                if run.topology is not None:
                    pd.DataFrame([
                        {
                            "polyhedron_id": poly.id,
                            "center_atom": poly.center_atom,
                            "center_element": poly.center_element,
                            "ligand_atom": lig["atom"],
                            "ligand_element": lig["element"],
                            "cell_x": lig["cell"][0],
                            "cell_y": lig["cell"][1],
                            "cell_z": lig["cell"][2],
                            "distance_A": lig["distance_A"],
                        }
                        for poly in run.topology.polyhedra for lig in poly.ligands
                    ]).to_csv(f"{prefix}__polyhedra.csv", index=False)
                    pd.DataFrame([
                        {
                            "atom1": c.atom1, "atom2": c.atom2,
                            "tx": c.cell[0], "ty": c.cell[1], "tz": c.cell[2],
                            "distance_A": c.distance_A, "interaction_class": c.category,
                            "shared_polyhedron": ";".join(c.polyhedron_ids),
                        } for c in run.topology.intra_contacts
                    ]).to_csv(f"{prefix}__intra_xx_contacts.csv", index=False)
                    Path(f"{prefix}__topology_validation.json").write_text(json.dumps({
                        "manifest_validation": run.manifest_validation,
                        "topology_warnings": run.topology.warnings,
                        "topology_config": asdict(run.topology.config),
                    }, indent=2))

        self._parsed = True

        if write_raw:
            self.config.to_json(raw_dir / "analysis_config_used.json")
            self._write_parsed_summary(raw_dir / "parsed_summary.json")

        return self.runs

    # -------------------------------------------------------------------------
    # 2. Diagnosis
    # -------------------------------------------------------------------------

    def run_diagnosis(self) -> dict[str, Path]:
        """Diagnose files, settings, charge spilling, interactions, overlaps."""
        self._ensure_parsed()

        out = self.out_dir / "01_diagnosis"
        out.mkdir(parents=True, exist_ok=True)

        file_rows = []
        quality_rows = []
        inventory_rows = []
        overlap_frames = []

        file_keys = [
            "lobsterout", "lobster_log", "lobsterin", "structure",
            "ICOHPLIST", "ICOBILIST", "ICOOPLIST",
            "COHPCAR", "COBICAR", "COOPCAR",
            "DOSCAR", "CHARGE", "GROSSPOP", "bandOverlaps", "manifest",
        ]

        for run in self.runs:
            file_row = {"run_id": run.run_id, "folder": str(run.folder)}
            for key in file_keys:
                path = run.files.get(key)
                file_row[f"has_{key}"] = bool(path)
                file_row[f"path_{key}"] = str(path) if path else ""
            file_rows.append(file_row)

            q = {
                "run_id": run.run_id,
                "folder": str(run.folder),
                "formula": run.structure.get("formula", ""),
                "n_sites": run.structure.get("n_sites", math.nan),
                **self._flatten_prefixed(run.lobsterin, "lobsterin"),
                **self._flatten_prefixed(run.lobsterout, "lobsterout"),
                "n_total_interactions_parsed": len(run.total_interactions),
                "n_orbital_interactions_parsed": len(run.orbital_interactions),
                "energy_populations_parsed": ",".join(sorted(run.energy_curves)),
                "topology_available": run.topology is not None,
                "n_polyhedra": len(run.topology.polyhedra) if run.topology else 0,
                "topology_warnings": " | ".join(run.topology.warnings) if run.topology else "",
                "manifest_present": bool(run.manifest),
                "manifest_structure_hash_matches": run.manifest_validation.get("structure_hash_matches"),
                "n_topology_distance_mismatch": (
                    int((~run.total_interactions["topology_distance_ok"]).sum())
                    if (not run.total_interactions.empty and "topology_distance_ok" in run.total_interactions.columns and run.topology is not None)
                    else 0
                ),
            }
            q["diagnosis_flag"] = self._diagnosis_flag(q)
            quality_rows.append(q)

            if not run.total_interactions.empty:
                inv = (
                    run.total_interactions
                    .groupby(["element_pair", "interaction_class"], dropna=False)
                    .agg(
                        n=("interaction_no", "count"),
                        distance_min_A=("distance_A", "min"),
                        distance_mean_A=("distance_A", "mean"),
                        distance_max_A=("distance_A", "max"),
                        mean_minus_ICOHP=("minus_ICOHP", "mean"),
                        sum_minus_ICOHP=("minus_ICOHP", "sum"),
                        mean_ICOBI=("ICOBI", "mean"),
                        mean_ICOOP=("ICOOP", "mean"),
                    )
                    .reset_index()
                )
                inv.insert(0, "run_id", run.run_id)
                inventory_rows.append(inv)

            if not run.band_overlap_quality.empty:
                temp = run.band_overlap_quality.copy()
                temp.insert(0, "run_id", run.run_id)
                overlap_frames.append(temp)

        file_inventory = pd.DataFrame(file_rows)
        run_quality = pd.DataFrame(quality_rows)
        interaction_inventory = (
            pd.concat(inventory_rows, ignore_index=True) if inventory_rows else pd.DataFrame()
        )
        overlap_quality = (
            pd.concat(overlap_frames, ignore_index=True) if overlap_frames else pd.DataFrame()
        )

        paths = {
            "file_inventory": out / "file_inventory.csv",
            "run_quality": out / "run_quality.csv",
            "interaction_inventory": out / "interaction_inventory.csv",
            "band_overlap_quality": out / "band_overlap_quality.csv",
            "technical_summary": out / "technical_summary.txt",
        }

        file_inventory.to_csv(paths["file_inventory"], index=False)
        run_quality.to_csv(paths["run_quality"], index=False)
        interaction_inventory.to_csv(paths["interaction_inventory"], index=False)
        overlap_quality.to_csv(paths["band_overlap_quality"], index=False)
        paths["technical_summary"].write_text(
            self._technical_summary(run_quality, interaction_inventory, overlap_quality)
        )

        return paths

    # -------------------------------------------------------------------------
    # 3. Integrated analysis
    # -------------------------------------------------------------------------

    def run_integrated_analysis(self) -> dict[str, Path]:
        """Write scalar integrated COHP/COBI/COOP tables for bar plots."""
        self._ensure_parsed()

        out = self.out_dir / "02_integrated_analysis"
        out.mkdir(parents=True, exist_ok=True)

        total = self._concat("total_interactions")
        orbital = self._concat("orbital_interactions")

        if total.empty:
            raise RuntimeError("No integrated interactions found. Check ICOHPLIST/ICOBILIST/ICOOPLIST.")

        pair_stats = self._aggregate_integrated(total, ["run_id", "element_pair", "interaction_class"])
        class_stats = self._aggregate_integrated(total, ["run_id", "interaction_class"])

        bar_long = self._bar_long(pair_stats, level="element_pair")
        bar_long = pd.concat([bar_long, self._bar_long(class_stats, level="interaction_class")], ignore_index=True)

        orbital_stats = pd.DataFrame()
        if not orbital.empty:
            orbital_stats = (
                orbital
                .groupby(["run_id", "element_pair", "interaction_class", "orbital_pair"], dropna=False)
                .agg(
                    n=("interaction_no", "count"),
                    minus_ICOHP_sum=("minus_ICOHP", "sum"),
                    minus_ICOHP_mean=("minus_ICOHP", "mean"),
                    ICOBI_sum=("ICOBI", "sum"),
                    ICOBI_mean=("ICOBI", "mean"),
                    ICOOP_sum=("ICOOP", "sum"),
                    ICOOP_mean=("ICOOP", "mean"),
                )
                .reset_index()
                .sort_values(["run_id", "interaction_class", "minus_ICOHP_sum"], ascending=[True, True, False])
            )

        paths = {
            "total_interactions": out / "integrated_total_interactions.csv",
            "orbital_interactions": out / "integrated_orbital_interactions.csv",
            "bar_pair_stats": out / "bar_pair_stats.csv",
            "bar_class_stats": out / "bar_class_stats.csv",
            "plot_bar_integrated_long": out / "plot_bar_integrated_long.csv",
            "plot_bar_orbital_contributions": out / "plot_bar_orbital_contributions.csv",
        }

        total.to_csv(paths["total_interactions"], index=False)
        orbital.to_csv(paths["orbital_interactions"], index=False)
        pair_stats.to_csv(paths["bar_pair_stats"], index=False)
        class_stats.to_csv(paths["bar_class_stats"], index=False)
        bar_long.to_csv(paths["plot_bar_integrated_long"], index=False)
        orbital_stats.to_csv(paths["plot_bar_orbital_contributions"], index=False)

        return paths

    # -------------------------------------------------------------------------
    # 4. Energy-resolved analysis
    # -------------------------------------------------------------------------

    def run_energy_resolved_analysis(
        self,
        group_by: tuple[str, ...] = ("interaction_class", "element_pair"),
        include_individual_total_curves: bool = False,
    ) -> dict[str, Path]:
        """Average energy-resolved curves over interaction classes/pairs."""
        self._ensure_parsed()

        out = self.out_dir / "03_energy_resolved_analysis"
        out.mkdir(parents=True, exist_ok=True)

        mean_rows: list[dict[str, Any]] = []
        individual_rows: list[dict[str, Any]] = []
        all_label_frames = []

        for run in self.runs:
            if not run.energy_labels.empty:
                all_label_frames.append(run.energy_labels)

            if run.total_interactions.empty:
                continue

            interaction_info = run.total_interactions[
                [
                    "interaction_no",
                    "atom1", "atom2",
                    "element_pair",
                    "interaction_class",
                    "distance_A",
                    "distance_reconstructed_A", "distance_error_A",
                    "topology_source", "shared_polyhedron",
                    "ICOHP", "minus_ICOHP", "ICOBI", "ICOOP",
                ]
            ].drop_duplicates("interaction_no")

            for population, curves in run.energy_curves.items():
                labels = run.energy_labels
                if labels.empty:
                    continue

                labels = labels[
                    (labels["population"] == population)
                    & (labels["is_total_interaction"])
                ].merge(interaction_info, on="interaction_no", how="left")

                if labels.empty:
                    continue

                energy = curves["energy_eV_vs_EF"].to_numpy(dtype=float)
                mask = self._energy_mask(energy)
                energy = energy[mask]

                sign = float(self.config.energy.plot_sign.get(population, 1.0))

                for group_col in group_by:
                    if group_col not in labels.columns:
                        continue

                    for group_value, group in labels.dropna(subset=[group_col]).groupby(group_col):
                        raw_stack, int_stack = self._curve_stack(curves, population, group["label_index"].astype(int), mask)
                        if raw_stack.size == 0:
                            continue

                        mean_raw = np.nanmean(raw_stack, axis=0)
                        mean_int = np.nanmean(int_stack, axis=0)

                        for e, rv, iv in zip(energy, mean_raw, mean_int):
                            mean_rows.append(
                                {
                                    "run_id": run.run_id,
                                    "population": population,
                                    "grouping": group_col,
                                    "group": group_value,
                                    "n_interactions_averaged": len(group),
                                    "energy_eV_vs_EF": e,
                                    "mean_raw": rv,
                                    "mean_integrated_raw": iv,
                                    "plot_sign": sign,
                                    "mean_plot": sign * rv,
                                    "mean_integrated_plot": sign * iv,
                                }
                            )

                if include_individual_total_curves:
                    for _, label in labels.iterrows():
                        raw_stack, int_stack = self._curve_stack(curves, population, [int(label["label_index"])], mask)
                        if raw_stack.size == 0:
                            continue
                        for e, rv, iv in zip(energy, raw_stack[0], int_stack[0]):
                            individual_rows.append(
                                {
                                    "run_id": run.run_id,
                                    "population": population,
                                    "interaction_no": int(label["interaction_no"]),
                                    "atom1": label["atom1"],
                                    "atom2": label["atom2"],
                                    "element_pair": label["element_pair"],
                                    "interaction_class": label["interaction_class"],
                                    "distance_A": label["distance_A"],
                                    "energy_eV_vs_EF": e,
                                    "raw": rv,
                                    "integrated_raw": iv,
                                    "plot_sign": sign,
                                    "plot": sign * rv,
                                    "integrated_plot": sign * iv,
                                }
                            )

        labels_all = pd.concat(all_label_frames, ignore_index=True) if all_label_frames else pd.DataFrame()
        mean_long = pd.DataFrame(mean_rows)
        individual_long = pd.DataFrame(individual_rows)
        mean_wide = self._energy_wide(mean_long)

        paths = {
            "energy_curve_labels": out / "energy_curve_labels.csv",
            "energy_mean_curves_long": out / "energy_mean_curves_long.csv",
            "energy_mean_curves_wide": out / "energy_mean_curves_wide.csv",
            "energy_individual_total_curves_long": out / "energy_individual_total_curves_long.csv",
        }

        labels_all.to_csv(paths["energy_curve_labels"], index=False)
        mean_long.to_csv(paths["energy_mean_curves_long"], index=False)
        mean_wide.to_csv(paths["energy_mean_curves_wide"], index=False)
        individual_long.to_csv(paths["energy_individual_total_curves_long"], index=False)

        return paths

    # -------------------------------------------------------------------------
    # 5. LobsterPy wrapper
    # -------------------------------------------------------------------------

    def run_lobsterpy_analysis(
        self,
        population: Population = "COHP",
        all_bonds: bool = True,
        orbital_resolved: bool = False,
        cutoff_icohp: float = 0.10,
        noise_cutoff: float = 0.10,
    ) -> dict[str, Path]:
        """Run LobsterPy CLI, if installed."""
        self._ensure_parsed()

        out = self.out_dir / "04_lobsterpy_analysis"
        out.mkdir(parents=True, exist_ok=True)

        exe = shutil.which("lobsterpy")
        rows = []

        if exe is None:
            path = out / "lobsterpy_status.csv"
            pd.DataFrame(
                [{"status": "lobsterpy executable not found", "hint": "pip install lobsterpy"}]
            ).to_csv(path, index=False)
            return {"lobsterpy_status": path}

        for run in self.runs:
            car_key = f"{population}CAR"
            list_key = {"COHP": "ICOHPLIST", "COBI": "ICOBILIST", "COOP": "ICOOPLIST"}[population]

            structure = run.files.get("structure")
            car_file = run.files.get(car_key)
            list_file = run.files.get(list_key)
            charge = run.files.get("CHARGE")

            missing = [name for name, p in [("structure", structure), (car_key, car_file), (list_key, list_file)] if not p]

            prefix = out / self._safe(run.run_id)
            stdout = prefix.with_suffix(f".lobsterpy_{population.lower()}.txt")
            stderr = prefix.with_suffix(f".lobsterpy_{population.lower()}.stderr.txt")

            if missing:
                rows.append({"run_id": run.run_id, "status": "skipped_missing_files", "missing": ",".join(missing)})
                continue

            cmd = [
                exe, "description",
                "-fcohp", str(car_file),
                "-ficohp", str(list_file),
                "-fstruct", str(structure),
                "--cutofficohp", str(cutoff_icohp),
                "--noisecutoff", str(noise_cutoff),
            ]

            if charge:
                cmd.extend(["-fcharge", str(charge)])
            if all_bonds:
                cmd.append("-allb")
            if orbital_resolved:
                cmd.append("-orbresol")
            if population == "COBI":
                cmd.append("--cobis")
            elif population == "COOP":
                cmd.append("--coops")

            result = subprocess.run(cmd, text=True, capture_output=True, check=False)
            stdout.write_text(result.stdout or "")
            stderr.write_text(result.stderr or "")

            rows.append(
                {
                    "run_id": run.run_id,
                    "population": population,
                    "status": "ok" if result.returncode == 0 else f"failed_{result.returncode}",
                    "stdout": str(stdout),
                    "stderr": str(stderr),
                    "command": " ".join(cmd),
                }
            )

        path = out / f"lobsterpy_status_{population.lower()}.csv"
        pd.DataFrame(rows).to_csv(path, index=False)
        return {"lobsterpy_status": path}

    def run_all(
        self, include_lobsterpy: bool = False, include_individual_curves: bool = False
    ) -> dict[str, dict[str, Path]]:
        self.parse_files(write_raw=True)
        outputs = {
            "diagnosis": self.run_diagnosis(),
            "integrated": self.run_integrated_analysis(),
            "energy_resolved": self.run_energy_resolved_analysis(
                include_individual_total_curves=include_individual_curves
            ),
        }
        if include_lobsterpy:
            outputs["lobsterpy"] = self.run_lobsterpy_analysis()
        return outputs

    # =========================================================================
    # File discovery
    # =========================================================================

    def _find_run_folders(self) -> list[Path]:
        if not self.root.exists():
            raise FileNotFoundError(self.root)

        candidates = [self.root] if self._detect_files(self.root).get("ICOHPLIST") else []

        iterator = self.root.rglob("*") if self.recursive else self.root.glob("*")
        for p in iterator:
            if p.is_file() and p.name.upper().startswith("ICOHPLIST"):
                candidates.append(p.parent)

        return sorted(set(candidates))

    def _make_run_id(self, folder: Path) -> str:
        if folder == self.root:
            return folder.name
        return str(folder.relative_to(self.root))

    def _detect_files(self, folder: Path) -> dict[str, Path | None]:
        patterns = {
            "lobsterout": ["lobsterout*", "LOBSTEROUT*"],
            "lobster_log": ["lobster.log", "LOBSTER.LOG"],
            "lobsterin": ["lobsterin*", "LOBSTERIN*"],
            "structure": ["POSCAR.lobster*", "POSCAR*lobster*", "CONTCAR*", "POSCAR*"],
            "ICOHPLIST": ["ICOHPLIST*.lobster", "ICOHPLIST*"],
            "ICOBILIST": ["ICOBILIST*.lobster", "ICOBILIST*"],
            "ICOOPLIST": ["ICOOPLIST*.lobster", "ICOOPLIST*"],
            "COHPCAR": ["COHPCAR*.lobster", "COHPCAR*"],
            "COBICAR": ["COBICAR*.lobster", "COBICAR*"],
            "COOPCAR": ["COOPCAR*.lobster", "COOPCAR*"],
            "DOSCAR": ["DOSCAR*.lobster", "DOSCAR*"],
            "CHARGE": ["CHARGE*.lobster", "CHARGE*"],
            "GROSSPOP": ["GROSSPOP*.lobster", "GROSSPOP*"],
            "bandOverlaps": ["bandOverlaps*.lobster", "bandOverlaps*"],
            "manifest": ["lobsterin*.bonding.json", "bonding_manifest.json", "lobster_bonding_manifest.json"],
        }
        return {key: self._first(folder, pats) for key, pats in patterns.items()}

    @staticmethod
    def _first(folder: Path, patterns: Iterable[str]) -> Path | None:
        for pat in patterns:
            hits = sorted(p for p in folder.glob(pat) if p.is_file())
            if hits:
                return hits[0]
        return None

    # =========================================================================
    # Parsers
    # =========================================================================

    def _parse_poscar(self, path: Path | None) -> dict[str, Any]:
        if path is None:
            return {"file": "", "formula": "", "n_sites": math.nan, "parse_error": "missing structure"}
        try:
            return read_poscar(path).summary()
        except Exception as exc:
            return {"file": str(path), "formula": "", "n_sites": math.nan, "parse_error": repr(exc)}

    def _parse_lobsterin(self, path: Path | None) -> dict[str, Any]:
        text = self._read(path)
        out = {
            "file": str(path) if path else "",
            "COHPstartEnergy": math.nan,
            "COHPendEnergy": math.nan,
            "basisSet": "",
            "basisFunctions": [],
            "cohpGenerator": [],
            "cohpBetween": [],
            "orbitalWise_requested": False,
        }
        if not text:
            return out

        for raw in text.splitlines():
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            key = parts[0]

            if key == "COHPstartEnergy" and len(parts) > 1:
                out["COHPstartEnergy"] = self._float(parts[1])
            elif key == "COHPendEnergy" and len(parts) > 1:
                out["COHPendEnergy"] = self._float(parts[1])
            elif key == "basisSet" and len(parts) > 1:
                out["basisSet"] = " ".join(parts[1:])
            elif key == "basisFunctions":
                out["basisFunctions"].append(line)
            elif key == "cohpGenerator":
                out["cohpGenerator"].append(line)
            elif key == "cohpBetween":
                out["cohpBetween"].append(line)

            if "orbitalWise" in line:
                out["orbitalWise_requested"] = True

        out["n_cohpGenerator"] = len(out["cohpGenerator"])
        out["n_cohpBetween"] = len(out["cohpBetween"])
        out["basisFunctions"] = " ; ".join(out["basisFunctions"])
        out["cohpGenerator"] = " ; ".join(out["cohpGenerator"])
        out["cohpBetween"] = " ; ".join(out["cohpBetween"])
        return out

    def _parse_lobsterout(self, path: Path | None) -> dict[str, Any]:
        text = self._read(path)
        out = {
            "file": str(path) if path else "",
            "version": "",
            "finished_normally": False,
            "abs_charge_spilling_percent": math.nan,
            "total_charge_spilling_percent": math.nan,
            "electrons_recovered": math.nan,
            "electrons_total": math.nan,
            "interactions_found_reported": math.nan,
            "has_warning": False,
            "has_error": False,
            "warnings_notes_errors": "",
        }
        if not text:
            return out

        if m := re.search(r"LOBSTER\s+v([0-9A-Za-z_.+-]+)", text):
            out["version"] = m.group(1)
        out["finished_normally"] = bool(re.search(r"finished\s+in", text, re.I))

        if m := re.search(r"abs\.\s*charge\s+spilling:\s*([0-9.+\-EeDd]+)%", text, re.I):
            out["abs_charge_spilling_percent"] = self._float(m.group(1))
        if m := re.search(r"total\s+charge\s+spilling:\s*([0-9.+\-EeDd]+)%", text, re.I):
            out["total_charge_spilling_percent"] = self._float(m.group(1))
        if m := re.search(r"number\s+of\s+electrons\s+recovered\s+by\s+projection:\s*([0-9.+\-EeDd]+)\s+of\s+([0-9.+\-EeDd]+)", text, re.I):
            out["electrons_recovered"] = self._float(m.group(1))
            out["electrons_total"] = self._float(m.group(2))
        if m := re.search(r"found\s+(\d+)\s+interactions", text, re.I):
            out["interactions_found_reported"] = int(m.group(1))

        flagged = []
        for line in text.splitlines():
            if any(x in line.upper() for x in ["WARNING", "ERROR", "NOTE"]):
                flagged.append(line.strip())

        out["has_warning"] = any("WARNING" in x.upper() for x in flagged)
        out["has_error"] = any("ERROR" in x.upper() for x in flagged)
        out["warnings_notes_errors"] = " | ".join(flagged)
        return out

    def _parse_and_merge_integrated(self, run: LobsterRun) -> tuple[pd.DataFrame, pd.DataFrame]:
        raw: dict[str, pd.DataFrame] = {}
        for population, list_key, value_col in [
            ("COHP", "ICOHPLIST", "ICOHP"),
            ("COBI", "ICOBILIST", "ICOBI"),
            ("COOP", "ICOOPLIST", "ICOOP"),
        ]:
            df = self._parse_integrated_list(run.files.get(list_key), value_col=value_col)
            if not df.empty:
                raw[population] = df

        if not raw:
            return pd.DataFrame(), pd.DataFrame()

        base = raw.get("COHP", next(iter(raw.values())))

        total = base[base["is_total_interaction"]].copy()
        orbital = base[~base["is_total_interaction"]].copy()

        for population, value_col in [("COBI", "ICOBI"), ("COOP", "ICOOP")]:
            other = raw.get(population)
            if other is None:
                total[value_col] = math.nan
                orbital[value_col] = math.nan
                continue

            other_total = other[other["is_total_interaction"]][["interaction_no", value_col]]
            total = total.merge(other_total, on="interaction_no", how="left")

            other_orbital = other[~other["is_total_interaction"]][["interaction_no", "atom1", "atom2", value_col]]
            orbital = orbital.merge(other_orbital, on=["interaction_no", "atom1", "atom2"], how="left")

        for df in [total, orbital]:
            if df.empty:
                continue
            df.insert(0, "run_id", run.run_id)
            df.insert(1, "folder", str(run.folder))
            df["element_pair"] = df.apply(lambda r: self._element_pair(r["el1"], r["el2"]), axis=1)
            structural = df.apply(lambda r: self._classify_interaction(r, run.topology), axis=1, result_type="expand")
            for col in structural.columns:
                df[col] = structural[col].values
            df["distance_error_A"] = (df["distance_reconstructed_A"] - df["distance_A"]).abs()
            df["topology_distance_ok"] = df["distance_error_A"] <= (
                run.topology.config.distance_tolerance_A if run.topology else math.nan
            )
            df["minus_ICOHP"] = -df["ICOHP"] if "ICOHP" in df.columns else math.nan
            df["cell_translation"] = df[["tx", "ty", "tz"]].astype(str).agg(" ".join, axis=1)

        if not orbital.empty:
            orbital["orbital_pair"] = orbital.apply(
                lambda r: f"{r['el1']}_{r['orb1']}-{r['el2']}_{r['orb2']}", axis=1
            )

        return total, orbital

    def _parse_integrated_list(self, path: Path | None, value_col: str) -> pd.DataFrame:
        text = self._read(path)
        if not text:
            return pd.DataFrame()

        rows = []
        for line in text.splitlines():
            parts = line.split()
            if len(parts) < 8 or not parts[0].isdigit():
                continue

            try:
                no = int(parts[0])
                atom1, atom2 = parts[1], parts[2]
                dist = self._float(parts[3])
                tx, ty, tz = map(int, parts[4:7])
                value = self._float(parts[7])
            except Exception:
                continue

            el1, i1, orb1 = self._parse_atom(atom1)
            el2, i2, orb2 = self._parse_atom(atom2)

            rows.append(
                {
                    "source_file": str(path),
                    "interaction_no": no,
                    "atom1": atom1,
                    "atom2": atom2,
                    "el1": el1,
                    "i1": i1,
                    "orb1": orb1,
                    "el2": el2,
                    "i2": i2,
                    "orb2": orb2,
                    "distance_A": dist,
                    "tx": tx,
                    "ty": ty,
                    "tz": tz,
                    value_col: value,
                    "is_total_interaction": orb1 is None and orb2 is None,
                }
            )

        return pd.DataFrame(rows)

    def _parse_energy_file(self, path: Path | None, population: Population) -> tuple[pd.DataFrame, pd.DataFrame] | None:
        text = self._read(path)
        if not text:
            return None

        lines = text.splitlines()
        label_lines = []
        data_start = None

        for idx, line in enumerate(lines[2:], start=2):
            stripped = line.strip()
            if stripped.split() and re.match(r"^[-+]?\d", stripped.split()[0]):
                data_start = idx
                break
            if stripped:
                label_lines.append(stripped)

        if data_start is None:
            return None

        try:
            data = np.loadtxt(path, skiprows=data_start)
        except Exception:
            rows = []
            for line in lines[data_start:]:
                try:
                    rows.append([float(x.replace("D", "E")) for x in line.split()])
                except Exception:
                    pass
            if not rows:
                return None
            data = np.array(rows, dtype=float)

        if data.ndim == 1:
            data = data.reshape(1, -1)

        label_rows = []
        for idx, label in enumerate(label_lines):
            row = self._parse_energy_label(label, idx)
            row["population"] = population
            label_rows.append(row)

        max_pairs = (data.shape[1] - 1) // 2
        columns: dict[str, Any] = {"energy_eV_vs_EF": data[:, 0]}

        for idx in range(min(max_pairs, len(label_lines))):
            columns[f"label_{idx}__{population}_raw"] = data[:, 1 + 2 * idx]
            columns[f"label_{idx}__I{population}_raw"] = data[:, 2 + 2 * idx]

        return pd.DataFrame(label_rows), pd.DataFrame(columns)

    def _parse_energy_label(self, label: str, label_index: int) -> dict[str, Any]:
        out = {
            "label_index": label_index,
            "raw_label": label,
            "is_average": label.strip().lower() == "average",
            "is_total_interaction": False,
            "is_orbital_interaction": False,
            "interaction_no": math.nan,
            "atom1": "",
            "atom2": "",
            "el1": "",
            "el2": "",
            "i1": math.nan,
            "i2": math.nan,
            "orb1": "",
            "orb2": "",
        }

        if out["is_average"]:
            return out

        m = re.match(r"^No\.(\d+):(.+?)->(.+?)\((.*?)\)\s*$", label)
        if not m:
            return out

        atom1 = m.group(2).replace("[", "_").replace("]", "")
        atom2 = m.group(3).replace("[", "_").replace("]", "")
        el1, i1, orb1 = self._parse_atom(atom1)
        el2, i2, orb2 = self._parse_atom(atom2)

        out.update(
            {
                "interaction_no": int(m.group(1)),
                "atom1": atom1,
                "atom2": atom2,
                "el1": el1,
                "el2": el2,
                "i1": i1,
                "i2": i2,
                "orb1": orb1 or "",
                "orb2": orb2 or "",
                "distance_label_A": self._float(m.group(4)),
                "is_total_interaction": orb1 is None and orb2 is None,
                "is_orbital_interaction": orb1 is not None or orb2 is not None,
            }
        )
        return out

    def _parse_band_overlaps(self, path: Path | None) -> pd.DataFrame:
        """Heuristic parser for bandOverlaps.lobster.

        If LOBSTER warns that a k-point could not be orthonormalized, this file is
        worth inspecting. The exact format may vary, so this parser just finds
        numeric square blocks and reports deviation from identity.
        """
        text = self._read(path)
        if not text:
            return pd.DataFrame()

        blocks = []
        current = []
        label = "block_0"
        counter = 0

        def flush():
            nonlocal current, label, counter
            if not current:
                return
            ncols = max(set(map(len, current)), key=list(map(len, current)).count)
            arr = np.array([r for r in current if len(r) == ncols], dtype=float)
            row = {
                "block_index": counter,
                "label": label,
                "n_rows": arr.shape[0],
                "n_cols": arr.shape[1],
                "is_square": arr.shape[0] == arr.shape[1],
                "max_abs_diag_minus_1": math.nan,
                "max_abs_offdiag": math.nan,
                "frobenius_norm_minus_identity": math.nan,
            }
            if row["is_square"] and arr.size:
                offdiag = arr - np.diag(np.diag(arr))
                row["max_abs_diag_minus_1"] = float(np.max(np.abs(np.diag(arr) - 1.0)))
                row["max_abs_offdiag"] = float(np.max(np.abs(offdiag)))
                row["frobenius_norm_minus_identity"] = float(np.linalg.norm(arr - np.eye(arr.shape[0])))
            blocks.append(row)
            current = []
            counter += 1

        float_re = re.compile(r"[-+]?\d*\.\d+(?:[EeDd][-+]?\d+)?|[-+]?\d+(?:[EeDd][-+]?\d+)?")

        for line in text.splitlines():
            lower = line.lower()
            if any(x in lower for x in ["k-point", "kpoint", "spin", "matrix"]):
                flush()
                label = line.strip()
                continue

            nums = [self._float(x) for x in float_re.findall(line)]
            if len(nums) >= 2:
                current.append(nums)
            else:
                flush()

        flush()
        return pd.DataFrame(blocks)

    # =========================================================================
    # Chemistry / aggregation
    # =========================================================================

    def _classify_interaction(
        self, row: pd.Series, topology: TopologyResult | None
    ) -> dict[str, Any]:
        if topology is not None:
            return classify_contact(
                topology,
                int(row["i1"]) if pd.notna(row.get("i1")) else None,
                int(row["i2"]) if pd.notna(row.get("i2")) else None,
                (int(row.get("tx", 0)), int(row.get("ty", 0)), int(row.get("tz", 0))),
                row.get("el1"), row.get("el2"),
            )

        # Safe fallback when no structure could be parsed. Element-group classes
        # remain usable, but topology-dependent X-X contacts are never guessed.
        el1, el2 = row.get("el1"), row.get("el2")
        for rule in self.config.interaction_rules:
            if self._rule_matches(rule, el1, el2):
                if frozenset([el1, el2]).issubset(set(self.config.element_groups.get("anion", []))):
                    cls = rule.unknown_group_name or "X-X_unclassified"
                else:
                    cls = rule.name
                return {
                    "interaction_class": cls,
                    "topology_source": "element-only-fallback",
                    "shared_polyhedron": "",
                    "distance_reconstructed_A": math.nan,
                }
        return {
            "interaction_class": "other",
            "topology_source": "element-only-fallback",
            "shared_polyhedron": "",
            "distance_reconstructed_A": math.nan,
        }

    def _topology_config_from_analysis_config(self) -> TopologyConfig:
        base = {
            "framework_cations": self.config.element_groups.get("framework_cation", ["P", "As", "Sb"]),
            "anions": self.config.element_groups.get("anion", ["S", "Se"]),
            "mobile_cations": self.config.element_groups.get("mobile_cation", ["Na"]),
        }
        base.update(self.config.topology or {})
        return TopologyConfig.from_dict(base)

    def _rule_matches(self, rule: InteractionRule, el1: str | None, el2: str | None) -> bool:
        if not el1 or not el2:
            return False

        pair = frozenset([el1, el2])

        for a, b in rule.element_pairs:
            if pair == frozenset([a, b]):
                return True

        el_groups = self.config.element_to_groups
        groups1 = el_groups.get(el1, set())
        groups2 = el_groups.get(el2, set())

        for g1, g2 in rule.group_pairs:
            if (g1 in groups1 and g2 in groups2) or (g1 in groups2 and g2 in groups1):
                return True

        return False

    @staticmethod
    def _aggregate_integrated(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
        return (
            df.groupby(group_cols, dropna=False)
            .agg(
                n=("interaction_no", "count"),
                distance_min_A=("distance_A", "min"),
                distance_mean_A=("distance_A", "mean"),
                distance_max_A=("distance_A", "max"),
                ICOHP_mean=("ICOHP", "mean"),
                ICOHP_sum=("ICOHP", "sum"),
                minus_ICOHP_mean=("minus_ICOHP", "mean"),
                minus_ICOHP_sum=("minus_ICOHP", "sum"),
                ICOBI_mean=("ICOBI", "mean"),
                ICOBI_sum=("ICOBI", "sum"),
                ICOOP_mean=("ICOOP", "mean"),
                ICOOP_sum=("ICOOP", "sum"),
            )
            .reset_index()
        )

    @staticmethod
    def _bar_long(df: pd.DataFrame, level: str) -> pd.DataFrame:
        if df.empty:
            return pd.DataFrame()

        metrics = [
            "minus_ICOHP_mean", "minus_ICOHP_sum",
            "ICOBI_mean", "ICOBI_sum",
            "ICOOP_mean", "ICOOP_sum",
            "distance_mean_A",
        ]

        label_col = level
        rows = []
        for _, r in df.iterrows():
            for metric in metrics:
                rows.append(
                    {
                        "run_id": r.get("run_id", ""),
                        "analysis_level": level,
                        "bar_label": r.get(label_col, ""),
                        "metric": metric,
                        "value": r.get(metric, math.nan),
                        "n_interactions": r.get("n", math.nan),
                    }
                )
        return pd.DataFrame(rows)

    def _curve_stack(self, curves: pd.DataFrame, population: str, label_indices: Iterable[int], mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        raw, integ = [], []
        for idx in label_indices:
            raw_col = f"label_{idx}__{population}_raw"
            int_col = f"label_{idx}__I{population}_raw"
            if raw_col in curves.columns and int_col in curves.columns:
                raw.append(curves[raw_col].to_numpy(dtype=float)[mask])
                integ.append(curves[int_col].to_numpy(dtype=float)[mask])
        if not raw:
            return np.empty((0, 0)), np.empty((0, 0))
        return np.vstack(raw), np.vstack(integ)

    def _energy_mask(self, energy: np.ndarray) -> np.ndarray:
        mask = np.ones_like(energy, dtype=bool)
        if self.config.energy.plot_energy_min is not None:
            mask &= energy >= self.config.energy.plot_energy_min
        if self.config.energy.plot_energy_max is not None:
            mask &= energy <= self.config.energy.plot_energy_max
        return mask

    @staticmethod
    def _energy_wide(df: pd.DataFrame) -> pd.DataFrame:
        if df.empty:
            return pd.DataFrame()

        wide_parts = []
        for metric in ["mean_plot", "mean_integrated_plot"]:
            temp = df.pivot_table(
                index="energy_eV_vs_EF",
                columns=["run_id", "population", "grouping", "group"],
                values=metric,
                aggfunc="first",
            )
            temp.columns = [f"{metric}__" + "__".join(map(str, c)) for c in temp.columns]
            wide_parts.append(temp)

        return pd.concat(wide_parts, axis=1).reset_index()

    # =========================================================================
    # Diagnostics helpers
    # =========================================================================

    @staticmethod
    def _diagnosis_flag(row: dict[str, Any]) -> str:
        flags = []

        if not row.get("lobsterout_finished_normally", False):
            flags.append("not_finished_normally")

        spill = row.get("lobsterout_abs_charge_spilling_percent", math.nan)
        if pd.isna(spill):
            flags.append("missing_abs_charge_spilling")
        elif spill > 5.0:
            flags.append("high_abs_charge_spilling_gt_5pct")
        elif spill > 2.0:
            flags.append("moderate_abs_charge_spilling_gt_2pct")

        if row.get("lobsterout_has_error", False):
            flags.append("lobsterout_error")
        if row.get("lobsterout_has_warning", False):
            flags.append("lobsterout_warning")
        if row.get("n_total_interactions_parsed", 0) == 0:
            flags.append("no_total_interactions_parsed")
        if not row.get("topology_available", False):
            flags.append("topology_unavailable")
        if row.get("n_topology_distance_mismatch", 0) > 0:
            flags.append("topology_distance_mismatch")
        if row.get("manifest_present", False) and row.get("manifest_structure_hash_matches") is False:
            flags.append("stale_manifest_ignored")

        return ";".join(flags) if flags else "ok"

    @staticmethod
    def _technical_summary(run_quality: pd.DataFrame, inventory: pd.DataFrame, overlap: pd.DataFrame) -> str:
        lines = ["LOBSTER technical summary", ""]

        if not run_quality.empty:
            lines.append("Run quality:")
            for r in run_quality.itertuples(index=False):
                lines.append(
                    f"- {getattr(r, 'run_id')}: "
                    f"flag={getattr(r, 'diagnosis_flag')}, "
                    f"spilling={getattr(r, 'lobsterout_abs_charge_spilling_percent', math.nan)} %, "
                    f"n_total={getattr(r, 'n_total_interactions_parsed')}, "
                    f"populations={getattr(r, 'energy_populations_parsed')}"
                )

        if not inventory.empty:
            lines.extend(["", "Interaction inventory:"])
            for run_id, g in inventory.groupby("run_id"):
                entries = []
                for r in g.itertuples(index=False):
                    entries.append(
                        f"{getattr(r, 'element_pair')}/{getattr(r, 'interaction_class')}: "
                        f"n={getattr(r, 'n')}, mean -ICOHP={getattr(r, 'mean_minus_ICOHP'):.4f}"
                    )
                lines.append(f"- {run_id}: " + "; ".join(entries))

        lines.append("")
        if overlap.empty:
            lines.append("Band overlap matrices: no bandOverlaps file parsed.")
        else:
            lines.append("Band overlap matrices: parsed; inspect band_overlap_quality.csv.")

        return "\n".join(lines)

    @staticmethod
    def _flatten_prefixed(d: dict[str, Any], prefix: str) -> dict[str, Any]:
        out = {}
        for k, v in d.items():
            if isinstance(v, (list, dict)):
                v = json.dumps(v, default=str)
            out[f"{prefix}_{k}"] = v
        return out

    def _write_parsed_summary(self, path: Path) -> None:
        summary = []
        for run in self.runs:
            summary.append(
                {
                    "run_id": run.run_id,
                    "folder": str(run.folder),
                    "files": {k: str(v) if v else "" for k, v in run.files.items()},
                    "formula": run.structure.get("formula", ""),
                    "n_total_interactions": len(run.total_interactions),
                    "n_orbital_interactions": len(run.orbital_interactions),
                    "energy_populations": sorted(run.energy_curves),
                    "topology_available": run.topology is not None,
                    "n_polyhedra": len(run.topology.polyhedra) if run.topology else 0,
                    "manifest_validation": run.manifest_validation,
                }
            )
        path.write_text(json.dumps(summary, indent=2))

    # =========================================================================
    # Small helpers
    # =========================================================================

    def _ensure_parsed(self) -> None:
        if not self._parsed:
            self.parse_files(write_raw=True)

    def _concat(self, attr: str) -> pd.DataFrame:
        frames = [getattr(run, attr) for run in self.runs if not getattr(run, attr).empty]
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    @staticmethod
    def _safe(s: str) -> str:
        return re.sub(r"[^A-Za-z0-9_.-]+", "__", str(s)).strip("_") or "run"

    @staticmethod
    def _read(path: Path | None) -> str:
        if path is None or not path.exists():
            return ""
        return path.read_text(errors="replace")

    @staticmethod
    def _float(x: str | Any) -> float:
        try:
            return float(str(x).replace("D", "E"))
        except Exception:
            return math.nan

    @staticmethod
    def _parse_atom(token: str) -> tuple[str | None, int | None, str | None]:
        clean = token.strip().replace("[", "_").replace("]", "")
        m = re.match(r"^([A-Z][a-z]?)(\d+)(?:_(.*))?$", clean)
        if not m:
            return None, None, None
        return m.group(1), int(m.group(2)), m.group(3)

    @staticmethod
    def _element_pair(el1: str | None, el2: str | None) -> str:
        if not el1 or not el2:
            return "unknown"
        return "-".join(sorted([el1, el2]))


# =============================================================================
# CLI
# =============================================================================


def main() -> None:
    parser = argparse.ArgumentParser(description="System-agnostic LOBSTER bonding analyzer.")
    parser.add_argument("root", type=Path, nargs="?", help="Folder containing one or more LOBSTER runs.")
    parser.add_argument("--out", type=Path, default=Path("lobster_analysis"))
    parser.add_argument("--config", type=Path, help="JSON config file.")
    parser.add_argument("--write-example-config", type=Path, help="Write example JSON config and exit.")
    parser.add_argument("--parse", action="store_true")
    parser.add_argument("--diagnosis", action="store_true")
    parser.add_argument("--integrated", action="store_true")
    parser.add_argument("--energy", action="store_true")
    parser.add_argument("--lobsterpy", action="store_true")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--all-with-lobsterpy", action="store_true")
    parser.add_argument("--individual-curves", action="store_true")
    parser.add_argument(
        "--manifest-mode", choices=["auto", "ignore", "require"], default="auto",
        help="auto: use matching sidecar config; ignore: reconstruct without it; require: fail if missing/stale",
    )

    args = parser.parse_args()

    if args.write_example_config:
        default_na3ps4_like_config().to_json(args.write_example_config)
        print(f"Wrote example config to: {args.write_example_config}")
        return

    if args.root is None:
        raise SystemExit("Please provide a root folder, or use --write-example-config.")

    config = LobsterAnalysisConfig.from_json(args.config) if args.config else default_na3ps4_like_config()
    analyzer = LobsterBondAnalyzer(args.root, args.out, config=config, manifest_mode=args.manifest_mode)

    if args.all_with_lobsterpy:
        analyzer.run_all(include_lobsterpy=True, include_individual_curves=args.individual_curves)
    elif args.all:
        analyzer.run_all(include_lobsterpy=False, include_individual_curves=args.individual_curves)
    else:
        if args.parse:
            analyzer.parse_files()
        if args.diagnosis:
            analyzer.run_diagnosis()
        if args.integrated:
            analyzer.run_integrated_analysis()
        if args.energy:
            analyzer.run_energy_resolved_analysis(include_individual_total_curves=args.individual_curves)
        if args.lobsterpy:
            analyzer.run_lobsterpy_analysis()
        if not any([args.parse, args.diagnosis, args.integrated, args.energy, args.lobsterpy]):
            analyzer.run_all(include_lobsterpy=False, include_individual_curves=args.individual_curves)

    print(f"Done. Output written to: {analyzer.out_dir}")


if __name__ == "__main__":
    main()
