"""Canonical cell and site tables shared by static structures and MD frames."""

from __future__ import annotations

import json
from collections.abc import Iterable

import numpy as np
import pandas as pd

from cat_shared.artifacts import TableArtifact

from .io import _structure_class
from .model import StructureFrame


def _plain(value):
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_plain(item) for item in value]
    if hasattr(value, "item"):
        return value.item()
    return value


def _number(value) -> float:
    return np.nan if value is None else float(value)


def structure_artifacts(frames: Iterable[StructureFrame]) -> tuple[TableArtifact, TableArtifact]:
    """Create the two loss-aware tables needed to reconstruct periodic structures."""
    cells, sites = [], []
    seen = set()
    for item in frames:
        item.validate()
        identity = (item.system, item.frame)
        if identity in seen:
            raise ValueError(f"Duplicate structure frame {identity}")
        seen.add(identity)
        structure = item.structure
        lattice = np.asarray(structure.lattice.matrix, dtype=float)
        cells.append(
            {
                "system": item.system,
                "frame": item.frame,
                "time_fs": _number(item.time_fs),
                "energy_eV": _number(item.energy_eV),
                "temperature_K": _number(item.temperature_K),
                "a_x": lattice[0, 0],
                "a_y": lattice[0, 1],
                "a_z": lattice[0, 2],
                "b_x": lattice[1, 0],
                "b_y": lattice[1, 1],
                "b_z": lattice[1, 2],
                "c_x": lattice[2, 0],
                "c_y": lattice[2, 1],
                "c_z": lattice[2, 2],
                "pbc_a": bool(structure.lattice.pbc[0]),
                "pbc_b": bool(structure.lattice.pbc[1]),
                "pbc_c": bool(structure.lattice.pbc[2]),
                "charge_e": _number(structure.charge),
                "volume_A3": float(structure.volume),
                "properties_json": json.dumps(_plain(structure.properties), sort_keys=True),
                "source": item.source,
            }
        )
        velocities = None if item.velocities_A_fs is None else np.asarray(item.velocities_A_fs)
        forces = None if item.forces_eV_A is None else np.asarray(item.forces_eV_A)
        for site_index, site in enumerate(structure):
            fractional = np.asarray(site.frac_coords, dtype=float)
            cartesian = np.asarray(site.coords, dtype=float)
            common = {
                "system": item.system,
                "frame": item.frame,
                "site": site_index,
                "label": site.label or f"site-{site_index}",
                "frac_x": fractional[0],
                "frac_y": fractional[1],
                "frac_z": fractional[2],
                "cart_x_A": cartesian[0],
                "cart_y_A": cartesian[1],
                "cart_z_A": cartesian[2],
                "velocity_x_A_fs": np.nan if velocities is None else velocities[site_index, 0],
                "velocity_y_A_fs": np.nan if velocities is None else velocities[site_index, 1],
                "velocity_z_A_fs": np.nan if velocities is None else velocities[site_index, 2],
                "force_x_eV_A": np.nan if forces is None else forces[site_index, 0],
                "force_y_eV_A": np.nan if forces is None else forces[site_index, 1],
                "force_z_eV_A": np.nan if forces is None else forces[site_index, 2],
                "site_properties_json": json.dumps(_plain(site.properties), sort_keys=True),
            }
            for species, occupancy in site.species.items():
                sites.append(
                    {
                        **common,
                        "species": str(species),
                        "element": species.symbol,
                        "occupancy": float(occupancy),
                    }
                )
    return (
        TableArtifact(
            "structure_cells",
            pd.DataFrame(cells),
            ("system", "frame"),
            {
                "time_fs": "fs",
                "energy_eV": "eV",
                "temperature_K": "K",
                "a_x": "Angstrom",
                "a_y": "Angstrom",
                "a_z": "Angstrom",
                "b_x": "Angstrom",
                "b_y": "Angstrom",
                "b_z": "Angstrom",
                "c_x": "Angstrom",
                "c_y": "Angstrom",
                "c_z": "Angstrom",
                "charge_e": "e",
                "volume_A3": "Angstrom^3",
            },
            {"grain": "one periodic cell per structure or trajectory frame"},
        ),
        TableArtifact(
            "structure_sites",
            pd.DataFrame(sites),
            ("system", "frame", "site", "species"),
            {
                "cart_x_A": "Angstrom",
                "cart_y_A": "Angstrom",
                "cart_z_A": "Angstrom",
                "velocity_x_A_fs": "Angstrom/fs",
                "velocity_y_A_fs": "Angstrom/fs",
                "velocity_z_A_fs": "Angstrom/fs",
                "force_x_eV_A": "eV/Angstrom",
                "force_y_eV_A": "eV/Angstrom",
                "force_z_eV_A": "eV/Angstrom",
            },
            {
                "site_indexing": "zero_based",
                "grain": "one species component per crystallographic site and frame",
            },
        ),
    )


def _optional_vector(rows: pd.DataFrame, names: list[str]):
    values = rows.groupby("site", sort=True)[names].first().to_numpy(float)
    return None if np.isnan(values).all() else values


def frames_from_tables(
    cells: pd.DataFrame,
    sites: pd.DataFrame,
    *,
    system: str | None = None,
    frame: int | None = None,
) -> list[StructureFrame]:
    """Reconstruct Pymatgen structures from canonical tables."""
    Structure = _structure_class()
    selected_cells = cells.copy()
    if system is not None:
        selected_cells = selected_cells[selected_cells["system"] == system]
    if frame is not None:
        selected_cells = selected_cells[selected_cells["frame"] == frame]
    if selected_cells.empty:
        raise ValueError("No structure cell matches the requested system/frame")
    output = []
    for _, cell in selected_cells.sort_values(["system", "frame"]).iterrows():
        rows = sites[(sites["system"] == cell["system"]) & (sites["frame"] == cell["frame"])]
        if rows.empty:
            raise ValueError(f"No site rows for {cell['system']} frame {cell['frame']}")
        species, coordinates, labels, properties = [], [], [], []
        for _, group in rows.sort_values(["site", "species"]).groupby("site", sort=True):
            reference = group.iloc[0]
            if not np.allclose(
                group[["frac_x", "frac_y", "frac_z"]],
                reference[["frac_x", "frac_y", "frac_z"]].to_numpy(float),
            ):
                raise ValueError("Species components on one site have different coordinates")
            species.append(dict(zip(group["species"], group["occupancy"], strict=True)))
            coordinates.append(reference[["frac_x", "frac_y", "frac_z"]].to_numpy(float))
            labels.append(str(reference["label"]))
            properties.append(json.loads(reference.get("site_properties_json", "{}") or "{}"))
        property_names = sorted(set().union(*(item.keys() for item in properties)))
        site_properties = {
            name: [item.get(name) for item in properties]
            for name in property_names
            if all(name in item for item in properties)
        }
        lattice = cell[["a_x", "a_y", "a_z", "b_x", "b_y", "b_z", "c_x", "c_y", "c_z"]]
        charge = None if pd.isna(cell.get("charge_e")) else float(cell["charge_e"])
        structure = Structure(
            np.asarray(lattice, dtype=float).reshape(3, 3),
            species,
            coordinates,
            charge=charge,
            site_properties=site_properties,
            labels=labels,
            properties=json.loads(cell.get("properties_json", "{}") or "{}"),
        )
        time = None if pd.isna(cell.get("time_fs")) else float(cell["time_fs"])
        energy = None if pd.isna(cell.get("energy_eV")) else float(cell["energy_eV"])
        temperature = None if pd.isna(cell.get("temperature_K")) else float(cell["temperature_K"])
        output.append(
            StructureFrame(
                str(cell["system"]),
                int(cell["frame"]),
                structure,
                time,
                energy,
                temperature,
                _optional_vector(
                    rows,
                    ["velocity_x_A_fs", "velocity_y_A_fs", "velocity_z_A_fs"],
                ),
                _optional_vector(rows, ["force_x_eV_A", "force_y_eV_A", "force_z_eV_A"]),
                str(cell.get("source", "")),
            )
        )
    return output
