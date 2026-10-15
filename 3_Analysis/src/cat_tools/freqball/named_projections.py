"""Named projections shared by tensor spectra, descriptors and mode tables."""

import re
from collections.abc import Mapping

import numpy as np
import pandas as pd
from cat_shared.artifacts import TableArtifact

from cat_tools.phonopy.projections import Projection

from .analysis import spectral_descriptors
from .directions import normalize_rows
from .tensor import directional_spectra


def cartesian_direction(value, lattice) -> tuple[float, float, float]:
    """Convert a Cartesian vector, direct-lattice axis or Miller plane normal."""
    frame = "cartesian"
    if isinstance(value, Mapping):
        if len(value) != 1 or next(iter(value)) not in {"cartesian", "fractional", "miller"}:
            raise ValueError("Direction needs exactly one of cartesian, fractional or miller")
        frame, value = next(iter(value.items()))
    vector = np.asarray(value, dtype=float)
    if vector.shape != (3,) or not np.isfinite(vector).all():
        raise ValueError("Direction must be a finite three-vector")
    if frame == "fractional":
        vector = vector @ np.asarray(lattice)
    elif frame == "miller":
        vector = vector @ np.linalg.inv(lattice).T
    return tuple(float(x) for x in normalize_rows(vector)[0])


def resolve_projections(config, *, symbols, selected, lattice) -> list[Projection]:
    """Resolve primitive indices once, with identical meaning in every output.

    Only tensor-selected sites are available. An explicit empty or unavailable
    selection is an error, never an implicit request for every atom.
    """
    if config is None:
        return []
    if not isinstance(config, Mapping):
        raise TypeError("projections must be a mapping")
    unknown = set(config) - {"elements", "sites", "directions", "custom"}
    if unknown:
        raise ValueError(f"Unknown projection settings: {sorted(unknown)}")
    selected = tuple(int(i) for i in selected)
    allowed = set(selected)
    items = []
    if config.get("elements", False):
        for symbol in dict.fromkeys(symbols[i] for i in selected):
            items.append({"name": symbol, "elements": [symbol]})
    if config.get("sites", False):
        items.extend({"name": f"site_{i}_{symbols[i]}", "atoms": [i]} for i in selected)
    items.extend(
        {"name": name, "direction": vector} for name, vector in config.get("directions", {}).items()
    )
    items.extend(config.get("custom", []))
    projections = []
    names = set()
    for item in items:
        unknown = set(item) - {"name", "atoms", "elements", "direction"}
        if unknown:
            raise ValueError(f"Unknown custom projection settings: {sorted(unknown)}")
        name = item.get("name")
        if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name):
            raise ValueError("Projection names must start with a letter and use letters/digits/_")
        if name in names or name == "mobile":
            raise ValueError(f"Duplicate or reserved projection name: {name}")
        names.add(name)
        if "atoms" in item and "elements" in item:
            raise ValueError(f"{name}: choose atoms or elements")
        if "atoms" in item:
            atoms = tuple(item["atoms"])
            if any(isinstance(i, bool) or not isinstance(i, (int, np.integer)) for i in atoms):
                raise ValueError(f"{name}: atoms must be zero-based integer primitive indices")
        elif "elements" in item:
            elements = set(item["elements"])
            missing = elements - {symbols[i] for i in selected}
            if missing:
                raise ValueError(
                    f"{name}: elements absent from tensor selection: {sorted(missing)}"
                )
            atoms = tuple(i for i in selected if symbols[i] in elements)
        else:
            atoms = selected
        if not atoms or len(set(atoms)) != len(atoms) or not set(atoms) <= allowed:
            raise ValueError(
                f"{name}: atoms must be nonempty, unique and included in atoms selection"
            )
        direction = cartesian_direction(item["direction"], lattice) if "direction" in item else None
        projections.append(Projection(name, atoms, direction))
    return projections


def projection_metadata(projections):
    return [
        {
            "name": p.name,
            "primitive_indices": list(p.atoms),
            "cartesian_direction": None if p.direction is None else list(p.direction),
        }
        for p in projections
    ]


def projection_artifacts(result, projections, *, low_cutoff=2.0, inverse_floor=0.1):
    """Sum selected tensor sites; a direction-free projection uses the trace."""
    if not projections:
        return []
    index = {int(p): i for i, p in enumerate(result.primitive_indices)}
    spectra, descriptors = [], []
    for projection in projections:
        tensor = result.tensor[[index[i] for i in projection.atoms]]
        values = (
            tensor[:, :3].sum(axis=(0, 1))
            if projection.direction is None
            else directional_spectra(tensor, projection.direction)[0].sum(axis=0)
        )
        spectra.append(
            pd.DataFrame(
                {
                    "system": result.system,
                    "projection": projection.name,
                    "frequency_index": np.arange(len(values)),
                    "frequency_THz": result.frequency_THz,
                    "pdos_per_THz": values,
                }
            )
        )
        metrics = spectral_descriptors(
            result.frequency_THz,
            values[None, :],
            low_frequency_cutoff=low_cutoff,
            inverse_frequency_floor=inverse_floor,
        )
        descriptors.append(
            {
                "system": result.system,
                "projection": projection.name,
                **{key: value[0] for key, value in metrics.items()},
            }
        )
    metadata = {
        "projections": projection_metadata(projections),
        "aggregation": "sum_atoms",
        "projection_weighting": "eigenvector",
        "low_frequency_cutoff_THz": low_cutoff,
        "inverse_frequency_floor_THz": inverse_floor,
    }
    units = {
        key: "THz^2" if key.endswith("_THz2") else "THz"
        for key in descriptors[0]
        if key.endswith(("_THz", "_THz2"))
    }
    units.update({"inverse_frequency_weight": "THz^-2", "inverse_frequency_fraction": "THz^-2"})
    return [
        TableArtifact(
            "freqball_projection_spectra",
            pd.concat(spectra, ignore_index=True),
            ("projection", "frequency_index"),
            {"frequency_THz": "THz", "pdos_per_THz": "THz^-1"},
            metadata,
        ),
        TableArtifact(
            "freqball_projection_descriptors",
            pd.DataFrame(descriptors),
            ("projection",),
            units,
            metadata,
        ),
    ]
