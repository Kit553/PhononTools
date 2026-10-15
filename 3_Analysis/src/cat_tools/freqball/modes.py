"""Mode-resolved structural-direction weights for mesh and path fatbands."""

from __future__ import annotations

import re
from collections.abc import Callable, Sequence

import numpy as np
from cat_shared.artifacts import TableArtifact

from cat_tools.phonopy.modes import ModeData, mode_artifact

from .families import FamilySet
from .named_projections import projection_metadata


def _eigenvectors(values: np.ndarray, atom_count: int) -> np.ndarray:
    array = np.asarray(values)
    modes = 3 * atom_count
    if array.ndim == 4 and array.shape[-2:] == (atom_count, 3):
        return array
    if array.ndim != 3 or array.shape[1:] != (modes, modes):
        raise ValueError(f"Unexpected eigenvector shape {array.shape}")
    return array.transpose(0, 2, 1).reshape(len(array), modes, atom_count, 3)


def _family_weight(
    eigenvectors: np.ndarray,
    selected: np.ndarray,
    families: FamilySet,
    predicate: Callable,
) -> np.ndarray:
    total = np.zeros(eigenvectors.shape[:2], dtype=float)
    found = False
    for primitive_index in selected:
        references = [
            item.cartesian
            for item in families.atom(int(primitive_index)).references
            if predicate(item)
        ]
        if not references:
            continue
        found = True
        projected = np.einsum(
            "qbc,rc->qbr",
            eigenvectors[:, :, primitive_index],
            np.asarray(references),
            optimize=True,
        )
        total += np.max(abs(projected) ** 2, axis=-1)
    return total if found else np.full_like(total, np.nan)


def _safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")


def directional_mode_artifact(
    data: ModeData,
    *,
    system: str,
    symbols: Sequence[str],
    masses: Sequence[float],
    selected_atoms: Sequence[int],
    families: FamilySet,
    name: str | None = None,
    projections=(),
) -> TableArtifact:
    """Add absolute and conditional local-family weights to a standard mode table."""
    if data.eigenvectors is None:
        raise RuntimeError("FreqBall mode analysis requires eigenvectors")
    selected = np.asarray(selected_atoms, dtype=int)
    eigenvectors = _eigenvectors(data.eigenvectors, len(symbols))
    mobile = np.sum(abs(eigenvectors[:, :, selected]) ** 2, axis=(2, 3))
    conduction = _family_weight(
        eigenvectors, selected, families, lambda item: item.category == "conduction"
    )
    tetrahedral = _family_weight(
        eigenvectors, selected, families, lambda item: item.category == "tetrahedral"
    )
    base = mode_artifact(
        data,
        system=system,
        symbols=symbols,
        masses=masses,
        include_participation_ratio=True,
        projections=projections,
    )
    frame = base.data.copy()
    frame["mobile_participation"] = mobile.reshape(-1)
    frame["conduction_direction_weight"] = conduction.reshape(-1)
    frame["tetrahedral_direction_weight"] = tetrahedral.reshape(-1)
    frame["conduction_direction_fraction"] = np.divide(
        conduction,
        mobile,
        out=np.full_like(mobile, np.nan),
        where=mobile > 0,
    ).reshape(-1)
    frame["tetrahedral_direction_fraction"] = np.divide(
        tetrahedral,
        mobile,
        out=np.full_like(mobile, np.nan),
        where=mobile > 0,
    ).reshape(-1)
    frame["direction_score"] = (
        frame["conduction_direction_fraction"] - frame["tetrahedral_direction_fraction"]
    )
    frame["dominant_direction_family"] = np.where(
        frame["direction_score"].isna(),
        "unassigned",
        np.where(frame["direction_score"] >= 0, "conduction", "tetrahedral"),
    )
    subtypes = sorted(
        {
            item.family
            for primitive_index in selected
            for item in families.atom(int(primitive_index)).references
        }
    )
    for subtype in subtypes:
        values = _family_weight(
            eigenvectors,
            selected,
            families,
            lambda item, subtype=subtype: item.family == subtype,
        )
        frame[f"{_safe_name(subtype)}_direction_weight"] = values.reshape(-1)
    return TableArtifact(
        name or f"freqball_{data.source}_modes",
        frame,
        base.keys,
        base.units,
        {
            **base.metadata,
            "selected_primitive_indices": selected.tolist(),
            "family_source": families.source,
            "projections": projection_metadata(projections),
            "weight_definition": "sum_atoms(max_references(|direction dot eigenvector|^2))",
        },
    )
