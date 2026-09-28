"""Optional LobsterPy cross-checks; never a source for canonical plot data."""

from __future__ import annotations

import numpy as np
import pandas as pd

from cat_shared.artifacts import TableArtifact

from .classify import element_pair
from .model import LobsterFiles, Population


def _canonical_means(integrals: pd.DataFrame) -> dict[str, float]:
    total = integrals[integrals["orbital"] == "total"].copy()
    per_bond = total.groupby(["bond_label", "element_pair"], as_index=False).agg(
        value=("integrated_value", "sum")
    )
    return per_bond.groupby("element_pair")["value"].mean().to_dict()


def benchmark_artifact(
    files: LobsterFiles,
    population: Population,
    integrals: pd.DataFrame,
    *,
    cutoff: float = 0.0,
    noise_cutoff: float = 1e-8,
    which_bonds: str = "all",
    type_charge: str = "Mulliken",
    absolute_tolerance: float = 1e-5,
    relative_tolerance: float = 1e-4,
) -> TableArtifact:
    """Run LobsterPy independently and compare its pair means to direct extraction."""
    sources = files.populations[population]
    if files.structure is None or not sources.complete:
        raise FileNotFoundError("LobsterPy benchmarking requires structure, COXXCAR, and ICOXXLIST")
    try:
        from lobsterpy.coxx.analyze import Analysis
    except ImportError as exc:
        raise ImportError(
            "Install the optional benchmark with: pip install 'cat-tooling[bonding]'"
        ) from exc
    analysis = Analysis.from_files(
        structure_path=files.structure,
        icoxxlist_path=sources.integrals,
        coxxcar_path=sources.curves,
        charge_path=files.charge,
        are_cobis=population == "COBI",
        are_coops=population == "COOP",
        cutoff_icohp=cutoff,
        noise_cutoff=noise_cutoff,
        which_bonds=which_bonds,
        type_charge=type_charge if files.charge is not None else "Valences",
        orbital_resolved=False,
        summed_spins=True,
    )
    direct = _canonical_means(integrals)
    lobsterpy = {
        element_pair(*str(pair).split("-", 1)): float(values["ICOHP_mean"])
        for pair, values in analysis.final_dict_bonds.items()
    }
    rows = []
    for pair in sorted(set(direct) | set(lobsterpy)):
        left, right = direct.get(pair, np.nan), lobsterpy.get(pair, np.nan)
        delta = right - left if np.isfinite(left) and np.isfinite(right) else np.nan
        matches = bool(
            np.isfinite(left)
            and np.isfinite(right)
            and np.isclose(left, right, atol=absolute_tolerance, rtol=relative_tolerance)
        )
        rows.append(
            {
                "system": files.system,
                "run": files.run,
                "population": population,
                "element_pair": pair,
                "canonical_mean": left,
                "lobsterpy_mean": right,
                "difference": delta,
                "relative_difference": delta / left if np.isfinite(delta) and left != 0 else np.nan,
                "status": "match" if matches else "mismatch",
            }
        )
    return TableArtifact(
        "lobsterpy_benchmark",
        pd.DataFrame(rows),
        ("system", "run", "population", "element_pair"),
        {
            "canonical_mean": "population dependent",
            "lobsterpy_mean": "population dependent",
            "difference": "population dependent",
            "relative_difference": "1",
        },
        {
            "purpose": "diagnostic_only",
            "canonical_source": "bond_integrals",
            "cutoff": cutoff,
            "noise_cutoff": noise_cutoff,
            "which_bonds": which_bonds,
            "absolute_tolerance": absolute_tolerance,
            "relative_tolerance": relative_tolerance,
        },
    )
