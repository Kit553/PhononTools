"""Normalize the v6 and CAT FreqBall input dialects without loading Phonopy."""

from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path


def normalize_config(document: Mapping) -> dict:
    """Retain v6 physical settings while using the canonical tensor output."""
    config = deepcopy(dict(document))
    if "phonopy" not in config:
        return config
    if "system" in config:
        raise ValueError("Do not mix v6 phonopy and CAT system configuration sections")
    phonopy = config.pop("phonopy")
    sampling = config.pop("sampling")
    selection = config.pop("selection", {})
    output = config.pop("output", {})
    config["system"] = {
        "path": ".",
        **{
            key: phonopy[key]
            for key in ("phonopy_yaml", "force_sets", "born", "masses")
            if key in phonopy
        },
    }
    config["atoms"] = config.pop("ions")
    config["spectrum"] = {
        "mesh": phonopy["mesh"],
        "frequency_range_THz": [sampling["freq_min"], sampling["freq_max"]],
        "pitch_THz": sampling["freq_pitch"],
        "sigma_THz": sampling["sigma"],
        # v6 uses Phonopy's default shifted even mesh.
        "gamma_center": False,
    }
    config["analysis"] = {
        "directions": {
            "count": sampling["n_directions"],
            "low_frequency_cutoff_THz": selection.get("low_freq_cutoff", 2.0),
        },
        "family_spectra": True,
    }
    family = config.pop("direction_families_file", None)
    tolerance = config.pop("direction_family_match_tolerance_A", 0.02)
    if family:
        config["families"] = {"file": family, "match_tolerance_A": tolerance}
    config["output"] = {
        "tensor": output.get("file", "freqball.h5"),
        "directory": "freqball-analysis",
        "spectra_dtype": output.get("spectra_dtype", "float32"),
    }
    return config


def anchor_inputs(config: Mapping, base: Path) -> dict:
    """Make external inputs portable when writing a configuration elsewhere."""
    result = deepcopy(dict(config))

    def absolute(value):
        path = Path(value).expanduser()
        return str((path if path.is_absolute() else base / path).resolve())

    result["system"]["path"] = absolute(result["system"].get("path", "."))
    family = result.get("families")
    if family and family.get("file"):
        family["file"] = absolute(family["file"])
    modes = result.get("modes") or []
    for mode in [modes] if isinstance(modes, dict) else modes:
        if mode.get("from_table"):
            mode["from_table"] = absolute(mode["from_table"])
    return result
