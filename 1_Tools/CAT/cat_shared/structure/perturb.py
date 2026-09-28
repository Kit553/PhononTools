"""Random atomic displacements for independent relaxation starts."""

from __future__ import annotations

import numpy as np


def perturb_structure(structure, amplitude: float = 0.02, *, seed: int | None = None):
    """Return a copy with each mobile atom displaced by ``amplitude`` angstrom.

    Directions are uniform in the allowed Cartesian subspace. VASP selective
    dynamics flags constrain fractional coordinates, including in skew cells.
    Fully fixed sites remain unchanged. No wrapping or standardization is done.
    """
    if not np.isfinite(amplitude) or amplitude <= 0:
        raise ValueError("amplitude must be finite and positive (angstrom)")
    if seed is not None and (not isinstance(seed, (int, np.integer)) or seed < 0):
        raise ValueError("seed must be a non-negative integer")
    allowed = np.asarray(
        structure.site_properties.get("selective_dynamics", [[True] * 3] * len(structure)),
        dtype=bool,
    )
    if allowed.shape != (len(structure), 3):
        raise ValueError("selective_dynamics must contain three flags per site")
    if not allowed.any():
        raise ValueError("structure has no mobile coordinates to perturb")
    rng = np.random.default_rng(seed)
    result = structure.copy()
    for index, mask in enumerate(allowed):
        if not mask.any():
            continue
        # Orthonormal basis for the span of the allowed direct lattice vectors.
        basis, _ = np.linalg.qr(structure.lattice.matrix[mask].T)
        direction = basis @ rng.normal(size=int(mask.sum()))
        direction /= np.linalg.norm(direction)
        result.translate_sites(
            [index], amplitude * direction, frac_coords=False, to_unit_cell=False
        )
    return result
