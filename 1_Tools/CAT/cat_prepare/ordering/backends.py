"""Thin icet adapters; internal proxy symbols represent categorical occupations.

Motif SQS optimizes correlations of motif identities at their centers. It is not
an independently random alloy on the individual ligand sites. All other scaffold
positions remain present to retain the geometry's symmetry constraints.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import linear_sum_assignment

from .generate import Candidate, pair_features
from .model import build_model


def _imports():
    try:
        from ase import Atoms
        from ase.data import chemical_symbols
        from icet import ClusterSpace
        from icet.tools import enumerate_structures
        from icet.tools.structure_generation import generate_sqs_from_supercells
    except ImportError as exc:
        raise ImportError(
            "Install ordering backends with: pip install 'cat-tooling[ordering]'"
        ) from exc
    return Atoms, chemical_symbols, ClusterSpace, enumerate_structures, generate_sqs_from_supercells


def proxy_lattice(model):
    Atoms, symbols, *_ = _imports()
    tags = {}
    for gi, group in enumerate(model.groups):
        for unit in group.units:
            tags[unit[0]] = [("state", gi, label) for label in sorted(group.states)]
            for index in unit[1:]:
                tags[index] = [("ligand", gi)]
    all_tags = []
    for i, site in enumerate(model.structure):
        all_tags.append(tags.get(i, [("fixed", site.specie.symbol)]))
    keys = sorted({tag for choices in all_tags for tag in choices}, key=repr)
    if len(keys) >= len(symbols):
        raise ValueError("icet proxy alphabet exceeds 118 distinct categorical states")
    mapping = {tag: symbols[i + 1] for i, tag in enumerate(keys)}
    choices = [[mapping[tag] for tag in ts] for ts in all_tags]
    atoms = Atoms(
        [cs[0] for cs in choices],
        scaled_positions=model.structure.frac_coords,
        cell=model.structure.lattice.matrix,
        pbc=True,
    )
    return atoms, choices, mapping


def _assignment(model, atoms, mapping):
    distances = model.structure.lattice.get_all_distances(
        model.structure.frac_coords, atoms.get_scaled_positions()
    )
    ii, jj = linear_sum_assignment(distances)
    if len(ii) != len(model.structure) or distances[ii, jj].max() > 1e-4:
        raise ValueError("icet output does not map back to the requested parent scaffold")
    species = atoms.get_chemical_symbols()
    result = []
    for gi, group in enumerate(model.groups):
        reverse = {mapping[("state", gi, label)]: label for label in group.states}
        labels = [reverse[species[jj[u[0]]]] for u in group.units]
        if {k: labels.count(k) for k in group.states} != group.counts:
            raise ValueError(f"icet output violates exact composition on {group.name}")
        result.append(labels)
    return result


def generate_icet_enumeration(
    parent,
    specifications,
    sizes,
    *,
    max_candidates=10000,
    max_enumerated=100000,
    cutoff=6.0,
    symprec=1e-5,
):
    """Symmetry-inequivalent decorations and supercell shapes from icet.

    Fractions are restricted through unique per-group proxy symbols. The explicit
    bounds fail rather than silently claiming a truncated enumeration is complete.
    """
    *_, enumerate_structures, _sqs = _imports()
    if not sizes or any(type(n) is not int or n < 1 for n in sizes):
        raise ValueError("icet_sizes must be positive integer multiples of the input cell")
    base = build_model(parent, specifications, validate_counts=False)
    atoms, choices, mapping = proxy_lattice(base)
    restrictions = {}
    for gi, group in enumerate(base.groups):
        for label, fraction in group.fractions.items():
            c = fraction * len(group.units) / len(parent)
            restrictions[mapping[("state", gi, label)]] = (max(0, c - 1e-10), min(1, c + 1e-10))
    output = []
    for number, decorated in enumerate(
        enumerate_structures(
            atoms,
            sizes,
            choices,
            concentration_restrictions=restrictions,
            niggli_reduce=False,
            symprec=symprec,
        )
    ):
        if number >= max_enumerated:
            raise ValueError("icet enumeration exceeds max_enumerated; reduce the search")
        matrix = decorated.cell.array @ np.linalg.inv(parent.lattice.matrix)
        if not np.allclose(matrix, np.rint(matrix), atol=1e-5):
            raise ValueError("icet returned a nonintegral parent-cell transformation")
        model = build_model(parent, specifications, np.rint(matrix).astype(int))
        assignments = _assignment(model, decorated, mapping)
        if len(output) >= max_candidates:
            raise ValueError("icet enumeration exceeds max_candidates; reduce the search")
        output.append(
            Candidate(
                model.decorate(assignments),
                [
                    {
                        "method": "icet-enumerate",
                        "matrix": model.matrix,
                        "assignments": {
                            g.name: a for g, a in zip(model.groups, assignments, strict=True)
                        },
                        "multiplicity_is_exact": False,
                    }
                ],
                pair_features(model, assignments, cutoff),
            )
        )
    if not output:
        raise ValueError("icet found no configurations at the requested sizes and fractions")
    return output


def generate_sqs(
    model, *, cutoffs, seed=0, steps=10000, symprec=1e-5, pair_cutoff=6.0, target_vector=None
):
    """Actual icet cluster-correlation optimization, with diagnostics saved.

    An optional target vector permits prescribed short-range correlations; such a
    result is identified as a target-correlation structure rather than an SQS.
    """
    _, _, ClusterSpace, _, sqs = _imports()
    from icet.tools.structure_generation import (
        _get_sqs_cluster_vector,
        compare_cluster_vectors,
        generate_target_structure_from_supercells,
    )

    if not cutoffs or any(not np.isfinite(x) or x <= 0 for x in cutoffs) or steps < 1:
        raise ValueError("SQS requires positive cluster cutoffs and n_steps")
    if not any(sum(n > 0 for n in g.counts.values()) > 1 for g in model.groups):
        raise ValueError("SQS needs at least two occupied states on an active sublattice")
    atoms, choices, mapping = proxy_lattice(model)
    cs = ClusterSpace(atoms, cutoffs, choices, symprec=symprec)
    by_symbol = {}
    for gi, group in enumerate(model.groups):
        by_symbol.update({mapping[("state", gi, k)]: v for k, v in group.fractions.items()})
    concentrations = {}
    for sub in cs.get_sublattices(atoms):
        if len(sub.chemical_symbols) > 1:
            concentrations[sub.symbol] = {s: by_symbol[s] for s in sub.chemical_symbols}
        else:
            concentrations[sub.symbol] = {sub.chemical_symbols[0]: 1.0}
    random_vector = _get_sqs_cluster_vector(cs, concentrations)
    kwargs = {
        "cluster_space": cs,
        "supercells": [atoms],
        "target_concentrations": concentrations,
        "n_steps": steps,
        "random_seed": seed,
    }
    if target_vector is None:
        decorated = sqs(**kwargs)
        target = random_vector
    else:
        target = np.asarray(target_vector, dtype=float)
        if target.shape != random_vector.shape or not np.isfinite(target).all():
            raise ValueError(f"target_vector must contain {len(random_vector)} finite values")
        decorated = generate_target_structure_from_supercells(
            **kwargs, target_cluster_vector=target
        )
    actual = cs.get_cluster_vector(decorated)
    assignments = _assignment(model, decorated, mapping)
    metadata = {
        "method": "sqs" if target_vector is None else "target-correlations",
        "backend": "icet",
        "seed": seed,
        "steps": steps,
        "matrix": model.matrix,
        "cutoffs_A": list(cutoffs),
        "cluster_vector": actual.tolist(),
        "target_vector": target.tolist(),
        "random_target_vector": random_vector.tolist(),
        "cluster_space": str(cs),
        "objective": float(compare_cluster_vectors(actual, target, cs.as_list)),
        "max_abs_correlation_error": float(np.max(np.abs(actual - target))),
        "proxy_symbols": {repr(k): v for k, v in mapping.items()},
        "target_concentrations": concentrations,
        "assignments": {g.name: a for g, a in zip(model.groups, assignments, strict=True)},
        "scope": "motif centers"
        if any(len(g.units[0]) > 1 for g in model.groups)
        else "atomic occupations",
        "multiplicity_is_exact": False,
    }
    return Candidate(
        model.decorate(assignments), [metadata], pair_features(model, assignments, pair_cutoff)
    )
