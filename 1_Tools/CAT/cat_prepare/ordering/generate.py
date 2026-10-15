"""Bounded enumeration, periodic motif arrangements, and symmetry reduction."""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations, product
from math import prod

import numpy as np
import spglib
from scipy.optimize import linear_sum_assignment


@dataclass
class Candidate:
    structure: object
    provenance: list[dict] = field(default_factory=list)
    descriptors: dict = field(default_factory=dict)


def _multisets(counts):
    """Bounded fixed-composition multiset permutations, without factorial duplicates."""
    labels = sorted(k for k, n in counts.items() if n)
    assignment = [None] * sum(counts.values())

    def visit(level, available):
        label = labels[level]
        if level == len(labels) - 1:
            for i in available:
                assignment[i] = label
            yield tuple(assignment)
            return
        for chosen in combinations(available, counts[label]):
            for i in chosen:
                assignment[i] = label
            occupied = set(chosen)
            yield from visit(level + 1, tuple(i for i in available if i not in occupied))

    yield from visit(0, tuple(range(len(assignment))))


def symmetry_permutations(model, symprec=1e-5):
    """Use spglib operations that preserve all variable groups and motif membership."""
    s = model.structure
    colors, lookup = [], {}
    variable_roles = {}
    for g, group in enumerate(model.groups):
        for unit in group.units:
            for role, i in enumerate(unit):
                # All equivalent ligands in the supported homogeneous motifs share a role.
                variable_roles[i] = ("variable", g, min(role, 1))
    for i, site in enumerate(s):
        tag = variable_roles.get(i, ("fixed", site.specie.symbol))
        if tag not in lookup:
            lookup[tag] = len(lookup) + 1
        colors.append(lookup[tag])
    data = spglib.get_symmetry((s.lattice.matrix, s.frac_coords, colors), symprec=symprec)
    if data is None:
        raise ValueError("spglib could not determine parent symmetry")
    colors = np.array(colors)
    permutations = []
    for rotation, translation in zip(data["rotations"], data["translations"], strict=True):
        transformed = s.frac_coords @ rotation.T + translation
        distances = s.lattice.get_all_distances(transformed, s.frac_coords)
        distances[colors[:, None] != colors[None, :]] = 1e10
        ii, jj = linear_sum_assignment(distances)
        if distances[ii, jj].max() > symprec * 2:
            continue
        mapping = dict(zip(ii, jj, strict=True))
        group_perms = []
        for group in model.groups:
            unit_lookup = {frozenset(u): i for i, u in enumerate(group.units)}
            perm = [unit_lookup.get(frozenset(mapping[i] for i in u)) for u in group.units]
            if None in perm:
                break
            group_perms.append(perm)
        if len(group_perms) == len(model.groups):
            permutations.append(group_perms)
    if not permutations:
        raise ValueError("No valid occupational symmetry operations (including identity)")
    return permutations


def canonical_assignment(assignments, permutations):
    keys = []
    for perms in permutations:
        key = []
        for labels, perm in zip(assignments, perms, strict=True):
            transformed = [None] * len(labels)
            for i, j in enumerate(perm):
                transformed[j] = labels[i]
            key.extend(transformed)
        keys.append(tuple(key))
    return min(keys)


def pair_features(model, assignments, cutoff):
    """Periodic directed center-pair fractions, including all images within cutoff.

    These describe occupational correlations, not formation energies. For motifs
    the centers represent complete units. Cross-group correlations are included.
    """
    from pymatgen.core import Structure

    centers, groups, labels = [], [], []
    for group, assignment in zip(model.groups, assignments, strict=True):
        centers.extend(model.structure[u[0]].frac_coords for u in group.units)
        groups.extend([group.name] * len(group.units))
        labels.extend(assignment)
    center_structure = Structure(model.structure.lattice, ["H"] * len(centers), centers)
    ci, pi, _, _ = center_structure.get_neighbor_list(cutoff)
    totals, counts = {}, {}
    for i, j in zip(ci, pi, strict=True):
        pair = tuple(sorted(((groups[i], labels[i]), (groups[j], labels[j]))))
        key = f"pair:{pair[0][0]}/{pair[0][1]}--{pair[1][0]}/{pair[1][1]}"
        group_key = tuple(sorted((groups[i], groups[j])))
        counts[key] = counts.get(key, 0) + 1
        totals[group_key] = totals.get(group_key, 0) + 1
    features = {}
    for ga in model.groups:
        for gb in model.groups:
            for a in ga.states:
                for b in gb.states:
                    pair = tuple(sorted(((ga.name, a), (gb.name, b))))
                    key = f"pair:{pair[0][0]}/{pair[0][1]}--{pair[1][0]}/{pair[1][1]}"
                    total = totals.get(tuple(sorted((ga.name, gb.name))), 0)
                    features[key] = counts.get(key, 0) / total if total else 0.0
    return features


def generate_fixed(
    model,
    method,
    *,
    seed=0,
    samples=8,
    steps=1000,
    cutoff=6.0,
    symprec=1e-5,
    max_raw=100000,
    axis=2,
):
    """Generate ordered structures by enumeration, swaps or correlation-based sampling.

    Enumeration covers the supplied supercell. Other methods sample candidate orderings without establishing completeness or energetic stability.
    """
    if cutoff <= 0 or samples < 1 or steps < 0 or max_raw < 1:
        raise ValueError("Positive cutoff, samples and max_raw, and nonnegative steps required")
    if axis not in (0, 1, 2):
        raise ValueError("layer_axis must be 0, 1 or 2")
    permutations = symmetry_permutations(model, symprec)
    raw_count = prod(g.multiplicity for g in model.groups)
    rng = np.random.default_rng(seed)
    if method == "enumerate":
        if raw_count > max_raw:
            raise ValueError(
                f"Enumeration requires {raw_count} raw assignments > max_raw={max_raw}; "
                "use icet-enumerate, sampling/SQS, or an explicitly larger bound"
            )
        iterator = product(*[list(_multisets(g.counts)) for g in model.groups])
    elif method in ("random", "clustered", "alternating", "layered"):
        # Precompute periodic same-group edges once, instead of rebuilding the
        # neighbor graph for every Monte Carlo trial.
        from pymatgen.core import Structure

        edges = []
        for group in model.groups:
            centers = Structure(
                model.structure.lattice,
                ["H"] * len(group.units),
                [model.structure[u[0]].frac_coords for u in group.units],
            )
            a, b, _, _ = centers.get_neighbor_list(cutoff)
            edges.append((a, b))

        def score_assignment(assignment):
            return sum(
                float(np.mean(np.array(labels)[a] == np.array(labels)[b]))
                for labels, (a, b) in zip(assignment, edges, strict=True)
                if len(a)
            )

        def patterns():
            for _ in range(samples):
                assignment = []
                for group in model.groups:
                    labels = [k for k, n in group.counts.items() for _ in range(n)]
                    if method == "layered":
                        shift = rng.random()
                        positions = [
                            (model.structure[u[0]].frac_coords[axis] + shift) % 1
                            for u in group.units
                        ]
                        order = np.argsort(positions, kind="stable")
                        decorated = [None] * len(labels)
                        for i, label in zip(order, labels, strict=True):
                            decorated[i] = label
                        labels = decorated
                    else:
                        rng.shuffle(labels)
                    assignment.append(labels)
                if method in ("clustered", "alternating"):
                    sign = 1 if method == "clustered" else -1
                    score = sign * score_assignment(assignment)
                    best_score, best_assignment = score, [list(a) for a in assignment]
                    active = [i for i, labels in enumerate(assignment) if len(set(labels)) > 1]
                    for step in range(steps if active else 0):
                        g = int(rng.choice(active))
                        i, j = rng.choice(len(assignment[g]), 2, replace=False)
                        if assignment[g][i] == assignment[g][j]:
                            continue
                        assignment[g][i], assignment[g][j] = assignment[g][j], assignment[g][i]
                        trial = sign * score_assignment(assignment)
                        temperature = 0.05 * (1 - step / max(steps, 1)) + 1e-6
                        if trial >= score or rng.random() < np.exp((trial - score) / temperature):
                            score = trial
                            if score > best_score:
                                best_score, best_assignment = score, [list(a) for a in assignment]
                        else:
                            assignment[g][i], assignment[g][j] = assignment[g][j], assignment[g][i]
                    assignment = best_assignment
                yield assignment

        iterator = patterns()
    else:
        raise ValueError(f"Unknown fixed-cell generation method: {method}")
    unique = {}
    for assignments in iterator:
        key = canonical_assignment(assignments, permutations)
        if key in unique:
            unique[key].provenance[0]["observed_multiplicity"] += 1
            continue
        metadata = {
            "method": method,
            "matrix": model.matrix,
            "seed": seed,
            "assignments": {
                g.name: list(a) for g, a in zip(model.groups, assignments, strict=True)
            },
            "observed_multiplicity": 1,
            "multiplicity_is_exact": method == "enumerate",
            "raw_assignment_space": raw_count,
            "symmetry_operations": len(permutations),
        }
        unique[key] = Candidate(
            model.decorate(assignments), [metadata], pair_features(model, assignments, cutoff)
        )
    return list(unique.values())


def reduce_candidates(candidates, *, stol=1e-4, ltol=1e-5, angle_tol=0.01):
    """Pymatgen matching retains chemical identities and permits primitive-cell matching."""
    from pymatgen.analysis.structure_matcher import StructureMatcher

    matcher = StructureMatcher(
        ltol=ltol,
        stol=stol,
        angle_tol=angle_tol,
        primitive_cell=True,
        scale=False,
        attempt_supercell=False,
    )
    unique, buckets = [], {}
    for candidate in candidates:
        key = candidate.structure.composition.fractional_composition.reduced_formula
        for index in buckets.get(key, []):
            if matcher.fit(unique[index].structure, candidate.structure):
                unique[index].provenance.extend(candidate.provenance)
                break
        else:
            buckets.setdefault(key, []).append(len(unique))
            unique.append(candidate)
    return unique


def select_candidates(candidates, count, *, method="diverse", seed=0):
    """Select within each composition using correlation diversity, never fake energies."""
    if count < 1 or method not in ("diverse", "random", "all"):
        raise ValueError("selection requires count > 0 and method diverse, random or all")
    buckets = {}
    for i, c in enumerate(candidates):
        key = tuple(sorted(c.structure.composition.fractional_composition.as_dict().items()))
        buckets.setdefault(key, []).append(i)
    selected = []
    rng = np.random.default_rng(seed)
    feature_keys = sorted({k for c in candidates for k in c.descriptors})
    for indices in buckets.values():
        if method == "all" or len(indices) <= count:
            selected.extend(indices)
            continue
        if method == "random":
            selected.extend(sorted(rng.choice(indices, count, replace=False).tolist()))
            continue
        x = np.array([[candidates[i].descriptors.get(k, 0) for k in feature_keys] for i in indices])
        if x.shape[1] == 0:
            raise ValueError("Diverse selection needs correlation descriptors")
        # Keep both edges of the sampled descriptor cloud, then fill its gaps.
        chosen = [int(np.argmax(np.linalg.norm(x - x.mean(axis=0), axis=1)))]
        while len(chosen) < count:
            distances = np.min(
                np.linalg.norm(x[:, None, :] - x[chosen][None, :, :], axis=2), axis=1
            )
            distances[chosen] = -1
            chosen.append(int(np.argmax(distances)))
        selected.extend(indices[i] for i in chosen)
    return sorted(selected)
