# Ordering command and model reference

`cat-ordering` generates fixed-composition atomic substitutions, coupled local motifs, extended ordered/clustered structures, and genuine icet SQS candidates. It writes structures and selects candidates for subsequent calculations. It does not run a DFT code. The same interface is available as `cat-tools ordering` and `cat-prepare ordering`.

Install from the repository root:

```bash
uv sync --extra ordering
# or: python -m pip install -e '.[ordering]'
cat-ordering --help
```

The implementation uses pymatgen for structures/matching and full phase diagrams, spglib for parent symmetry, icet/ASE for supercell enumeration and cluster-based SQS, and SciPy/HiGHS for composition-balanced mixing hulls. See also the [inventory of the JCH reference material](JCH_INVENTORY.md).

For a numbered procedure using your own inputs and working directory, follow the [ordering workflow](../../../2_Workflows/00_Common/05_Ordering/README.md).

An [editable system-agnostic template](../examples/ordering/template.yaml) is available for defining a new study; replace its element placeholders before use.

## Start with the supplied examples

From the repository root:

```bash
cat-ordering generate 1_Tools/CAT/examples/ordering/atomic.yaml --output /tmp/cat-atomic
cat-ordering generate 1_Tools/CAT/examples/ordering/motifs.yaml --output /tmp/cat-motifs
cat-ordering generate 1_Tools/CAT/examples/ordering/icet.yaml --output /tmp/cat-icet
cat-ordering generate 1_Tools/CAT/examples/ordering/sqs.yaml --output /tmp/cat-sqs
cat-ordering hull 1_Tools/CAT/examples/ordering/hull-toy.yaml --output /tmp/cat-hull
```

Choose fresh output directories: existing paths are never overwritten. YAML paths are relative to the YAML file; command-line output paths are relative to the shell. The supplied cubic Na3PS4 CIF is an idealized geometry generated with space group 217, a=7 Å and S coordinate x=0.18. It is a software example, **not a relaxed DFT reference**. Replace it with a validated scaffold for scientific work. Repeat the workflow for each relevant polymorph, using distinct `id_prefix` values.

## Specify occupational degrees of freedom

Every sublattice selects parent sites by `elements: [P, Sb]` or zero-based `indices: [6, 7]`. Selections repeat into each requested supercell. Atomic `states: [P, Sb]` permit independent substitutions. `fractions` specifies exact fractions of that selected sublattice, not fractions of all atoms. Fractions must sum to one and be commensurate with each requested cell; CAT never rounds them. The `compositions` list can override named sublattice fractions to scan a grid.

For correlated whole-unit substitutions, specify a motif:

```yaml
sublattices:
  - name: tetrahedra
    elements: [P]
    motif:
      neighbor_elements: [S]
      coordination: 4
      cutoff_A: 2.7
    states:
      PSe4: {center: P, neighbor: Se}
      SbS4: {center: Sb, neighbor: S}
    fractions: {PSe4: 0.5, SbS4: 0.5}
```

This binds each center to its four periodic neighbors, including tetrahedra that cross cell boundaries. A PSe4 unit gets a P center and four Se ligands together. The motif states can be any elements and any coordination number. All ligands within each supported motif state share one species; independent atomic sublattices cover mixed-ligand configurations. Overlapping selections, shared ligands, missing neighbors and ambiguous periodic memberships are rejected. Shared-polyhedron networks require a different constraint model and are not silently treated as independent motifs. Ordered parent scaffolds and element substitutions are supported; vacancies and positional/orientational disorder are not implemented in this first interface.

## Choose generators and cells

`supercells` accepts diagonal triples or full integer 3×3 matrices with positive determinant. Cell vectors and atomic positions are inherited without relaxation.

| Method | Meaning |
| --- | --- |
| `enumerate` | Exhaustive multiset decorations at the exact specified compositions in each explicit cell; reduced using valid parent symmetry operations. |
| `icet-enumerate` | icet enumeration of symmetry-inequivalent decorations and supercell shapes at `icet_sizes`, measured relative to the supplied input cell. |
| `random` | Reproducible fixed-composition random samples; these are not SQS. |
| `alternating` | Swap search favoring unlike neighboring units within each variable group. |
| `clustered` | Swap search favoring like neighboring units within each variable group. |
| `layered` | Contiguous occupation blocks ordered along fractional `layer_axis` (0, 1 or 2), with sampled periodic origins. |
| `sqs` | icet simulated-annealing optimization of cluster correlations against the ideal random target. |

`alternating`/`clustered` are heuristic pattern searches, not energy minimization; frustrated lattices may not permit perfect alternation. For PS4/SbSe4-type motifs, these searches act on motif-center neighbors. For independently mixed atomic groups, they bias each group separately; the output descriptors also measure cross-group correlations. `pair_cutoff_A` determines the scored neighbor range, including periodic-image multiplicities. Select a meaningful shell: if the cutoff contains no variable-center neighbors, the bias cannot distinguish patterns.

`samples`, `pattern_steps`, `seed`, `layer_axis`, `symprec`, `max_raw`, `max_candidates` and `max_enumerated` make search limits explicit. Exceeding an enumeration limit raises an error instead of returning an apparently complete truncated list. Runtime still grows rapidly; use small exact cells, icet shape enumeration, and larger sampled/SQS cells together. A finite pool cannot certify the true ground state across all possible supercells.

## SQS and prescribed short-range correlations

```yaml
methods: [sqs]
supercells: [[2, 2, 2]]
sqs:
  cutoffs_A: [7.5, 4.0]
  steps: 10000
  repeats: 3
```

Cutoffs specify successive cluster orders (pairs, triplets, ...). Each repeat uses a recorded seed. The manifest saves the cluster space, achieved and ideal target vectors, maximum absolute mismatch, and icet's objective (smaller is better). Finite-cell SQS generally do not match all target correlations exactly. Inspect errors and converge cell size, cutoffs and search effort rather than assuming the word SQS guarantees quality. Optimize an atomic SQS for an independently random atomic alloy; optimize a motif SQS for random placement of intact correlated units. Those are different statistical ensembles.

Internally, categorical state labels are mapped to unique element-like proxy symbols for icet. They are decoded back before any physical structure is written; the encoding is saved in the metadata. The full scaffold remains present so its geometry constrains symmetry. No pseudopotentials or physical interpretation are attached to these internal symbols. The adapter uses icet's target-vector helper and is constrained to icet 3.2–4.x; integration tests exercise icet 4.0.

For specified short-range order, `sqs.target_vector` can replace the random target. Copy the ordering/length from the saved cluster space and provide the complete vector. Such outputs are labeled `target-correlations`, not SQS. A target may be unrealizable; the saved residuals reveal that limitation. Pure endpoints should be generated with `enumerate`, not an SQS search without any swaps.

## Symmetry reduction, selection and provenance

Exact fixed-cell enumeration uses spglib operations that preserve fixed atoms, sublattice identities and whole-motif memberships. Its within-cell orbit counts are exact. Sampling hit counts are explicitly **not** thermodynamic degeneracies. Pymatgen StructureMatcher then removes chemically equivalent candidates across methods and cells, permitting primitive-cell reduction but no lattice rescaling. Reported generated symmetry uses the configured `symprec` (default 1e-5 Å) and `angle_tolerance=-1` (spglib automatic), consistent with enumeration. Collected relaxed symmetry uses 1e-3 Å and the same automatic angle criterion; both outputs record these tolerances. Recheck relaxed symmetry across tolerances before assigning a phase. Pymatgen's default 5-degree angle tolerance can otherwise promote a tetragonally strained cell to cubic. Tolerances are intentionally tight for generated, unrelaxed structures. They do not constitute a method for identifying thermally distorted structures.

`selection: {method: diverse, per_composition: 12}` uses farthest-point sampling in normalized periodic pair-fraction descriptors, independently at each chemical composition. `random` and `all` are also supported. This offers different local arrangements for calculation; it neither predicts energies nor guarantees that every hull ground state was selected. For small pools, calculate all candidates.

Output includes:

- `manifest.json`: resolved source/hash, input settings, package versions, candidate IDs, structures, space groups, descriptors and generation provenance.
- `candidates/ID/structure.json`, `POSCAR`, `metadata.json`: physical structures for each unique candidate. JSON preserves original site order; `POSCAR` is grouped by species. `poscar_to_structure_site` maps `POSCAR` indices to JSON indices.
- `candidates.csv`, `selected.txt`: inspect the pool and chosen calculation inputs.
- `energies.template.csv`: explicit cell-energy/convergence/method fields for selected candidates; initially blank because generation does not supply energies.

Changed species require appropriate `POTCAR` assembly and consistent numerical settings for every calculation. Relax each selected candidate and follow with an ordinary static calculation on the relaxed cell. Track final symmetry instead of trusting the original phase label: different parents may relax to the same phase.

## Collect energies and construct the mixing hull

Put results in `CALCULATIONS/ID/vasprun.xml`, then:

```bash
cat-ordering collect /tmp/cat-atomic/manifest.json CALCULATIONS \
  --energy-model PBEsol-520-validated-kdensity --output atomic-energies.json
cat-ordering hull hull-calculated.yaml --output mixing-hull
```

The collector requires converged ordinary `NSW=0, IBRION=-1` runs, unchanged composition, and matching supported global electronic settings/potential titles. It records final geometry/symmetry and actual k-points. Potential titles do not prove identical potential files, and identical numerical meshes are not appropriate for arbitrary cell sizes: potential identity, k-point density, smearing choice and energy convergence still need scientific validation. DFT+U collection is rejected until a per-species compatibility policy is supplied; audited entries can be provided directly instead. All selected calculations are required; missing files fail instead of silently changing the hull dataset.

Alternatively, fill the CSV template and give `manifest: .../manifest.json` plus `energies: .../energies.csv`, or provide inline `entries`. Each entry needs a unique `id`, the **actual cell composition**, `energy_eV` for that same cell, a shared `energy_model` identifier, and `converged: true`. JSON collections already contain their compositions. `energies` may be a list of collections from different phases; use distinct ID prefixes. Never pair a reduced formula with a supercell total energy.

The mixing-hull configuration names `terminal_ids` for the chosen endmembers. SciPy/HiGHS minimizes the energy over nonnegative, composition-balanced mixtures. The four Na3PnX4 endmembers span two independent substitution coordinates and are linearly dependent; the solver handles this directly. Entries outside the terminal composition space are rejected.

`hull.json`/`hull.csv` report energy per atom, energy above the supplied hull, on-hull status and decomposition weights. Weights are **atomic fractions**, not formula-unit fractions. `mixing_energy_eV_atom` uses the lowest-energy balanced mixture of the specified terminal entries as its reference. For a binary this is the endpoint tie line; for four dependent endpoints it is their convex envelope, not a bilinear interpolation. Several equally valid decompositions may exist.

The hull is a zero-temperature energy hull of the supplied configurations. An SQS is a representative random configuration, not necessarily a ground state. Vibrational/configurational entropy and finite-temperature free energies are not implicitly added. The synthetic `hull-toy.yaml` only demonstrates normalization and the two-dimensional composition geometry; it must not be used as materials data.

## Later expansion to a full thermodynamic hull

The same energy-entry format supports `mode: full`, backed by pymatgen PhaseDiagram. Add consistently referenced elemental and competing-compound entries, remove `terminal_ids`, and retain the same energy normalization. Missing elemental references are rejected. The completeness of the competing-phase dataset and any DFT compatibility/corrections remain explicit responsibilities. A mixing hull within Na3PnX4 alone cannot establish stability against other compounds.

## References and validation

- [icet structure enumeration](https://icet.materialsmodeling.org/advanced_topics/structure_enumeration.html)
- [icet SQS generation and objective](https://icet.materialsmodeling.org/advanced_topics/sqs_generation.html)
- [icet tools API](https://icet.materialsmodeling.org/moduleref/tools.html)
- [pymatgen phase diagrams](https://pymatgen.org/pymatgen.analysis)
- [Supercell paper, Okhotnikov et al. (2016)](https://doi.org/10.1186/s13321-016-0129-3)

`python -m unittest discover -s tests -p test_ordering.py -v` tests exact symmetry orbits, periodic motifs, stoichiometry, icet enumeration/SQS, reproducibility, export mapping, energy normalization, dependent endmembers and hull decompositions.

## Supercell as an optional enumeration backend

The [Supercell program](https://orex.github.io/supercell/) generates ordered occupations from disordered CIFs, including vacancies, charge balancing, symmetry-equivalent merging and output sampling. For a specified finite cell it can provide representative configurations and their symmetry multiplicities. It is useful as an independent check of small-cell enumeration or for importing experimental partial occupancies. Enumeration can grow combinatorially; writing only a sample does not guarantee that enumeration itself becomes inexpensive.

CAT currently uses icet/spglib/pymatgen. There is **no Supercell CLI adapter or CIF-pool importer in `cat-ordering` yet**. A future adapter should preserve cell, composition, symmetry settings, representative multiplicities and generator provenance in the existing manifest, then reuse selection and energy collection. Do not try `methods: [supercell]`: it is not an implemented method.

For Na3(P,Sb)(S,Se)4, conventional formal charges assign P/Sb +5 and S/Se -2. Consequently, at identical coordinates and lattice, these label substitutions leave every point charge unchanged: a formal-charge Ewald ranking cannot separate their ordering energies. Bonding and size effects still change the DFT energies after relaxation. Charge balancing/electrostatics can be more useful for heterovalent substitutions or vacancies; even then electrostatics is only a screening model.

Independent P/Sb and S/Se occupations do not enforce intact SbS4/PSe4 units. For that question, use CAT's coupled motif model. Supercell's random output sampling is not a correlation-optimized SQS; use the [icet SQS procedure](https://icet.materialsmodeling.org/advanced_topics/sqs_generation.html) for that ensemble. Neither enumeration nor SQS provides DFT energies or proves that the ground state has been found. Finite-cell multiplicities are not by themselves a complete finite-temperature entropy model.
