# Bonding analysis with CAT

`cat-bonding` reads completed LOBSTER calculations and exports CSV or Parquet tables for COHP, COBI, COOP, charges, DOS and electronic bands. It can classify contacts from chemistry or periodic structure topology, check requested interactions against the output, and compare integrated values with LobsterPy. [CATPlot](../../../../4_Plotting/README.md) turns the exported tables into figures.

Start with [one calculation](#start-with-one-calculation), then [plot the results](#plot-with-catplot). The remaining sections describe [inputs](#input-discovery), [configuration](#configuration-reference), [classification](#classifying-interactions), [outputs](#output-tables-and-conventions), [checks](#checking-the-results), [LobsterPy](#optional-lobsterpy-comparison), and the [Python API](#python-api-and-module-map).

## Install

Use Python 3.12 or newer. From the PhononTools repository root, in your chosen virtual environment:

```bash
python -m pip install -e '.[bonding]' -e ./4_Plotting
cat-bonding --help
catplot --help
```

The `bonding` extra installs Pymatgen for parsing and LobsterPy for the optional comparison. CATPlot is the local package in `4_Plotting`; no external plotting download is needed. Its Python import and primary command are both `catplot`.

For an existing uv setup, use `uv sync --extra bonding`, activate `.venv`, then `uv pip install -e ./4_Plotting`. A later `uv sync` may remove the separately installed plotting package; reinstall it if needed. CSV needs no extra I/O engine; CAT's dependencies already include PyArrow for Parquet.

These commands analyze existing output. Preparing and running VASP/LOBSTER is covered by the [LOBSTER input configuration guide](../../../../1_Tools/README.md#lobster-configuration-and-framework-detection).

## Start with one calculation

Keep the completed files together and put the analysis YAML beside the calculation:

```text
my-bonding-project/
├── calculation/        # COHPCAR.lobster, ICOHPLIST.lobster, POSCAR, etc.
├── bonding.yaml
└── plots/              # created when you copy/render the recipes below
```

Save this as `bonding.yaml`, adapting the chemistry rules to your system:

```yaml
system:
  path: calculation
  name: Na3PS4
  recursive: false

extraction:
  populations: [COHP, COBI, COOP]
  spin_mode: summed
  orbitals: true

interactions:
  element_groups:
    pn: [P, As, Sb]
    x: [S, Se]
    mobile: [Na]
  rules:
    - {name: Pn-X, groups: [pn, x]}
    - {name: Na-X, groups: [mobile, x]}
    - {name: X-X, groups: [x, x]}

benchmark:
  enabled: false

output:
  directory: bonding-analysis
  format: csv
```

Run from the project directory:

```bash
cat-bonding bonding.yaml
```

Paths in the YAML are relative to the YAML file, not your shell directory. Missing populations are skipped: a COHP-only calculation does not need dummy COBI/COOP files. A discovered file that cannot be parsed causes an error. Start by inspecting `bonding-analysis/lobster_runs.csv` and `bonding-analysis/bond_input_validation.csv`, then plot the extracted tables.

For several calculations, place them in separate subdirectories beneath `system.path` and set `recursive: true`. All tables use one `system.name`; `run` is the relative folder name, such as `phase-A` or `substitution-01`. A calculation directly in `system.path` has `run = .`. Comparison plots normally use `run` to distinguish calculations.

## Input discovery

Each detected run is read independently. Canonical filenames are preferred, followed by `.gz` variants and then sorted matching renamed/numbered files. Keep one calculation's outputs per directory: discovery selects the first matching file, not the newest or the most converged one.

| Files | What they enable |
| --- | --- |
| `COHPCAR.lobster`, `COBICAR.lobster`, `COOPCAR.lobster` | Energy-dependent population and integrated curves; each is optional |
| `ICOHPLIST.lobster`, `ICOBILIST.lobster`, `ICOOPLIST.lobster` | Integrated values, bond translations, multiplicities and interaction summaries |
| `POSCAR.lobster`, otherwise `CONTCAR`, otherwise `POSCAR` | Structure; needed for topology, DOS, bands and the LobsterPy comparison |
| `CHARGE.lobster` | Mulliken and Löwdin charges |
| `DOSCAR.lobster` plus structure | Total and projected DOS |
| `lobsterout` | LOBSTER version, warnings and spilling percentages |
| `bandOverlaps.lobster` | Maximum reported band-overlap deviation |
| `lobsterin` | Checks of requested populations and contacts against actual output |
| `lobsterin.bonding.json` or another `*.bonding.json` | Optional preparation metadata and topology settings |
| `FATBAND_*.lobster`, `KPOINTS`, structure, and a Fermi-energy source | Electronic band energies and projection weights |

Band extraction needs `vasprun.xml`, a discovered COXXCAR, or an explicit `bands.efermi_eV`. With no `vasprun.xml` and no supplied value, the first available COXXCAR supplies the Fermi energy. These must describe the same calculation/path. Discovery accepting `.gz` does not guarantee every downstream reader supports it: keep the structure and topology/input metadata uncompressed for this workflow.

## Configuration reference

The quick-start YAML is sufficient for basic extraction. The repository also has a [larger example](../../../examples/bonding.yaml); it enables the optional LobsterPy comparison, so review that choice before copying it.

| Setting | Default | Effect |
| --- | --- | --- |
| `system.path` | `.` | Root containing the LOBSTER calculations |
| `system.name` | Root folder name | Shared `system` label in all tables |
| `system.recursive` | `true` | Search child directories as well as the root |
| `extraction.populations` | `[COHP, COBI, COOP]` | Populations to extract; names are case-sensitive |
| `extraction.spin_mode` | `summed` | `summed` gives `spin: sum`; `resolved` gives available `up`/`down` channels |
| `extraction.orbitals` | `true` | Include available orbital-resolved population rows alongside totals |
| `extraction.energy_range_eV` | No restriction | Inclusive `[min, max]` filter for population curves and DOS; does not change ICOXXLIST values or band exports |
| `extraction.spilling_warning_percent` | `5.0` | Charge-spilling warning threshold in the run inventory; does not reject data |
| `bands.enabled` | `true` | Export FATBAND data when all prerequisites exist |
| `bands.efermi_eV` | Unset | Fermi-energy value passed to the band parser; useful without `vasprun.xml` |
| `benchmark.enabled` | `false` | Run the independent LobsterPy comparison |
| `output.directory` | `bonding-analysis` | Directory for tables and `manifest.yaml` |
| `output.format` | `csv` | `csv` or `parquet` |

`spin_mode` and `orbitals` apply to COHP/COBI/COOP. DOS and bands retain their available spin channels and projections independently. There are no separate YAML switches to disable charges or DOS when their files are present.

The supported command forms are:

```bash
cat-bonding bonding.yaml
cat-bonding analyze bonding.yaml
python -m cat_tools.bonding bonding.yaml
cat-bonding bonding.yaml --no-benchmark
cat-bonding benchmark.yaml --benchmark-only
```

`--benchmark-only` still writes inventory and input checks but skips ordinary analysis tables. Give it a separate output directory: rerunning an output directory replaces its generated tables and removes previously generated tables no longer in the result. Unrelated files are retained. `manifest.yaml` describes table keys, units, sources and checksums; it requires no manual maintenance.

The CLI also exposes `--out PATH` and `--format csv|parquet`. Currently, an override fails if the corresponding `output.directory` or `output.format` is also present in YAML. Edit the YAML instead, or omit that key when using its CLI flag. Relative `--out` paths resolve beside the configuration file too.

## Classifying interactions

### Chemistry rules

Rules are tested in order, independently of atom order; the first match wins. They can use element groups, as above, or literal pairs:

```yaml
interactions:
  rules:
    - {name: P-S, elements: [P, S]}
    - {name: framework-X, element_pairs: [[P, S], [Sb, Se]]}
```

`group_pairs` is also accepted for a list of group pairs. An unmatched contact uses its alphabetically ordered element pair, e.g. `Na-S`. Rules change the `interaction` label; the original atom labels and bond identity remain available.

### Periodic structure topology

Use topology when distinguishing contacts within and between coordination polyhedra or connected frameworks. It uses atom indices **and periodic translations**, not just atom-number ranges. Add `topology` under the existing `interactions` mapping:

```yaml
interactions:
  topology:
    framework_cations: [P, As, Sb]
    anions: [S, Se]
    mobile_cations: [Na]
    coordination: 4
    center_ligand_search_A: 3.4
    distance_tolerance_A: 0.03
```

This is an illustrative tetrahedral configuration, not a universal chemistry prescription. The selected structure must match LOBSTER's atom order. Use a `POSCAR` with an element-symbol line and a single scale factor. Topology settings are taken from explicit YAML first, then the optional preparation sidecar, then recognized element groups (`framework_cation`/`framework`, `anion`/`anions`, and optionally `mobile_cation`/`mobile`). Arbitrary group names such as `pn` and `x` do not automatically enable topology. `topology: {enabled: false}` explicitly disables it.

| Topology setting | Default / interpretation |
| --- | --- |
| `framework_cations`, `anions`, `mobile_cations` | Empty species lists; select the relevant chemistry |
| `center_atoms`, `ligand_atoms` | Empty lists; additional **one-based** atoms, unioned with species selections. Leave species lists empty to select only particular atoms |
| `coordination` | `null`: all ligands in the distance window. An integer selects that many nearest ligands and requires enough candidates |
| `center_ligand_min_A`, `center_ligand_search_A` | `0.0`, `3.4`: lower and upper center–ligand distances in Å |
| `distance_tolerance_A` | `0.03`: allowed difference between reconstructed and reported bond lengths |
| `framework_class`, `mobile_class` | `framework-ligand`, `mobile-ligand` |
| `intra_class`, `inter_class`, `unclassified_class` | `X-X_intra`, `X-X_inter`, `X-X_unclassified` |
| `inter_component` | Unset; optional connected-framework definition described below |
| `image_range` | `2`, retained for compatibility; current neighbor searches derive periodic bounds from the cell and distance cutoff |

Exact topology classification takes precedence over chemistry rules; rules still handle contacts outside the configured topology. Ordinary `X-X_intra` contacts share a polyhedron; ordinary `X-X_inter` contacts do not. This alone does not prove that they belong to disconnected molecules or chains.

For connected frameworks, add the following inside `interactions.topology`:

```yaml
inter_component:
  cutoff_A: 4.2
  selection: nearest
  tolerance_A: 0.03
```

This switches the inter-contact search to distinct periodic connected components. `selection: nearest` keeps nearest contacts with the specified distance tolerance; `all` keeps all eligible contacts up to the cutoff. Contacts within the same component that are neither framework bonds nor intra-polyhedron edges receive `X-X_same_component`. The cutoff controls the expected contact set for checks; it does not filter the population tables.

The current structure is always used to rebuild topology. A stale sidecar is reported through its structure hash; its chemistry settings may still seed that reconstruction. Supply explicit YAML if those settings no longer describe the system. Legacy `atom_index_groups` are accepted, but `split_by_atom_index_groups` no longer infers intra/inter contacts; it falls back to `unknown_group_name` or the rule name.

## Output tables and conventions

Each table has `system` and `run`. Filenames below use CSV; Parquet keeps the same names and columns. Only tables supported by the discovered inputs are written.

| Table | Rows and important columns |
| --- | --- |
| `lobster_runs.csv` | One row per run: input inventory, `status`, `warnings`, spilling percentages, `max_band_overlap_deviation` |
| `bond_integrals.csv` | One row per population/bond/orbital/spin: `integrated_value`, `length_A`, `multiplicity`, `interaction`, atom labels, translations |
| `bond_summary.csv` | One row per population/spin/interaction: `bond_count`, `contact_count`, `integrated_mean`, `integrated_sum`, `integrated_min`, `integrated_max`, `element_pairs` |
| `bond_curves.csv` | Population/bond/orbital/spin/energy rows: `energy_eV`, `value`, `integrated_value`, `series_kind`, bond metadata |
| `bond_curve_summary.csv` | Population/spin/interaction/energy rows: `value_mean`, `value_sum`, `integrated_mean`, `integrated_sum`, contact/bond counts |
| `charges.csv` | One row per site: `atom`, `element`, `mulliken_charge_e`, `loewdin_charge_e` |
| `electronic_dos.csv` | Energy/spin rows at several resolutions: `scope`, `site`, `element`, `orbital`, `DOS`, `integrated_DOS` |
| `electronic_bands.csv` | Band/spin/k-point rows: `distance`, `energy_eV`, `trace`, `point_label`, fractional `qx/qy/qz`, and `<atom>_<orbital>_weight` columns |
| `bond_input_validation.csv` | Always: requested populations versus discovered files, sidecar presence/hash checks |
| `bond_polyhedra.csv` | With topology: center–ligand membership, ligand images and distances |
| `bond_topology_contacts.csv` | With topology: exact intra-polyhedron X-X contacts and translations |
| `bond_request_validation.csv` | With topology and `lobsterin`: explicit contacts and generator counts checked against output |
| `bond_xx_validation.csv` | With topology and `lobsterin`: expected, requested and observed X-X contacts |
| `lobsterpy_benchmark.csv` | Optional: direct and LobsterPy element-pair means and their differences |

### Population values, signs and aggregation

The exported population values keep LOBSTER's sign. To display `−COHP` or `−ICOHP`, negate the chosen plotting column; do not negate COBI or COOP automatically. The exporter labels COHP/ICOHP values with `value_unit: eV` and COBI/COOP with `1`. Curve and DOS energies are exported on their LOBSTER relative-to-Fermi grids; bands explicitly subtract the parser's Fermi energy. Do not shift an exported `energy_eV` column a second time.

`bond_integrals.integrated_value` comes from ICOXXLIST at the Fermi level. `bond_curves.integrated_value` is the file's integrated curve at each sampled energy. Restricting the energy window clips rows; it does not reintegrate them.

Summaries use only total-bond rows, excluding orbital decompositions and the file's average curve. For values `v` and multiplicities `m`:

```text
contact_count   = sum(m)
integrated_sum  = sum(m * v)
integrated_mean = sum(m * v) / sum(m)
```

The curve summary applies the same weighting at each energy. `bond_count` counts distinct LOBSTER labels; `contact_count` includes multiplicity. Means compare per-contact values, while sums also depend on how many contacts were requested. Use comparable interaction requests when comparing different calculations.

In detailed curves, `series_kind` is `bond`, `orbital`, or `average`. Select the intended kind before plotting or aggregating. With no matching ICOXXLIST, curves fall back to multiplicity 1 and cannot recover exact translations from the list; classification can be rules-only or `unclassified`. Such summaries do not recover missing contact multiplicities.

### Sites, spin and projections

Bond-list `site1`/`site2`, topology atoms and LOBSTER labels are **one-based**. Charge/DOS `site`, band indices and `kpoint_index` are **zero-based**. DOS aggregate rows use `site = -1`.

DOS `scope` is `total`, `site_orbital`, `element_orbital`, or `element`; select one scope rather than summing them together. Only total DOS rows carry integrated DOS. DOS and bands use `up`/`down`, even when population extraction uses `spin: sum`. Spin-down DOS is exported without a plotting sign reversal.

Band `trace` identifies a spin/band pair within one run. Weight column names are sanitized atom/orbital names ending in `_weight`; inspect the CSV header to choose one. This exports electronic FATBAND data already present in the calculation; it does not generate FATBAND files or connect this analysis to FreqBall.

## Checking the results

`lobster_runs.status` reports LOBSTER warnings or charge spilling above the configured threshold. An `ok` row alone does not establish a trustworthy calculation: unavailable metrics remain empty, and the threshold does not judge basis quality for every system. Inspect the original warnings, spilling values and band overlaps alongside the data.

The input check reports `unavailable` without `lobsterin`, `warning` for missing or unexpected population outputs, and `ok` when those checks agree. It checks the actual input's population requests even if extraction selects fewer populations.

With topology, `topology_distance_ok`, `length_reconstructed_A` and `length_error_A` help detect structure/order/image mismatches. The request and X-X tables distinguish:

| Status or symptom | What to inspect |
| --- | --- |
| `missing_request` | A contact expected from the configured topology was not requested in `lobsterin` |
| `requested_missing_output`, `requested_contact_missing_output`, `output_missing` | Missing/incomplete lists, skipped populations or an unfinished calculation |
| `wrong_periodic_image_requested`, `wrong_periodic_image_in_output` | Atom pair matches but its periodic image does not; inspect translations and atom order |
| `generator_contact_count_mismatch` | Distance window and structural contact count versus summed output multiplicities |
| `unexpected_request`, `present_but_not_declared`, `unexpected_output_contact` | Input, configured topology and output describe different contact sets |
| Sidecar hash mismatch | Structure changed since input preparation; inspect the current settings and geometry |

Some statuses combine a wrong-image request with missing output. These are diagnostic tables; they do not repair inputs, delete contacts or rerun calculations.

Optional check settings are:

```yaml
validation:
  inter_xx_cutoff_A: 4.2
  validate_inter_xx: true
```

`validate_inter_xx: false` removes inter contacts from the expected set; it does not disable all checks. The cutoff falls back through explicit topology settings, matching sidecar settings and finally 4.2 Å when omitted. Set it explicitly when comparing a specific requested contact window.

## Optional LobsterPy comparison

This compares two analysis routes; it is **not a runtime/resource benchmark**. It requires a structure and both the curve and integrated-list file for a population.

```yaml
benchmark:
  enabled: true
  populations: [COHP]
  cutoff: 0.0
  noise_cutoff: 1.0e-8
  which_bonds: all
  type_charge: Mulliken
  absolute_tolerance: 1.0e-5
  relative_tolerance: 1.0e-4
```

These are the defaults except for `enabled`; populations default to those selected for extraction. Keep benchmark populations within that selection. Without CHARGE, the code uses LobsterPy's `Valences` mode. The direct mean sums resolved spins and takes an unweighted mean over bond labels by element pair; this differs from `bond_summary`'s multiplicity weighting. LobsterPy may also select a different set of bonds through its cutoffs and chemistry analysis.

The table contains `canonical_mean`, `lobsterpy_mean`, signed `difference = lobsterpy_mean - canonical_mean`, `relative_difference`, and `match`/`mismatch` from the supplied tolerances. Inspect selection and multiplicity before interpreting a mismatch as a parser error. Canonical plots use the ordinary population tables; this comparison never replaces their values.

## Plot with CATPlot

The [bonding recipes](../../../recipes/bonding) are editable YAML files. With the project layout from the quick start, run the following **from your project directory**, replacing `/path/to/PhononTools` with your checkout:

```bash
mkdir -p plots/recipes
cp /path/to/PhononTools/3_Analysis/recipes/bonding/*.plot.yaml plots/recipes/
catplot render plots/recipes/cohp-energy.plot.yaml
catplot render plots/recipes/icohp-substitution.plot.yaml
```

Recipe data paths resolve relative to the recipe file. With this layout, `../../bonding-analysis/...` points to the exported tables. Draft figures go to `plots/draft/`. Edit the recipe's data paths if your layout or table format differs.

| Recipe | Table and selection | Plot |
| --- | --- | --- |
| [cohp-energy](../../../recipes/bonding/cohp-energy.plot.yaml) | `bond_curve_summary`, COHP, one run, `spin: sum` | Mean −COHP versus energy, grouped by interaction |
| [icohp-substitution](../../../recipes/bonding/icohp-substitution.plot.yaml) | `bond_summary`, COHP, `spin: sum` | Mean −ICOHP per contact versus run, grouped by interaction |
| [bond-length](../../../recipes/bonding/bond-length.plot.yaml) | `bond_integrals`, COHP, total orbitals, one run, `spin: sum` | −ICOHP versus bond length |
| [dos](../../../recipes/bonding/dos.plot.yaml) | `electronic_dos`, element totals, one run, `spin: up` | Element-resolved DOS |
| [charges](../../../recipes/bonding/charges.plot.yaml) | `charges`, one run | Mulliken and Löwdin charge versus site |
| [fatbands](../../../recipes/bonding/fatbands.plot.yaml) | `electronic_bands`, one run, `spin: up` | Band dispersion; optional projection-weighted points |

The single-run recipes initially select `run: '.'`. For a recursive scan, replace that filter value with the actual subfolder name. The DOS and band recipes select `up`; add a separate layer or panel for `down` if present. For resolved population extraction, replace the `sum` filter with the channel you want.

Do not render every recipe blindly: charges, DOS or bands may be absent from your calculation. An empty plot usually means a filter selected no rows; a missing-column error usually means the recipe points to the wrong table or projection column.

### Change the quantity or inspect a particular bond

The COHP recipe uses a display transform rather than changing the source CSV:

```yaml
transforms:
  - {op: scale, column: value_mean, factor: -1, output: minus_COHP}
```

For COBI/COOP, change the population filter, remove the negation, map `y` directly to `value_mean`, and update the axis label. To show interaction sums rather than per-contact means, use `value_sum` or `integrated_sum` and label that choice.

For a particular bond, use `bond_curves.csv`, filter `population`, `run`, `spin`, `series_kind: bond`, `orbital: total`, and the chosen `bond_label`; plot `value` against `energy_eV`. CSV readers normally infer purely numeric labels as numbers, so filter `bond_label` with `value: 1`, for example. For an orbital decomposition, select `series_kind: orbital` instead and group by `orbital`.

CATPlot groups by one column and draws the supplied rows. It does not automatically separate runs, spins or bonds, or average duplicate bars. Filter to a single intended series, use separate panels/layers, or create a combined grouping column explicitly when comparing several dimensions. The supplied COHP filters prevent curves from different runs being joined together.

The band recipe groups by `trace` and preserves k-point order. Check path breaks before publication: the exporter supplies `point_label` but no branch column, so disconnected segments need separate layers or filtered data. Its commented scatter layer shows how to map an actual `_weight` column to point size; replace the example column with one in your output.

### Figure controls and Python use

```bash
catplot inspect line
catplot inspect scatter
catplot render plots/recipes/cohp-energy.plot.yaml \
  --set 'panels.0.axes.x.limits=[-8, 4]'
```

`style: {profile: paper, theme: default}` selects the default appearance; `layout: {type: single, edge_mm: 85}` sets panel size. Use the [CATPlot tutorial](../../../../4_Plotting/TUTORIAL.md) for stacked/grid layouts, legends, additional formats and annotations. Normal rendering writes draft figures. The optional `--final` mode additionally saves the recipe and a lock file with input hashes; it is not required for bonding analysis or routine plotting.

To render from Python:

```python
import matplotlib.pyplot as plt
from catplot.recipe import load_recipe
from catplot.render import render_to_files

path = "plots/recipes/cohp-energy.plot.yaml"
result = render_to_files(load_recipe(path), path)
plt.close(result.rendered.figure)
```

## Python API and module map

Run the same workflow from a script:

```python
from cat_tools.bonding import run_workflow

written_paths = run_workflow("bonding.yaml", no_benchmark=True)
```

Or extract a table directly without writing it:

```python
from cat_tools.bonding import InteractionClassifier, discover_runs
from cat_tools.bonding.extract import integral_artifact, summary_artifact

runs = discover_runs("calculation", system="sample", recursive=False)
classifier = InteractionClassifier.from_config({})
integrals = integral_artifact(runs[0], "COHP", classifier, orbitals=False)
summary = summary_artifact(integrals.data)
print(summary.data)
```

Extractor functions return a `TableArtifact`: `.data` is the Pandas DataFrame, with `.name`, `.keys`, `.units` and `.metadata` describing it. Call `.validate()` before consuming a manually assembled artifact; `run_workflow` validates during export.

| Module | Responsibility / callable interface |
| --- | --- |
| [discovery.py](discovery.py) | `discover_runs`: discover available inputs without parsing them |
| [model.py](model.py) | `PopulationFiles` and `LobsterFiles`: paths, run identity and completeness properties |
| [classify.py](classify.py) | `InteractionClassifier.from_config`, `.classify`, `.classify_details`, `.with_topology`; element-pair helpers |
| [extract.py](extract.py) | `integral_artifact`, `curve_artifact`, both summary functions, `charge_artifact`, `dos_artifact`, `band_artifact`, `run_artifact` |
| [validation.py](validation.py) | `ValidationSettings` and the input/request/X-X/polyhedra/contact artifact builders |
| [benchmark.py](benchmark.py) | `benchmark_artifact`: independent LobsterPy comparison |
| [workflow.py](workflow.py) | `run_workflow` and the CLI; discovery, extraction, aggregation and export |
| [topology.py](topology.py) | Compatibility exports of `cat_shared.bonding.topology`, including `TopologyConfig`, `build_polyhedra`, contact geometry and manifest helpers |
| [lobsterin.py](lobsterin.py) | Compatibility exports of the shared input parser and `cat_prepare.lobsterin` generator |
| [__init__.py](__init__.py), [__main__.py](__main__.py) | Public convenience imports and `python -m cat_tools.bonding` |

`parse_lobsterin` reads calculation intent. Input generation remains available via `cat-prepare lobsterin --config preparation.yaml -o lobsterin`; the older `cat-bonding lobsterin ...` and `cat-bonding generate-lobsterin ...` forms delegate to it. Preparation YAML and analysis YAML describe different tasks: follow the [preparation guide](../../../../1_Tools/README.md#lobster-configuration-and-framework-detection) for basis sets, interaction requests, topology derivation and optional sidecars.
