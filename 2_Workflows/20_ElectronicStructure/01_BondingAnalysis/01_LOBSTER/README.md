# VASP → LOBSTER → bonding tables and figures

**Use:** Prepare a basis-compatible VASP calculation, run LOBSTER, inspect projection quality and export COHP/COBI/COOP, charges, DOS and electronic fatbands.

**Files:** [lobsterin.seed.yaml](lobsterin.seed.yaml) · [bonding.yaml](bonding.yaml) · [INCAR.lobster.template](INCAR.lobster.template) · [full analysis reference](../../../../3_Analysis/src/cat_tools/bonding/README.md).

## 1. Choose a projection basis before the wavefunction run

Start from the accepted relaxed `POSCAR` and its matching `POTCAR`. Install CAT's `bonding` extra; LOBSTER itself is a separately installed executable with its own version/manual. Edit the seed for every species, valence/semicore orbital set and interaction range in your material. Basis functions must be compatible with the PAW datasets and the chosen LOBSTER basis library; a familiar basis-set name is not evidence of compatibility for every potential or functional.

```bash
cat-prepare lobsterin POSCAR --config lobsterin.seed.yaml -o lobsterin
```

`--config` selects the YAML containing the basis and interaction requests. Keep `lobsterin.bonding.json` beside the generated file. It records requested contacts, structure identity and topology. The Si/O seed is illustrative; material-specific examples are in [CAT examples](../../../../1_Tools/CAT/examples/).

`Generator` creates type-pair distance ranges across all selected species; explicit `Atom` contacts and framework selections preserve periodic cell images. Use `Generator: []` for strictly site-specific requests. The framework model supports center–ligand bonds, intra-polyhedron contacts and inter-component contacts; see [periodic topology preparation](../../../../1_Tools/README.md#inter-polyanion-and-inter-chain-contacts). This replaces the archived PairFinder invocation and avoids losing image information. LOBSTER atom numbers in input contacts are one-based; CAT output tables declare their own indexing in the manifest.

## 2. Prepare the static `WAVECAR` calculation

Adapt the linked `INCAR` template using your converged functional, cutoff, k mesh, magnetic/SOC treatment and electronic tolerance. `NSW=0`, `IBRION=-1` keep the geometry fixed; `LWAVE=.TRUE.` retains the wavefunctions. For the ordinary collinear route use the compatible `ISYM=0` or `-1` choice for the selected LOBSTER setup. `LORBIT=11` is useful for VASP orbital projections, but is not what creates LOBSTER's projection and is not a universal LOBSTER prerequisite.

Derive the basis band count with the installed Pymatgen instead of copying the legacy `NBANDS=64`:

```python
from pymatgen.io.lobster import Lobsterin
basis = Lobsterin.from_file('lobsterin')
basis.write_INCAR(incar_input='INCAR.lobster.template', incar_output='INCAR',
                  poscar_input='POSCAR', isym=0)
```

Inspect the resulting `INCAR`. This utility sets a small number of tags and `NBANDS`; it does not validate the full physical model. Check VASP's **actual** `NBANDS` after parallelization-related rounding and inspect LOBSTER's handling of available bands. Benchmark parallelization using the [common procedure](../../../00_Common/00_SettingUp/VASP/benchmarking.md); there is no blanket prohibition on band parallelization for every VASP/LOBSTER setup.

For ordinary bonding/DOS use a converged integration mesh. Execute VASP with your normal scheduler script. Require electronic convergence and preserve matching `POSCAR`, `POTCAR` identity, `INCAR`, `KPOINTS`, `WAVECAR`, `OUTCAR` and `vasprun.xml`. Do not substitute a `WAVECAR` from a different cell, basis count, potential set or k-point list.

## 3. Run LOBSTER and assess projection quality

Run the installed LOBSTER executable in that directory, for example:

```bash
lobster > lobster.stdout 2>&1
```

Use an allocation and thread count appropriate to that executable. Inspect native `lobsterout`, charge/total spilling, warnings and `bandOverlaps.lobster` when produced. Compare reconstructed DOS/bands and vary plausible bases where needed. Low spilling is necessary evidence of a useful projection, not a proof that every requested bond descriptor is converged. Record the basis, software version and actual input.

Projection reuse (`saveProjectionToFile` / `loadProjectionFromFile`, when supported) is only valid for matching wavefunctions and basis settings. Do not reuse a saved projection across changed geometry, potentials or k meshes. CAT's seed exposes `system.save` for saving projections and `system.write` for writing basis functions; advanced native keywords can be added to the final lobsterin according to the installed manual. Retain those edits with its preparation sidecar; regenerating the file will replace manual additions.

## 4. Extract the supported canonical data

Set `system.path` in `bonding.yaml` to the completed run, then:

```bash
cat-bonding bonding.yaml
```

| Native data | CAT result / check |
| --- | --- |
| ICOHPLIST / ICOBILIST / ICOOPLIST | Per-contact integrals, classifications and summaries |
| COHPCAR / COBICAR / COOPCAR | Total/orbital, summed/resolved-spin curves and window summaries |
| CHARGE | Mulliken/Löwdin charge table |
| `DOSCAR.lobster` | Total, site/orbital and element electronic DOS |
| FATBAND files with matching `KPOINTS`/structure/Fermi reference | Electronic band energies and projection columns |
| lobsterin, sidecar and actual geometry | Requested-vs-found contact and topology validation |
| lobsterout / bandOverlaps | Run inventory and projection-quality diagnostics |

Missing optional output is not manufactured. Inspect the inventory and validation tables before interpreting a missing interaction as zero. Configure `spin_mode: resolved` when comparing spin channels; `summed` combines them for population curves/integrals. DOS and FATBAND retain their available spin channels. A curve energy-window filter does not recalculate native ICOXXLIST integrals. Use consistent energy reference and sign: a plot of −COHP reverses the native COHP sign ; it does not change the stored data.

The optional LobsterPy comparison is diagnostic and does not generate CAT's primary tables. Enable it in a separate analysis configuration/output directory when useful. Thresholds such as 5% spilling are configurable diagnostic choices, not universal scientific acceptance criteria. Explicit topology settings take precedence over advisory sidecar settings; inspect stale-structure warnings.

## 5. Plot and compare

Copy the appropriate [bonding recipes](../../../../3_Analysis/recipes/bonding/), edit table/output paths and filters, then:

```bash
catplot render cohp-energy.plot.yaml
catplot render charges.plot.yaml
catplot render dos.plot.yaml
catplot render fatbands.plot.yaml
```

CATPlot is a separate optional installation. Select one run/spin where appropriate; for electronic fatbands replace the recipe's example orbital-weight column with an actual exported column. Keep identical energy windows, sign conventions and suitable shared scales when comparing compositions. Bonding integrals alone do not establish phase stability, ionic mobility or an activation barrier.

## 6. Electronic FATBAND preparation

A normal mesh-only `WAVECAR` does not contain an arbitrary plotted path. Prepare the combined integration mesh plus zero-weight path points expected by the installed LOBSTER workflow. Pymatgen provides `Lobsterin.write_KPOINTS` with `line_mode=True`; use the corresponding primitive cell and inspect the path/labels before running. Add native `createFatband` requests for the desired basis orbitals, generate the matching VASP wavefunctions and rerun LOBSTER. Large all-orbital requests can be expensive; start with the chemically relevant projections.

Preserve the exact `KPOINTS` and a compatible Fermi energy from the calculation. CAT's band reader requires FATBAND files, structure and `KPOINTS`; the Fermi energy comes from matching `vasprun.xml`, a supported population file or an explicit `bands.efermi_eV`. A pure line-mode VASP calculation described in the separate [bands guide](../../02_BandsAndDOS/README.md) is not automatically this LOBSTER input.

## 7. Advanced native outputs retained from the legacy scope

The legacy guide also listed analyses beyond CAT's current extractors:

| Native analysis | How to continue |
| --- | --- |
| k-point/bond-resolved spilling, k-space COHP | Request through the installed LOBSTER manual; preserve diagnostics alongside the standard quality tables |
| Gross populations, Madelung energies/site potentials | Keep native files and their charge/basis convention; CAT currently exports the supported charge table, not all these quantities |
| BWDF / COHP-weighted distributions | Use native output for disordered/extended systems; no current CAT distribution-table exporter |
| Density-of-energy decomposition | Native LOBSTER calculation and analysis; not interchangeable with pair COHP |
| Real-space wavefunctions/orbitals, Hamiltonian/overlap matrices | Native version-specific requests and external visualization/downstream model; not a CAT tight-binding workflow |
| Molecular/fragment orbitals and formation-energy analysis | Explicit native fragment definitions including periodic images; not inferred from CAT's polyhedron labels |

These remain available through LOBSTER where supported. Do not put arbitrary manual keywords into CAT's strict seed schema or assume its parser exposes every native feature. Preserve the archived [legacy keyword discussion](../../../../OLD_Workflows/05_Lobster/README.md) for provenance, but verify spelling, availability and scientific meaning against the manual shipped with your executable; that historical list is not a validated input.

Sources: [Pymatgen LOBSTER input utilities](https://pymatgen.org/pymatgen.io.lobster.html), [VASP orbital projections](https://vasp.at/wiki/LORBIT), and your installed LOBSTER manual. The maintained procedure improves basis-derived band counts, periodic contacts, provenance, output validation and reproducible plotting over the legacy guide.
