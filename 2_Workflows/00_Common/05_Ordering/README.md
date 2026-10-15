# Ordering and mixing hull

**Use:** Generate occupational configurations and build a mixing hull from converged energies.

**Files:** [Ordering template](../../../1_Tools/CAT/examples/ordering/template.yaml) · [convergence record](convergence-record.md) · [input examples](../../../1_Tools/CAT/examples/ordering/README.md).

## 1. Choose the inputs and working directory

Start wherever you want to keep this study, with `cat-ordering` available. You need an ordered parent structure and a definition of which sites may change occupation. Use your own filenames and directory layout; the short names below are examples. No project-path environment variable is required.

```bash
cat-ordering --help
```

If CAT is not installed, follow the [tools setup](../../../1_Tools/README.md#setup) with the `ordering` extra (`uv sync --extra ordering` from your checkout, or `python -m pip install -e '.[ordering]'` in your chosen environment).

The main sequence is **define occupations → generate → inspect/select → converge and calculate → collect → build and refine the hull**. Start with one parent; repeat for other parents or phases relevant to your question. Atomic mixing, coupled motifs and SQS are choices, not three mandatory stages.

Throughout this guide, replace `/path/to/...` and uppercase tokens such as `ID` before running commands. Relative shell paths refer to the directory where you run the command; structure paths inside YAML refer to that YAML's directory. Use the [tools reference](../../../1_Tools/CAT/ordering/README.md) for all model options and the [editable template](../../../1_Tools/CAT/examples/ordering/template.yaml) for a starting configuration. Concrete material examples and the [JCH inventory](../../../1_Tools/CAT/ordering/JCH_INVENTORY.md) remain in tools.

## 2. Inspect the parent structure

Use an existing ordered structure, or convert one into a format CAT can read. For example, from your chosen working directory:

```bash
cat-prepare structure convert /path/to/starting.cif parent.vasp
```

Skip conversion if your existing file is suitable; point the YAML directly to it. Resolve partial occupancies before writing an ordered VASP structure. Inspect species, coordinates, contacts, connectivity and the cell in your viewer.

For a numerical summary, supply your actual structure path after `python -`:

```bash
python - parent.vasp <<'PY'
import sys
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
s = Structure.from_file(sys.argv[1])
print('Composition:', s.composition, 'Sites:', len(s), 'Ordered:', s.is_ordered)
print('Lengths [A]:', s.lattice.abc, 'Angles [deg]:', s.lattice.angles)
for tolerance in (1e-5, 1e-3, 1e-2):
    sym = SpacegroupAnalyzer(s, symprec=tolerance, angle_tolerance=-1)
    print('symprec [A]:', tolerance, sym.get_space_group_symbol(), sym.get_space_group_number())
PY
```

Check tolerance-sensitive symmetry assignments against the geometry. The explicit angle setting uses spglib's automatic criterion; pymatgen's default 5-degree angular tolerance can promote a tetragonally strained cell to cubic. A directory or candidate label is not evidence of a phase. Keep distinct parent structures only when they represent intended alternatives.

## 3. Define the occupational model

Copy the linked template to your working directory, or write your own YAML:

```bash
cp /path/to/PhononTools/1_Tools/CAT/examples/ordering/template.yaml ordering.yaml
```

Edit `ordering.yaml` with your preferred editor. The template describes one binary sublattice solely to show the input structure; it is **not ready to run** until you replace its element placeholders and review its settings.

1. Set `structure` to your parent file and give the batch a unique `id_prefix`.
2. Select variable parent sites by real `elements` or zero-based `indices`. All other sites remain fixed in species during generation.
3. Set the allowed `states` and exact sublattice fractions. Add sublattices if your problem has independent occupations on more than one set of sites.
4. Define the `compositions` to sample, including the terminal compositions needed for your intended mixing space. A binary endpoint pair is only one possibility.
5. Choose `supercells` so every requested occupation has an integer number of sites. The generator does not round fractions. Review enumeration cost as the cell and number of variable sites grow.
6. Start with `methods: [enumerate]` and `selection: {method: all}` for a manageable exact pool. Set finite enumeration limits and a meaningful neighbor cutoff.

For **coupled polyhedral units**, replace independent atomic states with the [motif model](../../../1_Tools/CAT/ordering/README.md#specify-occupational-degrees-of-freedom). Choose centers, ligand selection, coordination and cutoff from your structure; verify periodic membership and absence of shared ligands. Independent label swaps do not impose a whole-unit constraint. Supported motif states use one ligand species per unit; arbitrary overlapping constraints need a different model.

If needed, prepare additional YAML files for other parents, supercell shapes, clustered/layered patterns or SQS. Use new ID prefixes for separate batches. The [generator reference](../../../1_Tools/CAT/ordering/README.md#choose-generators-and-cells) explains these options. No particular crystal system, number of phases or chemical family is required by this procedure.

## 4. Generate and inspect the configurations

From the directory containing `ordering.yaml`:

```bash
cat-ordering generate ordering.yaml --output pool-r01
cat pool-r01/candidates.csv
cat pool-r01/selected.txt
```

Use a fresh output name; generation refuses existing directories. It generates structures and metadata without running VASP. `candidates.csv` lists IDs, compositions, atom counts, symmetry and selection. `selected.txt` contains `POSCAR` paths **relative to the pool**, not bare IDs.

Inspect representative `pool-r01/candidates/ID/POSCAR` files in your viewer. Check stoichiometry, contacts, local arrangements and supercell size. For exact small pools, calculate all representatives if feasible. For larger pools, geometric diversity selection reduces cost but is not an energy ranking; record that approximation and revisit unselected candidates in step 11.

Reduction occurs within each generation call; different batches may overlap. Preserve manifests and IDs once calculations begin. An enumeration-limit error means the requested search was not completed; increasing the limit or changing to sampling is an explicit decision, not evidence of exhaustive coverage.

**Optional Supercell cross-check:** Supercell can independently enumerate and symmetry-merge a fixed-cell occupation problem, including partially occupied CIF inputs. CAT does not yet invoke it or import its output pool. Random sampling is not SQS optimization, and formal-charge ranking cannot distinguish isovalent swaps at unchanged coordinates. See the [Supercell reference](../../../1_Tools/CAT/ordering/README.md#supercell-as-an-optional-enumeration-backend).

## 5. Check SQS quality, if using SQS

Skip this step for a purely ordered search. For an SQS batch, inspect its manifest with the following command (replace the input path if needed):

```bash
python - pool-r01/manifest.json <<'PY'
import json, sys
for candidate in json.load(open(sys.argv[1]))['candidates']:
    for p in candidate['provenance']:
        if p['method'] == 'sqs':
            print(candidate['id'], 'seed:', p['seed'], 'steps:', p['steps'],
                  'objective:', p['objective'],
                  'max residual:', p['max_abs_correlation_error'])
PY
```

Increase optimization effort and independent repeats in a new YAML/output batch. Then compare larger cells, different shapes and expanded correlation cutoffs at the same composition. Inspect the achieved and target cluster vectors by shell; a single maximum residual does not establish quality. Compare objectives only within the same cluster space. Record residuals and any finite-cell limitations.

After calculation, test energies and relevant observables across sizes and seeds. An atomic SQS and a coupled-motif SQS represent different random ensembles. Neither replaces the ordered competitors when looking for a zero-temperature ground state. The [SQS reference](../../../1_Tools/CAT/ordering/README.md#sqs-and-prescribed-short-range-correlations) contains the settings and diagnostics.

## 6. Establish numerical accuracy before production

Keep a convergence record in whichever format/location you use. An optional [record template](convergence-record.md) is provided. State the energy scale to resolve and the acceptable numerical uncertainty **before** choosing settings.

Use eV/atom or a clearly defined formula-unit normalization throughout. For a formula unit containing `n` atoms, `delta_E_meV_per_atom = delta_E_meV_per_formula_unit / n`. A target larger than the splitting of interest cannot resolve that splitting. Converge relative and mixing energies, as well as the forces/stress that determine relaxed geometries; a total-energy plateau for one pure structure is insufficient.

Choose pilots covering terminal compositions, different bonding/electronic characters, relevant parent phases, mixed orderings and representative large cells. The pilot set depends on your system; there is no mandatory number of phases or endmembers. Use one consistent physical model and potential choice per element.

In a pilot directory containing `POSCAR`, `POTCAR`, `INCAR`, `KPOINTS` and the reviewed `submit.sh` required by the scan generator:

```bash
cat-tools vasp converge --encut --range 1.0 2.0 0.25 --no-submit
```

This prepares an illustrative cutoff scan in multiples of the largest `POTCAR` ENMAX. It does not submit. Execute the prepared points with your local launcher or scheduler, wait for completion, and inspect them:

```bash
cat-tools vasp check --prefix ENCUT_ --verbose
```

Update the pilot's parent `INCAR` with the accepted explicit cutoff before testing sampling. For example:

```bash
cat-tools vasp converge --kpoints --man 40 50 60 70 --no-submit
```

These are trial reciprocal-density values, not accepted meshes. Execute, inspect, and update the baseline before proceeding. Test electronic tolerance, smearing and ionic stopping criteria in the same way; the [convergence workflow](../01_Convergence/VASP/README.md) gives the commands and full acceptance checks. Recheck cutoff at the final mesh and perform a combined stricter calculation on identical geometries.

Record the energy convention. The hull collector uses `Vasprun.final_energy` (normally the sigma→0 energy); the scan summary prints free energy TOTEN. In a smearing scan these are not interchangeable. Extend scans without a convincing plateau, and verify convergence again on accepted relaxed structures.

Keep a common explicit cutoff and a validated reciprocal-sampling policy across compared structures. Different supercells generally need different numerical meshes to give comparable sampling. Freeze accepted settings only after the pilots pass; retain potential provenance, actual meshes and stricter-reference results.

## 7. Prepare the selected calculations where you want to run them

Generation and analysis can run locally. DFT can run locally, on a cluster, or through your existing job manager. CAT does not require one of these arrangements.

For one selected ID, an example layout is:

```bash
mkdir -p calculations/ID/relax_01
cp pool-r01/candidates/ID/POSCAR calculations/ID/relax_01/POSCAR
```

Replace both occurrences of `ID` with the same selected ID from the manifest. Repeat for the candidates you will calculate, or use your existing staging script. Keep the mapping from candidate IDs to calculation directories. Do not replace the meaning of an existing ID when extending the search.

Supply your validated relaxation `INCAR`, matching `POTCAR`, cell-specific `KPOINTS`, and launcher in each calculation directory. Assemble `POTCAR` for the actual `POSCAR` species order; a single shared `POTCAR` does not fit different species lists. For example, **inside one prepared relaxation directory**:

```bash
cat-tools vasp potcar --sources
cat-tools vasp kpoints POSCAR DENSITY
```

Replace `DENSITY` with the validated density and use your configured potential library. These commands write `POTCAR`/`KPOINTS`; skip them if matching, verified inputs are already present. Review actual meshes and per-species/site settings before running.

For an **optional cluster transfer**, substitute your login and destination paths:

```bash
rsync -av --dry-run pool-r01/ USER@HOST:/path/to/destination/pool-r01/
rsync -av pool-r01/ USER@HOST:/path/to/destination/pool-r01/
rsync -av --dry-run calculations/ USER@HOST:/path/to/destination/calculations/
rsync -av calculations/ USER@HOST:/path/to/destination/calculations/
```

Preview each transfer before executing it. Transfer only what the receiving machine needs, retain manifests, and preserve any newer results already there. If some calculations finished locally, submit only the remaining candidates. Use a compatible installation on the receiving machine; a copied workstation virtual environment is not portable. Relative YAML paths make later generation easier to relocate; collection only needs saved manifests and calculation outputs.

## 8. Relax and verify each candidate

From a prepared calculation directory, use your tested VASP launcher. On Slurm:

```bash
sbatch submit.sh
```

For local execution, use your tested local command instead. Choose ionic/cell degrees of freedom for the physical question. Unconstrained zero-pressure bulk commonly uses `IBRION=2`, `ISIF=3`, sufficient `NSW` and negative `EDIFFG`; a fixed-cell or strained calculation defines a different, constrained comparison.

From the parent of example `relax_01`, `relax_02`, ... directories:

```bash
cat-tools vasp check --prefix relax_ --verbose
cat-tools vasp check --prefix relax_ --relax diff
```

Use your actual run prefix. Require electronic convergence, accepted forces and residual stress, sensible structure, and stability upon a restart from `CONTCAR`. Reaching `NSW` or scheduler completion does not establish convergence. Continue or verify in a new run directory, retaining previous outputs; follow the [relaxation procedure](../02_Relaxation/VASP/README.md).

Inspect final geometries and symmetry. Different parents may relax to the same minimum; different IDs do not necessarily remain distinct phases. Test appropriate symmetry-breaking starts where a high-symmetry stationary point is suspect. Disabling symmetry alone does not introduce a distortion, and relaxation alone does not replace a phonon stability check.

## 9. Calculate consistent static energies

Once a candidate's relaxation is accepted, create its static directory. The collector expects `RESULTS_ROOT/ID/vasprun.xml`; use any `RESULTS_ROOT` location. The following uses `statics` relative to your chosen working directory:

```bash
mkdir -p statics/ID
cp /path/to/accepted-relaxation/CONTCAR statics/ID/POSCAR
cp /path/to/accepted-relaxation/{INCAR,POTCAR,submit.sh} statics/ID/
cd statics/ID
```

Replace `ID` and the relaxation source path. Use a new directory, not an existing result. Edit the copied `INCAR` to retain the accepted physical/electronic settings and set an ordinary SCF static calculation:

```text
IBRION = -1
NSW = 0
ISIF = 2
ISTART = 0
ICHARG = 2
```

Remove relaxation-only `EDIFFG`/IWAVPR. These fresh-start settings do not need orbital or charge restart files. Regenerate a mesh for the accepted cell if appropriate:

```bash
cat-tools vasp kpoints POSCAR DENSITY
```

Replace `DENSITY` with your validated value. Run with your local launcher or scheduler, then verify SCF convergence and final forces/stress. Return to relaxation if the static settings produce material residual forces. Test stricter settings for close competitors and their reference endpoints on the same geometries. Repeat for each selected ID. The [static workflow](../03_StaticReference/VASP/README.md) provides the complete checks.

## 10. Collect energies and build the mixing hull

Run wherever CAT can access the manifest and completed statics. There is no need to move them if they are already accessible. If retrieving results from a cluster, preview the transfer first, for example:

```bash
rsync -av --dry-run USER@HOST:/path/to/results/statics/ statics/
rsync -av USER@HOST:/path/to/results/statics/ statics/
```

From the directory where `pool-r01` and `statics` are accessible:

```bash
cat-ordering collect pool-r01/manifest.json statics \
    --energy-model validated-protocol-r01 --output energies-r01.json
```

Replace the paths and protocol label with yours. Every selected candidate needs a converged ordinary static result. The collector checks composition and supported electronic settings/potential titles; it cannot validate the preceding relaxation, potential file identity or k-point convergence for you. DFT+U currently requires the audited-entry route described in the tools reference. Use fresh collection filenames for revisions.

List the collected IDs, compositions, normalized energies and final symmetries:

```bash
python - energies-r01.json <<'PY'
import json, sys
from pymatgen.core import Composition
for r in json.load(open(sys.argv[1]))['entries']:
    c = Composition(r['composition'])
    print(r['id'], c.reduced_formula, r['energy_eV'] / c.num_atoms,
          'eV/atom', r.get('final_spacegroup', 'not recorded'))
PY
```

Create `hull.yaml` in your preferred editor. Replace the endpoint placeholders with actual IDs at the terminal compositions of your chosen mixing space:

```yaml
mode: mixing
energies: energies-r01.json
terminal_ids: [ENDMEMBER_1_ID, ENDMEMBER_2_ID]
```

Use two terminals for a binary line; add the required terminals for a larger composition space. Select the lowest accepted energy at each terminal composition among the phases you included. Use `energies: [energies-r01.json, energies-r02.json]` to combine batches, with unique IDs and one consistent energy protocol. Energy-file paths are relative to `hull.yaml`. Do not include multiple revisions of the same ID.

```bash
cat-ordering hull hull.yaml --output hull-r01
cat hull-r01/hull.csv
```

`mixing_energy_eV_atom` references the lowest balanced mixture of the specified terminals. `energy_above_hull_eV_atom` references all supplied configurations. Decomposition weights in `hull.json` are atomic fractions. The solver uses actual cell atom counts; never pair a reduced formula with an unreduced cell total energy. Dependent terminal compositions are supported, but entries must lie within their mixing space. Inspect final structures when interpreting apparently equal phases.

## 11. Converge the sampled hull and extend the study

Treat these as separate checks:

| Check | Extend the calculations | Acceptance evidence |
| --- | --- | --- |
| Numerical accuracy | Tighter cutoff/mesh/SCF and relaxation checks for near-hull structures and terminals | Relative and mixing energies, forces and stress meet your recorded error budget. |
| Ordering coverage | More compositions, cell sizes/shapes, motifs and previously unselected candidates | The sampled hull stabilizes over the tested search domain; report remaining coverage limits. |
| Random ensemble, if used | Additional SQS sizes/shapes and independent seeds | Correlation residuals and observables stabilize to a stated accuracy. |

For a new search, copy and edit your configuration. Change the ID prefix, requested cells/compositions or method as needed, and choose a fresh output path:

```bash
cp ordering.yaml ordering-r02.yaml
# Edit ordering-r02.yaml before running the next command.
cat-ordering generate ordering-r02.yaml --output pool-r02
```

Calculate the new candidates, collect them into `energies-r02.json`, and include both accepted datasets in a new hull configuration. Compare hull vertices, decompositions and energy differences. Near-ties within measured numerical uncertainty remain unresolved; the `on_hull` solver tolerance is not a scientific error bar. Finite sampling cannot prove a ground state across all configurations.

This yields a zero-temperature mixing/ordering energy hull within the specified family. SQS do not automatically add configurational entropy. To assess stability against other compounds later, add consistently referenced elemental and competing phase calculations and use `mode: full`; see the [full-hull reference](../../../1_Tools/CAT/ordering/README.md#later-expansion-to-a-full-thermodynamic-hull).
