# VASP Convergence

**Use:** Establish numerical accuracy at a reasonable cost before production. Converge the static electronic problem, relax, validate on the resulting structure, then test the target property. For phonons, energy convergence alone is insufficient: small errors in forces can change the inferred force constants substantially so be aware.

**Prerequisites:** Complete [VASP setup](../../00_SettingUp/VASP/README.md), activate CAT, and prepare a checked structure and matching potentials. The examples assume bulk, nonmagnetic PBEsol and CPU/MPI VASP with slurm. Adapt the physical model before scanning; convergence cannot fix the wrong magnetic state, functional, potential, or structure.

**Files:** [INCAR.static](INCAR.static) · [submit.sh](submit.sh) · [optional convergence record](convergence-record.md) · [original scripts and guide](Legacy/README.md).

**Navigation:** [Criteria](#1-decide-what-must-converge) · [Setup](#2-prepare-a-working-directory) · [Run and inspect](#3-generate-submit-and-read-a-scan) · [ENCUT](#4-converge-the-plane-wave-cutoff) · [k-points](#5-converge-the-k-point-mesh) · [Electronic settings](#6-test-the-remaining-electronic-settings) · [Relaxation](#7-converge-the-relaxed-structure) · [Validation](#8-validate-and-freeze-the-production-inputs) · [Library](#9-keep-an-optional-convergence-library) · [Property checks](#10-converge-the-target-property).

## 1. Decide What Must Converge

Choose tolerances for the quantity you will publish, not just for the total energy. A phase-energy comparison, a soft phonon mode, and an MD trajectory can require different settings. Write down the tolerances before choosing a plateau. If using the optional [convergence library](#9-keep-an-optional-convergence-library), create your record now and fill it as the scans finish. For an initial screen, 1 meV/atom in energy and 0.001 eV/angstrom in force differences can be useful **illustrative targets**, but they are not guarantees for any particular property.

Compare each candidate with a stricter reference on the **same geometry, atom order, and physical state**:

| Quantity | Comparison | Why it matters |
| --- | --- | --- |
| Energy | `abs(E - E_ref) / N_atoms` | Numerical noise relative to the energy differences of interest |
| Forces | `max_i norm(F_i - F_i_ref)` | Accuracy of relaxation and force constants |
| Stress | `max_ab abs(stress_ab - stress_ref_ab)` | Cell relaxation, strain and elastic response |

Use more than one stricter point to establish a plateau. Highly symmetric structures can have zero forces by symmetry even at poor settings: include a representative small displacement, and strain where relevant, when testing force/stress accuracy. Keep that test geometry fixed throughout each static comparison.

The initial one-parameter scans are an efficient screen, not proof that parameters are independent. Recheck the selected cutoff at the final *k* mesh, and test the *k* mesh together with smearing for metals. Preserve the same spin/charge state across comparisons; an electronic-state change needs investigation rather than averaging into a convergence curve.

## 2. Prepare a Working Directory

Work outside the repository's templates. With `PHONONTOOLS_ROOT` set as in the setup guide:

```bash
workflow="$PHONONTOOLS_ROOT/2_Workflows/00_Common/01_Convergence/VASP"
mkdir -p "$SCRATCH/calculations/my-system/convergence"
cd "$SCRATCH/calculations/my-system/convergence"
cp "$workflow/INCAR.static" INCAR
cp "$workflow/submit.sh" submit.sh
cp /path/to/checked/POSCAR POSCAR
cat-tools vasp potcar --sources
cat-tools vasp kpoints POSCAR 50
```

Replace the structure path. If reproducing an existing calculation, copy its matching POTCAR instead of rebuilding it. The builder reports each potential's ENMAX and the maximum; `--sources` additionally records source paths. Keep the potential release, variants and species order fixed throughout the tests.

Review `INCAR` and **edit `submit.sh` before generating any scan**: add the required account/partition, and set resources for the actual node. Its four active ranks, 8G memory and 30-minute limit are a small-cell pilot budget; exclusive allocation still reserves the node. Use the [benchmark exercise](../../00_SettingUp/VASP/benchmarking.md) to size representative points, including the largest cutoff/mesh, before the scan. Keep the allocation and VASP parallelization fixed while comparing numerical settings. Benchmark NCORE/KPAR separately if optimizing performance; changing them can also change the actual number of bands. See the setup guide's [parallelization discussion](../../00_SettingUp/VASP/README.md#7-match-parallelization-to-the-node-layout).
Add magnetism, SOC, +U, dispersion or other model settings where needed.

CAT expects these filenames in the working directory:

| File | Role |
| --- | --- |
| `POSCAR` | Identical starting geometry for every point in a scan |
| `POTCAR` | Identical potential datasets |
| `INCAR` | Current baseline, updated after each accepted scan |
| `KPOINTS` | Current baseline mesh; not required for the k-point scan itself |
| `submit.sh` | Required job script, including for `--no-submit` |

The generator copies only these inputs. It does not copy WAVECAR, CHGCAR, POTCAR.sources, or additional method-specific files. The template's `ISTART=0`, `ICHARG=2` therefore makes each point an independent start. If your method needs extra files, place them in every generated directory before submission. Preserve provenance in the parent directory.

CAT's *k*-point helper writes a Gamma-centered mesh, rounding density times each reciprocal-vector length **without `2*pi`**, with a minimum of one point per direction. For orthogonal cells, density 50 roughly means `N_i * L_i = 50 angstrom`; it is not 50 *k*-points per atom or per unit cell. It expects a single positive POSCAR scale factor. Normalize other POSCAR scale conventions first. Inspect skewed cells and manually restrict vacuum directions for slabs/wires; the isotropic density scan is primarily a bulk workflow.

## 3. Generate, Submit, and Read a Scan

`cat-tools vasp converge` scans one axis per invocation. **It submits jobs by default.** All generation examples below use `--no-submit` so you can inspect the inputs and resource requests first.

| Axis | `--range LOW HIGH STEP` | `--man VALUES...` |
| --- | --- | --- |
| `--encut` | Multiples of largest POTCAR ENMAX | Absolute cutoff in eV |
| `--kpoints` | Reciprocal-density values | Reciprocal-density values |
| `--ediff` | Electronic thresholds in eV | Usually logarithmically spaced thresholds |
| `--ediffg` | Ionic thresholds | Negative values for force criteria in eV/angstrom |
| `--sigma` | Smearing widths in eV | Smearing widths in eV |
| `--nbands` | Integer band counts, only useful for bonding atm | Integer band counts |

`--prefix NAME_` gives a new output prefix for repeats. Existing directories cause an error rather than being overwritten. `--set TAG=VALUE` is repeatable and applies to generated inputs only; the scanned parameter takes precedence. Calculation profiles are:

| Option | Tags set by CAT |
| --- | --- |
| `--calc static` (default) | `IBRION=-1`, `NSW=0`, `ISIF=2` |
| `--calc ions` | `IBRION=2`, `NSW=100`, `ISIF=2` |
| `--calc full` | `IBRION=2`, `NSW=100`, `ISIF=3` |
| `--calc keep` | Keeps supplied relaxation tags |

Explicit `--set` values override the profile. For example, `--calc full --set NSW=200` allows more ionic steps. CAT rejects an EDIFFG scan with the static profile; with `keep`, you must supply actual relaxation settings yourself.

### Submit inspected directories

After generating the ENCUT scan in section 4, inspect representative `INCAR`, `KPOINTS`, `.convgen_meta`, and `submit.sh` files. Then run this Bash block from the parent directory, changing `prefix` for each subsequent scan:

```bash
(
    set -euo pipefail
    shopt -s nullglob
    prefix=ENCUT_
    dirs=("${prefix}"*/)
    ((${#dirs[@]})) || { echo "No matching directories" >&2; exit 1; }
    for dir in "${dirs[@]}"; do
        [[ -f "$dir/.convgen_meta" && -f "$dir/submit.sh" ]] || exit 1
        [[ ! -e "$dir/.convgen_jobid" ]] || {
            echo "Already submitted: $dir; inspect before resubmitting" >&2
            exit 1
        }
    done
    for dir in "${dirs[@]}"; do
        (
            cd "$dir"
            response=$(sbatch --parsable submit.sh)
            jobid=${response%%;*}
            printf '%s\n' "$jobid" > .convgen_jobid
            printf 'JOBID=%s\n' "$jobid" >> .convgen_meta
            printf '%s: %s\n' "$dir" "$jobid"
        )
    done
)
```

The stored job IDs let CAT retrieve Slurm accounting. Run the block once per scan; if submission fails partway through, check Slurm and submit only the remaining directories. Do not remove job markers or rerun blindly. Alternatively, omit `--no-submit` on the **first** generation call once you trust the inputs. Removing it on a repeated call does not submit existing directories.

### Judge results, not just job completion!

```bash
cat-tools vasp check --prefix ENCUT_ --verbose
```

CAT reports scan/effective values, last OUTCAR `free energy TOTEN`, energy per atom, maximum force magnitude, maximum absolute stress component, and available resource data. `--verbose` adds Slurm state and allocation details. Wall/CPU timing can fall back to OUTCAR when accounting is unavailable. MaxRSS/task is a task-level peak, not total node memory; `kB` in the stress column means **kilobar**, with 1 kbar = 0.1 GPa.

**MaxF and MaxStress are absolute values, not differences from a reference.** A stationary maximum can hide changes in other atoms/components. For final acceptance, compare the full force vectors and stress tensors using the definitions in section 1, then record those differences separately. Non-zero forces on the fixed unrelaxed structure are expected.

## 4. Converge the Plane-Wave Cutoff

Start with the provisional mesh from section 2 and otherwise fixed settings:

```bash
cat-tools vasp converge --encut --range 1.0 2.5 0.25 --no-submit
# Inspect and submit with prefix=ENCUT_ using section 3; wait for completion.
cat-tools vasp check --prefix ENCUT_ --verbose
```

This preserves the original broad ENMAX-multiplier screen. It can be expensive: narrow it using previous validated experience, or extend it when necessary. For absolute values, use `--encut --man 400 450 500 550 600` instead, with values appropriate to your potentials. Range and manual modes deliberately use different units.

Select the lowest cutoff meeting your energy, force and stress tolerances against stricter points. Keep ENCUT explicit and identical for structures whose energies you compare. Cell optimization is particularly sensitive to basis incompleteness; infer Pulay-stress effects from cutoff dependence of stress/geometry rather than treating a printed pressure correction as a measured error. See [ENCUT](https://vasp.at/wiki/ENCUT) and [Pulay stress](https://vasp.at/wiki/Pulay_stress).

**Before continuing:** edit the parent `INCAR` to the selected ENCUT. CAT does not promote a result or modify the parent automatically. Leave the original scan directories intact as evidence.

## 5. Converge the k-Point Mesh

With the selected cutoff in the parent INCAR:

```bash
cat-tools vasp converge --kpoints --range 30 70 10 --no-submit
# Inspect and submit with prefix=KPOINTS_; wait for completion.
cat-tools vasp check --prefix KPOINTS_ --verbose
```

Inspect the actual meshes in the generated KPOINTS files: several density values can round to the same mesh and are then duplicate accuracy tests. Mesh parity, symmetry, and Fermi-surface sampling can produce non-monotonic behavior; look for stability across several distinct meshes. Record centering, shift, actual mesh and irreducible *k*-point count, not just the CAT density. For a custom Monkhorst-Pack or anisotropic sequence, prepare separate directories with explicitly chosen KPOINTS files instead.

**Before continuing:** regenerate the parent KPOINTS using your selected density. For example, **only if 60 was accepted**:

```bash
cat-tools vasp kpoints POSCAR 60
```

This overwrites the parent KPOINTS. Recheck the chosen ENCUT against a stricter cutoff at this mesh, using a fresh prefix such as `ENCUT_finalmesh_`. For metals, repeat relevant meshes at the selected smearing width in section 6.

## 6. Test the Remaining Electronic Settings

Keep the selected ENCUT and KPOINTS in the parent. Generate, inspect, submit, and check each relevant scan in turn; update the parent INCAR after accepting a setting. The ranges below are examples, not universal requirements.

### Electronic stopping criterion: EDIFF

```bash
cat-tools vasp converge --ediff --man 1e-4 1e-5 1e-6 1e-7 1e-8 --no-submit
# Submit with prefix=EDIFF_; wait for completion.
cat-tools vasp check --prefix EDIFF_ --verbose
```

EDIFF is an SCF stopping threshold in eV for the whole cell, not an energy-per-atom error bound. Choose it so SCF noise is below the force/property error budget. Confirm tight points really reach their threshold; increasing NELM alone does not repair an unstable electronic minimization. Force-difference calculations often need tighter electronic convergence than routine total-energy work. See [EDIFF](https://vasp.at/wiki/EDIFF).

### Occupations and smearing: SIGMA

```bash
cat-tools vasp converge --sigma --man 0.20 0.10 0.05 0.02 --no-submit
# Submit with prefix=SIGMA_; wait for completion.
cat-tools vasp check --prefix SIGMA_ --verbose
```

Keep ISMEAR fixed within this scan. For Gaussian smearing, inspect the extrapolated `energy(sigma->0)` as well as forces and stress; CAT's table uses TOTEN, which changes with smearing. Reducing SIGMA requires adequate *k*-point sampling, particularly in metals. Choose the occupation scheme for the system: avoid positive ISMEAR for gapped systems; tetrahedron methods do not provide a SIGMA scan, and metallic relaxation forces need special care. If Fermi-Dirac smearing represents a physical electronic temperature, lowering it changes that model, keep this in mind! See the [smearing guide](https://vasp.at/wiki/Smearing_technique).

### Number of bands: NBANDS, when relevant

Ground-state calculations often need only verification that enough empty bands are available and SCF convergence is sound. Projection, spectroscopy and excited-state methods can impose additional requirements; for LOBSTER, choose the band count consistently with the selected local basis and inspect projection quality.

Read the actual NBANDS in a successful OUTCAR first. If 100–400 bands is a physically appropriate range for **your** system, for example:

```bash
cat-tools vasp converge --nbands --range 100 400 50 --no-submit
# Submit with prefix=NBANDS_; wait for completion.
cat-tools vasp check --prefix NBANDS_ --verbose
```

Replace those counts based on electron count, spin and target observable; do not start below the occupied-state requirement, this will break your calculation. Record the **actual** NBANDS from OUTCAR because VASP can increase the requested count for its parallel layout. Converge the unoccupied-state observable itself where applicable, rather than total energy alone. See [NBANDS](https://vasp.at/wiki/NBANDS).

### Real-space projectors and other numerical choices

Keep `LREAL=.FALSE.` as the reference here. If production needs `LREAL=Auto`, benchmark it at the selected ENCUT/PREC and tune **ROPT explicitly** against that reference. Record one ROPT value per POTCAR species, in species order; its defaults depend on PREC. This is a separate test, not one of CAT's six scan axes. See the setup guide's LREAL/ROPT discussion and the [ROPT reference](https://vasp.at/wiki/ROPT). Likewise, assess changes to ADDGRID or precision settings before carrying them into production.

## 7. Converge the Relaxed Structure

Once the static problem is controlled, compare relaxations starting from the **same parent POSCAR**:

```bash
cat-tools vasp converge --ediffg --man -0.02 -0.01 -0.005 -0.002 \
    --calc full --no-submit
# Submit with prefix=EDIFFG_; wait for completion.
cat-tools vasp check --prefix EDIFFG_ --relax diff
cat-tools vasp check --prefix EDIFFG_ --relax abs
cat-tools vasp check --prefix EDIFFG_ --relax ref
```

Use `--calc ions` for fixed-cell internal relaxation. The `full` profile assumes all cell degrees of freedom should relax; for slabs, constraints, or other cell choices, supply suitable tags with `--set` or use `--calc keep` with a reviewed INCAR. Negative EDIFFG specifies a force threshold; positive values use an energy-change criterion and answer a different question. See [EDIFFG](https://vasp.at/wiki/EDIFFG).

`diff` compares each initial POSCAR to its own CONTCAR; `abs` reports final cell parameters; `ref` compares final structures with the automatically selected strictest scan point. You can specify a reference explicitly with `--relax ref DIRECTORY`. Verify that the reference itself converged, e.g. through very strict settings.

Check whether cell parameters, positions and the target property stabilize. EDIFFG is not an independent stress tolerance: inspect residual stress relative to the intended external pressure, as well as forces. A relaxation that stops at NSW needs continuation or diagnosis. Different minima across runs indicate a structural/state-selection problem in addition to numerical convergence which needs separate treatment.

For variable-cell runs, restart from the final CONTCAR with the selected settings and verify stability; rebuilding the basis for the new cell helps expose finite-basis effects. Do this in a new directory, preserving the original relaxation outputs.

## 8. Validate and Freeze the Production Inputs

On the accepted relaxed geometry, compare static calculations using selected and stricter settings. Both must use the **same final POSCAR**; comparing energies on two differently relaxed structures mixes numerical and structural effects.

Create a separate validation directory with that CONTCAR copied as POSCAR, the same POTCAR, your selected parent INCAR/KPOINTS, and submit.sh. For example, **if 540 eV was selected and 650 eV is an appropriate stricter check**:

```bash
cat-tools vasp converge --encut --man 540 650 --prefix ValidateENCUT_ --no-submit
```

Repeat for a denser mesh/tighter electronic threshold as needed, and include a combined stricter calculation to catch interacting errors. Regenerate a density-based mesh for a changed cell deliberately and record the actual result. Validate a displaced geometry too if the final structure has symmetry-enforced zero forces.

If validation fails, revisit the relevant scans and relax again when the new settings change forces or stress materially. Otherwise perform/confirm the final relaxation at the chosen production settings and save the final input set and convergence evidence. Related compounds can provide starting guesses, but changes in volume, potential variants, magnetic state, metallicity or supercell size require new checks. Production settings for compared phases must give consistent relative accuracy.

## 9. Keep an Optional Convergence Library

A folder of Markdown records provides a human-readable, editable library on a workstation or cluster. No database or extra package is required. Keep one record per phase/model/property and hardware context; start a new record when those change substantially.

For a personal library:

```bash
export CAT_CONVERGENCE_LIBRARY="$HOME/convergence-library"
mkdir -p "$CAT_CONVERGENCE_LIBRARY"
cp -n "$PHONONTOOLS_ROOT/2_Workflows/00_Common/01_Convergence/VASP/convergence-record.md" \
    "$CAT_CONVERGENCE_LIBRARY/my-system-PBEsol-phonons-PALMA.md"
```

Edit that copy with your usual text editor (vim my beloved <3). Optionally keep the export in `~/.bashrc`. `CAT_CONVERGENCE_LIBRARY` is only a convenience variable used in these instructions; CAT does not read it or silently reuse recorded settings. For a group library, point it to an agreed writable project directory and choose a unique descriptive filename. Existing records should remain attributable; please record corrections and superseding results!

The [template](convergence-record.md) records acceptance criteria, selected and reference settings, geometry/potential provenance, validation status, and timing context. Populate it as you proceed, including unsuccessful points. Save the raw tables beside the record, using a unique name for each scan:

```bash
cat-tools vasp check --prefix ENCUT_ --verbose \
    > "$CAT_CONVERGENCE_LIBRARY/my-system-PBEsol-phonons-PALMA-ENCUT.txt"
sha256sum POSCAR POTCAR INCAR KPOINTS submit.sh
```

Copy the hashes into the record or retain them with the archived inputs. Save final-input hashes again after validation. Keep portable relative paths where possible; a cluster-only absolute path is insufficient for someone receiving just the record.

Timing rows include wall time, measured CPU time, allocated core-hours, memory, iteration counts and repeated-run variation. Separate static timings from complete relaxations. Compare like hardware/builds and parallel layouts; a faster run caused by fewer iterations or a different restart policy is not evidence of better scaling.
The record helps another user choose an initial cutoff/mesh and resource request, but they should run a short selected-versus-stricter test on their own system before adoption. Share metadata and permitted inputs without bundling licensed potentials.

For performance tuning, strongly prefer the group's shared
[benchmark log](../../00_SettingUp/VASP/benchmark-log.md), with access for all
colleagues. Link the shared table from this convergence record rather
than keeping the evidence only in a personal library. Follow the
[logging procedure](../../00_SettingUp/VASP/benchmarking.md#5-share-results-and-improve-settings-together)
to add simple rows for settings, resources and timings, including failed runs.

## 10. Converge the Target Property

Stable electronic inputs are the starting point for the next workflow, not its endpoint. Depending on the calculation, test:

- Phonon supercell size and displacement amplitude; electronic force noise versus the displacement signal.
- Thermal transport q-mesh, third-order interaction cutoff and supercell size.
- NEB image count and force thresholds along the path.
- MD cell size, time step, equilibration and sampling duration.
- Spectroscopy empty-state manifold, sampling and method-specific cutoffs.

Measure convergence of the actual observable and append the result or link its separate record. This is what turns a reusable numerical starting point into a validated production workflow.

## References and Further Reading

The linked VASP pages above give tag definitions and method limitations. The original workflow's broader reading is retained here:

1. [Janssen et al., automated optimization and uncertainty quantification of plane-wave convergence parameters, npj Computational Materials 10, 263 (2024)](https://doi.org/10.1038/s41524-024-01388-2).
2. [Choudhary and Tavazza, k-point and cutoff convergence in high-throughput DFT, Computational Materials Science 161, 300–308 (2019)](https://doi.org/10.1016/j.commatsci.2019.02.006).
3. [Prandini et al., precision and efficiency of solid-state pseudopotentials, npj Computational Materials 4, 72 (2018)](https://doi.org/10.1038/s41524-018-0127-2).
4. [Bosoni et al., reproducible verification of DFT implementations, Nature Reviews Physics 6, 45–58 (2024)](https://doi.org/10.1038/s42254-023-00655-3).
5. [Lejaeghere et al., reproducibility of DFT calculations of solids, Science 351, aad3000 (2016)](https://doi.org/10.1126/science.aad3000).
