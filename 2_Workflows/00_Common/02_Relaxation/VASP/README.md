# VASP Relaxation

**Use:** Obtain a well-converged stationary structure under the chosen constraints, verify it survives a restart, and calculate a reproducible static reference. This stationary structure is then used as the input for all further characterization.

**Prerequisites:** Complete [setup](../../00_SettingUp/VASP/README.md) and [convergence](../../01_Convergence/VASP/README.md). Use the selected numerical settings, potentials and physical model from those checks. The templates assume non-magnetic bulk PBEsol with CPU/MPI VASP, please adapt it to your choice of system.

**Files:** [INCAR.relax](INCAR.relax) · [INCAR.static](INCAR.static) · [submit.sh](submit.sh) · [Legacy](Legacy/README.md).

**Navigation:** [Structure](#1-choose-and-inspect-the-starting-structure) · [Setup](#2-prepare-the-first-run) · [Settings](#3-choose-the-relaxed-degrees-of-freedom) · [Run](#4-submit-and-monitor) · [Acceptance](#5-check-forces-stress-and-structural-changes) · [Restart](#6-restart-and-verify-stability) · [Static reference](#7-calculate-the-final-static-reference).

## 1. Choose and Inspect the Starting Structure

Start from an experimental, a previously calculated, or a prepared model. Low-temperature experimental data can reduce the initial mismatch, but different functionals, pressures and magnetic states can favour different cells, so don't assume it to be the best choice for all cases. A database structure is often a good starting guess, but depending on the choice of setup, your minimum can be far off of theirs, always check for documentation of their setup! 
This setup guide covers CIF sources, VESTA, ASE and pymatgen conversion.

For a CIF conversion and an inspectable copy of a `POSCAR`:

```bash
cat-prepare structure convert starting.cif POSCAR
cat-prepare structure convert POSCAR starting.cif
```

These are alternative examples; choose distinct output names and keep your source files. Inspect species, ordering, occupancy treatment, interatomic distances and periodic images. Resolve partial occupancies before conversion to VASP input. Retain atom ordering through the workflow, especially with selective dynamics, MAGMOM, or downstream site-specific analysis.

A relaxation follows the starting basin and can preserve imposed symmetry. It does not guarantee the global minimum. If competing distortions or magnetic states matter, relax distinct candidates and compare them consistently. A high-symmetry stationary structure may still be a saddle point; phonons or other stability tests come later.

### Optional: Test a symmetry-breaking starting structure

If you are unsure about the true energy basin of your system, the [VASP structure-optimization guide](https://vasp.at/wiki/Structure_relaxation) recommends intentionally breaking symmetry in the starting structure rather than simply disabling symmetry with `ISYM=0`. CAT can prepare an independent trial:

```bash
cat-prepare structure perturb CONTCAR POSCAR.perturbed --amplitude 0.02 --seed 42
```

This moves each mobile atom by 0.02 Å in a random direction while preserving the cell, atom order and selective-dynamics constraints. The amplitude is an example, not a converged choice. Inspect interatomic distances, use the output as `POSCAR` in a separate relaxation directory, and compare several seeds against the unperturbed reference using consistent settings. Do not symmetry-standardize the perturbed structure before relaxing. Check the symmetry VASP detects in `OUTCAR`.

Atomic displacements only explore distortions compatible with the chosen cell and constraints. A one-atom primitive cell only undergoes a rigid translation; a larger supercell may be needed for symmetry-breaking modes, for these cases also check against the [ordering workflow](../../05_Ordering/README.md). To test cell-shape distortions, also consider `cat-tools vasp rescale` and suitable `ISIF` settings. Returning to the reference supports local stability but does not prove global stability or replace a phonon check.

## 2. Prepare the First Run

Keep successive runs in separate directories so inputs, geometry and results remain attributable:

```text
relaxation/
├── relax_01/       # initial relaxation
├── relax_02/       # verification or continuation from relax_01
└── static_01/      # final fixed-geometry reference
```

With CAT active and `PHONONTOOLS_ROOT` set:

```bash
workflow="$PHONONTOOLS_ROOT/2_Workflows/00_Common/02_Relaxation/VASP"
mkdir -p "$SCRATCH/calculations/my-system/relaxation"
cd "$SCRATCH/calculations/my-system/relaxation"
mkdir relax_01
cp "$workflow/INCAR.relax" relax_01/INCAR
cp "$workflow/submit.sh" relax_01/submit.sh
cp /path/to/checked/POSCAR relax_01/POSCAR
cp /path/to/converged/POTCAR relax_01/POTCAR
cp /path/to/converged/KPOINTS relax_01/KPOINTS
```

Replace the source paths. Copy the selected electronic/model settings into `relax_01/INCAR`, retaining the reviewed relaxation tags below. Do not replace a magnetic or +U calculation's INCAR wholesale with the nonmagnetic example. Be sure to set the tags according to your results from convergence tests.

If assembling potentials for a new setup, use `cat-tools vasp potcar --verbose` inside `relax_01`, then verify the release and variants against your convergence record. If generating the mesh, `cat-tools vasp kpoints POSCAR 50` is a starting point **only when density 50 has been validated for that cell**. It overwrites `KPOINTS` and creates a Gamma-centered mesh from reciprocal lengths without factor `2*pi`. For a changed cell or supercell, inspect the actual mesh and revalidate sampling rather than blindly copying a density or mesh from a different structure.

Review `submit.sh`: it retains an example PALMA-II VASP 6 module stack; verify availability. The four-rank, 8G, one-hour profile is for a small initial trial on an exclusive node, not a large-system resource prescription. Size longer runs using the [benchmark exercise](../../00_SettingUp/VASP/benchmarking.md). Confirm those modules, the account, partition, memory, wall time and rank layout for your allocation. Keep `NCORE`/`KPAR` consistent with the allocation; do not combine `NPAR` and `NCORE`.

## 3. Choose the Relaxed Degrees of Freedom

| Setting | Choice and consequence |
| --- | --- |
| `IBRION=2` | Conjugate-gradient structural optimization; a robust starting choice |
| `NSW=100` | Maximum ionic steps, not an accuracy criterion |
| `ISIF=2` | Relax positions at a fixed cell; calculate stress |
| `ISIF=3` | Relax positions, cell shape and volume for unconstrained bulk |
| `EDIFFG=-1E-03` | Example force-magnitude threshold of **0.001** eV/angstrom |
| `EDIFF` | Electronic noise must be sufficiently below the ionic/property error budget |

Use fixed-cell relaxation for an intentionally imposed cell, or an appropriate constrained procedure for slabs, interfaces and strain studies. Do not let an unconstrained bulk profile collapse vacuum or erase an intended strain. Selective dynamics controls allowed atomic motion; it does not by itself fix lattice degrees of freedom. Consult [ISIF](https://vasp.at/wiki/ISIF) for alternatives and version requirements.

Negative `EDIFFG` uses a force stopping criterion; positive `EDIFFG` uses energy changes between ionic steps. Choose the magnitude using the convergence workflow, and inspect constrained degrees of freedom separately. See [EDIFFG](https://vasp.at/wiki/EDIFFG).

The templates retain `PREC=Accurate`, `LASPH` and `ADDGRID`, with `LREAL=.FALSE.` as the reference. If production needs `LREAL=Auto`, use the tested per-species `ROPT` settings from convergence. Keep the chosen `ENCUT` explicit during cell optimization. A tight `EDIFFG` cannot compensate for inaccurate forces or basis-dependent stress!

Gaussian smearing is the example baseline. Use the occupation scheme and `SIGMA` validated for the system. Keep occupations consistent through verification. For metallic forces and stress, see the [smearing guidance](https://vasp.at/wiki/Smearing_technique).

## 4. Submit and Monitor

From your calculation directory you are ready to run—an initial relaxation, a continuation, or the final static calculation—submit with:

```bash
sbatch submit.sh
```

From the parent directory, use CAT to compare the latest results across your relaxation runs:

```bash
cat-tools vasp check --prefix relax_ --verbose
```

The table brings together the latest energy, maximum force, stress and available job information. Run it again as the calculations progress; after a run, add `--relax diff` to inspect how its structure changed. Choose a different prefix to summarize another group of directories. Manually prepared runs have no scan parameter, so those columns remain empty.

Use the trend in forces and energy to judge progress, and `OSZICAR`/`OUTCAR` to investigate steps that stall or oscillate. Section 5 describes the final acceptance checks.

### Watch the ionic steps

From the calculation directory or the parent of your numbered relaxation runs, export the recorded ionic steps and a gif with:

```bash
cat-tools vasp relaxation --out relaxation-view --gif --fps 6
```

Without a filename, CAT discovers and joins the segments automatically:

| Layout in the current directory | Order used |
| --- | --- |
| `OUTCAR` only | That single run |
| `OUTCAR.s1`, `OUTCAR.s2`, `OUTCAR.s10`, `OUTCAR` | Archives in numeric order, then the current `OUTCAR` |
| `OUTCAR.sN` files without a current `OUTCAR` | Archives in numeric order |
| No local `OUTCAR` or archives; `relax_1/`, `relax_2/`, `relax_10/` | Numbered directories in numeric order; each run's archives precede its `OUTCAR` |

Names such as `relax_01` also work. Other directory names are ignored. A local `OUTCAR`/archive chain takes precedence over child directories. Omit `--gif` to export just the trajectory files and step summary.

You can also supply a directory or multiple files explicitly:

```bash
cat-tools vasp relaxation /path/to/relaxation --gif
cat-tools vasp relaxation relax_01/OUTCAR --gif    # Only this file
cat-tools vasp relaxation OUTCAR.s* OUTCAR --gif
cat-tools vasp relaxation relax_*/OUTCAR --gif
```

Explicit files use the same numeric suffix/directory ordering; selecting just `OUTCAR` exports only that file. Each segment uses its own adjacent `POSCAR` by default; `--poscar PATH` supplies a shared species/order reference. Atom count and species sequence !must agree! Keep atom identity/order unchanged across restarts; permutations of identical elements cannot be detected automatically.

The combined trajectory has continuous step numbers plus source file, segment and local-step metadata. PhonoVis labels the source segment; `steps.csv` retains the same provenance. Repeated boundary geometries are kept, since a restart can recompute forces on the same geometry. Stride applies to the combined sequence. Segment numbering follows filenames, not timestamps, so ensure the names represent the intended chronology.

Install CAT's `structure` extra for export and `visualization` extra for viewing/GIF rendering. For example, from the repository, `uv sync --extra structure --extra visualization`.

The output directory contains:

- `relaxation.extxyz`: coordinates, full cell vectors, forces and ionic-step indices for PhonoVis/ASE.
- `relaxation.xyz` and `relaxation.vmd`: a VMD trajectory and loader, including per-frame cell lengths/angles.
- `steps.csv`: global step, maximum force, volume, segment number, local step and source file.

Open the trajectory in PhonoVis, or load it in VMD:

```bash
cat-phonovis trajectory relaxation-view/relaxation.extxyz
vmd -e relaxation-view/relaxation.vmd
```

For a GIF directly during export:

```bash
cat-tools vasp relaxation --out relaxation-view --gif --fps 6
```

Or render the exported trajectory later, including on a machine without VASP:

```bash
MPLBACKEND=Agg cat-phonovis trajectory relaxation-view/relaxation.extxyz \
    --gif relaxation-view/relaxation.gif --fps 6
```

PhonoVis shows atoms and the changing cell with fixed plot bounds, plus step number, maximum force and volume. The animation follows recorded geometries, without displacement amplification or interpolation. These are optimizer steps, **not physical time**. `--stride 5` exports every fifth block and always includes the last complete one. Coordinates are kept as written; periodic boundary crossings can therefore look like jumps, no handling of pbc (yet). VMD uses spheres rather than inferred bonds, avoiding misleading fixed bonds during rearrangements; its XYZ does not retain the full lattice orientation, which is preserved only in the extxyz.

An unfinished final block in a growing OUTCAR is skipped with a warning, allowing a snapshot of progress. Re-run the export to update it; outputs in the chosen directory are replaced. The exporter expects standard VASP POSITION/TOTAL-FORCE blocks and a matching species order, and does not infer SCF convergence.

## 5. Check Forces, Stress, and Structural Changes

Before accepting a run, establish that the electronic cycles converged, and the ionic optimizer reached the intended accuracy rather than its step limit. Look for VASP's structural-accuracy message, but also inspect the actual final quantities and constraints.

| Check | Where and how to interpret it |
| --- | --- |
| Forces | Final `POSITION ... TOTAL-FORCE (eV/Angst)` block; inspect vector norms and allowed components |
| Stress | Final `in kB` line, in XX, YY, ZZ, XY, YZ, ZX order; `kB` means kilobar |
| Cell and positions | `CONTCAR` versus starting `POSCAR`; inspect unexpected collapse, rearrangement or symmetry retention |
| Electronic state | Magnetization/occupations and SCF history; confirm the intended state survived |
| Termination | Scheduler state, `OUTCAR` timing footer, SCF and ionic convergence evidence |

For zero-pressure full-cell relaxation, check diagonal and shear stresses against your tolerance. For imposed pressure or constraints, judge the stress relative to that target and only the cell degrees of freedom you actually relaxed. A fixed-cell run need not have zero stress.

CAT's MaxF and MaxStress are absolute maxima, not errors against a reference or a pressure-adjusted residual. With selective dynamics, an all-component maximum may include constrained forces. Check cutoff dependence and restart behavior; see [Pulay stress](https://vasp.at/wiki/Pulay_stress).

For a periodic structure comparison preserving site identity:

```bash
cat-tools vasp compare relax_01/POSCAR relax_01/CONTCAR \
    --same-order --top 10 --csv relax_01/displacements.csv
cat-tools vasp check --prefix relax_ --relax diff
cat-tools vasp check --prefix relax_ --relax abs
```

`compare` handles periodic images, unlike a text diff. It assumes the same cell setting and does not search cell rotations or supercell mappings. Inspect lattice changes as well as atomic displacements: with a changing cell, displacement lengths depend on the reference lattice. Use `--align-origin` only when removing a global translation is justified. Compare constrained/site-specific structures with atom correspondence intact.

## 6. Restart and Verify Stability

First ensure the old job is no longer running and CONTCAR is complete and parseable. Preserve the old directory and create a new one:

```bash
mkdir relax_02
cp relax_01/CONTCAR relax_02/POSCAR
cp relax_01/{INCAR,POTCAR,KPOINTS,submit.sh} relax_02/
```

Review inputs, then submit `relax_02` using section 4. This intentionally starts the electronic problem afresh with the template's `ISTART`/`ICHARG`, while using the previous final geometry. It also rebuilds the basis for the current cell. If reusing `WAVECAR`/`CHGCAR` for speed, deliberately set compatible restart tags; those files are not copied by this example.

Keep the mesh initially fixed to isolate restart behaviour. If the cell changed enough to warrant a new mesh, test that separately. After completion:

```bash
cat-tools vasp compare relax_01/CONTCAR relax_02/CONTCAR --same-order --top 10
cat-tools vasp check --prefix relax_ --relax ref relax_02
cat-tools vasp check --prefix relax_ --verbose
```

The reference directory is explicit because manually prepared runs have no scan metadata for automatic selection. Accept the result when force/stress criteria hold and further movement/cell changes are within the tolerance required for your property. If not, diagnose numerical noise, cell-basis effects, insufficient steps, constraints or a changing electronic state, and continue in another numbered directory. Small movement alone does not prove convergence.

### CAT's in-place restart alternative

For an existing in-place workflow, the maintained replacement for `varc` is:

```bash
cat-tools vasp restart --finished-only relax_01
```

This **modifies the directory**: it archives selected outputs as `*.sN` and copies `CONTCAR` to `POSCAR`. Use it only after the job stops. `--finished-only` checks normal termination, not SCF or ionic convergence. An interrupted run with a verified usable `CONTCAR` can be handled without that option; do not force past a bad geometry.

The helper archives `OUTCAR`, `CONTCAR`, `vasprun.xml`, `out_VASP` and out when present. It saves the original `POSCAR` once as `POSCAR_init` if that file does not already exist, but does not retain every successive input `POSCAR`. It does **not** archive `INCAR`, `POTCAR`, `KPOINTS`, `submit.sh`, `OSZICAR`, this template's `vasp.out`/`vasp.err`, slurm logs or restart files. Preserve those separately before changing inputs or rerunning if desired; separate run directories above are the complete-history approach. `--backup-only` archives without promoting `CONTCAR`, and rerunning against identical archived outputs is a no-op, including after backup-only. The helper neither submits a job nor checks the scheduler.

## 7. Calculate the Final Static Reference

The dedicated [static-reference workflow](../../03_StaticReference/VASP/README.md)
also covers output selection and the handoff to property calculations. The steps below
provide the direct continuation from this relaxation.

Use the accepted final geometry in a new directory. Assuming `relax_02` passed verification:

```bash
mkdir static_01
cp relax_02/CONTCAR static_01/POSCAR
cp relax_02/{INCAR,POTCAR,KPOINTS,submit.sh} static_01/
```

Edit `static_01/INCAR`, retaining all selected electronic and physical-model settings, to set:

```text
SYSTEM = static-reference
NSW = 0
IBRION = -1
ISIF = 2
```

Remove relaxation-only `EDIFFG` and any explicit `IWAVPR` carried from an old input. Keeping `ISIF=2` allows the final stress check. [INCAR.static](INCAR.static) shows the complete baseline, but copying it wholesale could discard your selected model/settings. Set `LWAVE`/`LCHARG` according to downstream needs: the example disables both, whereas bonding or charge analysis may require them.

Submit `static_01` using section 4 and inspect:

```bash
cat-tools vasp check --prefix static_ --verbose
```

Verify the final static forces and stresses still satisfy the intended accuracy. If denser sampling or different occupations change the forces materially, another relaxation may be needed.

`OUTCAR` reports `free energy TOTEN`, energy without entropy, and `energy(sigma->0)`. CAT reports **TOTEN**. Use the same energy definition and smearing treatment when comparing calculations; the values need not coincide. Relative phase energies and convex hulls require consistent reference compositions, functional, potentials, numerical accuracy and energy treatment across all phases.

A converged relaxation establishes a stationary point under its constraints. It does not establish dynamical stability or include vibrational/thermal free energies. Continue with the relevant phonon, electronic or other property workflow and its convergence checks.
