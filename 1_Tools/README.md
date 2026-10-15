## CAT Preparation and Calculation Tools

**Use:** Prepare structures and calculation inputs, generate and audit VASP/Phono3Py jobs, prepare restarts, and inspect slurm logs. These tools can run locally or on the HPC, some of them might be helpful during analysis as well.

**Prerequisites:**

- Python 3.12 or newer for the installed CAT toolbox. Install it from the repository root, not from this folder.
- Bash 4 or newer and standard Linux utilities for the shell-based VASP and Phono3Py commands.
- Your own VASP potential library for `POTCAR` generation; it is not included in the repository.
- Slurm `sbatch` for submission and `sacct` for accounting data. Log scanning works without a live scheduler.
- The correct module environment on the nodes that execute generated Phono3Py jobs, use `caterpillar search <module_name>` when unsure.
- Install Python dependencies from the repository-level `pyproject.toml` and `uv.lock`; optional features use the documented extras.
- The `structure` extra for Pymatgen/ASE structure conversion; Windows SoftBV/BVPA and WSL for the BVPA launcher.

**Folder Setup:**

| Folder | Contents | Role |
|:-------|:---------|:-----|
| [CAT](CAT/README.md) | `cat_utilities`, `cat_prepare`, and `cat_shared` packages | Maintained implementations, reusable functions, and packaged shell scripts. |
| [General](General/README.md) | Standalone `SlurmScan.py` | Independent tool used by the `0_bashrc.d` scanner aliases. |
| [legacy](legacy/README.md) | Archived VASP, Phono3Py, and bonding implementations | Excluded from installation; use the maintained CAT commands except for legacy workflows. |

---

## Setup

From your repository checkout:

```bash
cd /path/to/PhononTools
uv sync --locked --extra structure
source .venv/bin/activate
cat-tools --help
```

Alternatively, install into an existing Python environment:

```bash
python -m pip install -e '.[structure]'
```

The single repository-level installation also provides the analysis commands in [3_Analysis](../3_Analysis/README.md). Add the `ordering`, `bonding`, `visualization`, or `mpi` extras as needed for those tasks; MPI also needs a working system MPI runtime e.g. `MPI4Py`.

`cat-tools` groups related commands. `cat-prepare` remains available as a shorter equivalent of `cat-tools prepare`. Help is availabe from:

```bash
cat-tools vasp --help
cat-tools phono3py --help
cat-tools prepare --help
cat-prepare --help
```

Run commands from the calculation directory unless you supply the command's input paths or `--root` option. `cat-tools` keeps your current working directory and forwards arguments and exit status to the selected utility. Shell scripts are run with Bash, not sourced into your interactive shell.

SlurmScan is independent of CAT and uses only the Python standard library:

```bash
python3 /path/to/PhononTools/1_Tools/General/SlurmScan.py --root /path/to/runs
```

This enables it to work independently of CAT.

The former VASP and Phono3Py implementations are archived in [legacy](legacy/README.md). Use of CAT commands is preferred for maintained versions.

---

## Command Overview

| Command | Purpose | Writes or runs |
|:--------|:--------|:---------------|
| `cat-tools vasp potcar` | Select and concatenate VASP potentials. | Replaces `POTCAR` and its source record after validation, replaces `apgu`. |
| `cat-tools vasp kpoints` | Generate a reciprocal-density mesh. | Overwrites `KPOINTS`, replaces `GenKPOINTS.sh`. |
| `cat-tools vasp converge` | Prepare a parameter convergence scan. | Creates directories and **submits by default**; use `--no-submit` for generation only, replaces `CovGen.sh`. |
| `cat-tools vasp check` | Summarize convergence results. | Prints a report; queries accounting when available, replaces `ConvCheck.sh`. |
| `cat-tools vasp restart` | Archive outputs and use the latest structure. | Writes `.sN` backups and updates `POSCAR`, unless `--backup-only`, replaces `varc`. |
| `cat-tools vasp relaxation` | Export recorded `OUTCAR` ionic steps. | Writes extxyz, VMD files, step summaries and optional GIF. |
| `cat-tools vasp compare` | Compare periodic structures. | Prints a report; optionally writes a CSV. |
| `cat-tools vasp rescale` | Apply homogeneous lattice strains. | Writes a structure unless `--dry-run`; replacement requires `--force`. |
| `cat-tools phono3py sample` | Batch selected grid points into jobs. | Writes job scripts; submission requires `--submit`. |
| `cat-tools phono3py audit` | Check outputs and prepare recovery/collection. | Writes audit reports; can generate rerun/collection scripts and submit with `--submit`. |
| `slurm_scan` | Classify slurm logs and resource usage. | Prints a report; `--failed` writes a directory list. |
| `cat-tools ordering ...` | Generate configurations, collect static energies and build mixing/full hulls. | Writes fresh output directories/files; never runs VASP or submits jobs. See the [ordering reference](CAT/ordering/README.md). |
| `cat-tools prepare ...` | Convert structures; prepare LOBSTER, BVPA, and mass scans. | Depends on the subcommand, described below. |

Every command supports `--help` :).

---

## VASP Commands

### `potcar` — potential selection

**Purpose:** Assemble a `POTCAR` in `POSCAR` species order, with a record of the selected potential directories. Newer version of the `apgu` command.

**Usage:** `cat-tools vasp potcar`

```bash
export VASP_POTCAR_DIR=/path/to/VASP_POTCAR
cat-tools vasp potcar
cat-tools vasp potcar --sources
cat-tools vasp potcar --verbose
```

The command reads element symbols from line 6 of a VASP 5+ `POSCAR`. It prompts for the potential library and, where several variants exist, the variant for each element. GW potentials are excluded. The available library subdirectories are `PAW_PBE_64`, `PAW_LDA_64`, `PAW_PBE_54`, `PAW_LDA_54`, `USPP_GGA`, `USPP_LDA`, and `PAW_PW91`. The `POTCAR` is the authoritative source of elements identity, when making changes be sure to check your `POTCAR`s.

| Environment variable | Default | Purpose |
|:---------------------|:--------|:--------|
| `VASP_POTCAR_DIR` | `$HOME/VASP_POTCAR` | Root containing the library subdirectories. |
| `POSCAR` | `POSCAR` | Input structure. |
| `POTCAR` | `POTCAR` | Output filename in the working directory. |

The output is assembled in temporary files and checked for the expected number of potential blocks. On success, it replaces `POTCAR` and optionally writes `POTCAR.sources` when providing the `--sources` flag, `--verbose` or `-v` prints them to terminal. Missing potentials or invalid input stop generation with a non-zero status. The command always prints each dataset’s ENMAX and the maximum in eV; `--verbose` is not required for this summary.

---

### `kpoints` — reciprocal-density mesh

**Purpose:** Write a Gamma-centered mesh from the reciprocal lattice-vector lengths. Newer version of the `GenKPOINTS` command.

**Usage:** `cat-tools vasp kpoints [POSCAR] [DENSITY]`

```bash
cat-tools vasp kpoints                 # POSCAR, density 50
cat-tools vasp kpoints CONTCAR 60
```

Each mesh dimension is the rounded product of the density and reciprocal-vector length, without the `2π` factor, with a minimum of one. `KPOINTS` is written in the current directory and an existing file is replaced. This helper expects a single positive `POSCAR` scale factor; use a normalized structure if the input uses another scaling convention.

---

### `converge` — convergence job generation

**Purpose:** Create a directory per parameter value, copy the input set, adjust the `INCAR` or `KPOINTS`, and optionally submit the jobs. Newer version of the `ConvGen` command.

**Usage:** `cat-tools vasp converge AXIS (--man VALUES... | --range LOW HIGH STEP) [OPTIONS]`

Required in the working directory: `INCAR`, `POSCAR`, `POTCAR`, and `submit.sh`. `KPOINTS` is also required unless scanning `--kpoints`.

```bash
cat-tools vasp converge --encut --man 400 450 500 --no-submit
cat-tools vasp converge --encut --range 1.0 1.4 0.1 --no-submit
cat-tools vasp converge --kpoints --man 30 40 50 --no-submit
cat-tools vasp converge --ediffg --man -0.03 -0.02 -0.01 --calc ions --no-submit
```

| Option | Meaning |
|:-------|:--------|
| `--encut`, `--kpoints`, `--ediff`, `--ediffg`, `--sigma`, `--nbands` | Select one scan axis. |
| `--man VALUES...` | Explicit values; `ENCUT` values are absolute eV. |
| `--range LOW HIGH STEP` | Inclusive range; for `ENCUT`, values multiply the largest `POTCAR` ENMAX. |
| `--calc static\|ions\|full\|keep` | Set the relaxation mode; default `static`. `keep` preserves supplied relaxation tags. |
| `--set TAG=VALUE` | Additional `INCAR` override; may be repeated. |
| `--prefix PREFIX` | Override the automatic directory prefix, useful in combination with the `check` function. |
| `--no-submit` | Generate only. Without this flag, each job is submitted using `sbatch`. |
| `-v`, `--verbose` | Print setup and submission details. |

K-point scan values are reciprocal densities, not literal mesh dimensions. An `EDIFFG` scan needs an ionic relaxation mode. Existing scan directories are rejected. Each directory receives `.convgen_meta`; submitted jobs also get `.convgen_jobid` for accounting lookup. Generation/submission is sequential, so earlier directories may exist if a later step fails. Adjust slurm ressources in the `submit.sh` beforehand.

---

### `check` — convergence results

**Purpose:** Summarize calculation results and resource usage across directories matching a scan prefix. Intended to be used after the `generate` command during convergence testing. Newer version of the `ConvCheck` command.

**Usage:** `cat-tools vasp check --prefix PREFIX [OPTIONS]`

```bash
cat-tools vasp check --prefix ENCUT
cat-tools vasp check --prefix EDIFFG --relax diff
cat-tools vasp check --prefix EDIFFG --relax abs
cat-tools vasp check --prefix EDIFFG --relax ref
cat-tools vasp check --prefix EDIFFG --relax ref reference_directory
```

`diff` compares each `POSCAR` with its `CONTCAR`. `abs` reports final cell parameters. `ref` compares final structures with the strictest setting or with the supplied reference directory. Without a relaxation option, the command reports static convergence results. `--verbose` adds details.

The checker uses the generator's metadata and VASP output. If job IDs and `sacct` are available, it includes accounting data; otherwise it uses `OUTCAR` timing where possible. It does not run or resubmit failed calculations.

---

### `restart` — archive and prepare a restart

**Purpose:** Preserve current VASP output and copy the latest `CONTCAR` to `POSCAR`. Newer version of the `varc` command.

**Usage:** `cat-tools vasp restart [OPTIONS] [DIRECTORY]`

```bash
cat-tools vasp restart --backup-only
cat-tools vasp restart /path/to/interrupted_run
cat-tools vasp restart --finished-only
```

Stop the calculation before using this command: it does not query the scheduler. Interrupted runs are accepted when `CONTCAR` is usable. `--finished-only` requires normal termination in `OUTCAR`; `--backup-only` leaves `POSCAR` unchanged; `--force` accepts an unusually short, non-empty `CONTCAR`, this can lead to unexpected behavior by VASP.

Existing `OUTCAR`, `CONTCAR`, `vasprun.xml`, `out_VASP`, and `out` files are copied to a common `.sN` suffix. The original `POSCAR` is saved as `POSCAR_init` if that file does not exist. If all current outputs already match a backup, the command exits successfully without creating another backup or changing `POSCAR`. Does not submit a restart job.

---

### `relaxation` — recorded ionic-step visualization

```bash
cat-tools vasp relaxation --out relaxation-view --gif --fps 6
cat-phonovis trajectory relaxation-view/relaxation.extxyz
vmd -e relaxation-view/relaxation.vmd
```

Export complete `OUTCAR` position/force blocks with their current cells. With no input argument, join local `OUTCAR.sN` numerically followed by `OUTCAR`; if neither exists, discover numbered `relax_N` directories. Explicit files or directories are also accepted; specifying just `OUTCAR` selects only that file. Local archives/current `OUTCAR` take precedence over numbered child directories. Joined frames retain global/local step numbers and source segments; duplicate boundary geometries are kept. A matching `POSCAR` beside `OUTCAR` supplies species/order; override with `--poscar PATH`. `--stride N` keeps every Nth frame plus the last. Writes extxyz (full cell and forces), VMD XYZ/loader, and `steps.csv` (step, maximum force, volume); `--gif` adds a PhonoVis animation. Existing outputs are replaced. Requires the `structure` extra; GIF/viewing also needs `visualization`.

The trajectory is recorded optimization progress, not a phonon oscillation or physical-time simulation. Incomplete final blocks are skipped with a warning. See the [relaxation workflow](../2_Workflows/00_Common/02_Relaxation/VASP/README.md#4-submit-and-monitor) for usage and interpretation.

---

### `compare` — periodic structural comparison

**Purpose:** Report lattice changes and same-species atomic displacements between two `POSCAR`/`CONTCAR` files.

**Usage:** `cat-tools vasp compare REFERENCE OTHER [OPTIONS]`

```bash
cat-tools vasp compare POSCAR CONTCAR --top 10
cat-tools vasp compare reference/POSCAR relaxed/CONTCAR --align-origin --csv displacements.csv
```

By default, atoms are matched by species to minimize the largest periodic displacement. `--same-order` compares corresponding indices instead. `--align-origin` tries global fractional shifts; `--reference-cell a|b|average` selects the lattice used for displacement lengths. `--top` controls the displayed rows, and `--csv` writes all matches.

This is the original matching utility. It assumes the same cell setting and does not search rotated cells or supercell mappings. The separate `cat-structure compare` command retains its Pymatgen matching method; choose the one appropriate for your comparison.

---

### `rescale` — QHA and lattice-strain structures

**Purpose:** Apply a homogeneous lattice transformation while keeping fractional coordinates fixed.

**Usage:** `cat-tools vasp rescale [POSCAR] TRANSFORMATION [OPTIONS]`

```bash
cat-tools vasp rescale POSCAR --iso 5 --dry-run
cat-tools vasp rescale POSCAR --iso 5 --dir
cat-tools vasp rescale CONTCAR --ca 3 -o POSCAR_ca
```

| Transformation | Meaning |
|:---------------|:--------|
| `--iso PERCENT` | Change cell volume isotropically by the given percentage. |
| `--a`, `--b`, `--c PERCENT` | Change one lattice-vector length by the given percentage. |
| `--ca PERCENT` | Change c/a while preserving volume; a and b scale equally. |

Select exactly one transformation. The default input is `POSCAR`; the default output appends a transformation label to its filename. `--dir` writes `POSCAR` in a labeled directory, and `-o` selects an explicit path. Existing output is refused unless `--force` is supplied. `--dry-run` prints the transformation without writing.

VASP 4/5 species layouts, Direct/Cartesian coordinates, selective-dynamics flags, trailing data, and positive/negative single scale factors are supported. Output uses scale factor 1.0. This command prepares individual structures.

---

## Phono3Py Commands

### `sample` — distributed grid-point jobs

**Purpose:** Generate slurm jobs for batches of selected Phono3Py grid points.

**Usage:** `cat-tools phono3py sample MODE M1 M2 M3 N_GP [OPTIONS] [GP ...]`

```bash
cat-tools phono3py sample rta 21 21 21 10 --isotope --write-gamma
cat-tools phono3py sample jdos 21 21 21 10 --ts "300" --num-freq-points 600
cat-tools phono3py sample jdos 21 21 21 3 --ts "1 100 300 600 1000" --num-freq-points 601
cat-tools phono3py sample ise 21 21 21 5 --gp-file missing_gps.txt --dry-run
```

Modes are `rta`, `wigner`, `lbte`, `jdos`, `spf`, `ise`, `rse`, `gruneisen`, and `stp`. `N_GP` is the number of grid points per job. Explicit grid points take precedence over `--gp-file`, which takes precedence over `ir_grid_points.yaml`. Indices are validated and deduplicated.

For **JDOS**, each generated job invokes Phono3Py once per requested temperature, sequentially within its grid-point batch. This is the default on all runtime versions: Phono3Py 2.7.0 calculates multiple temperatures but its writer saves only the first. CAT also expands `--tmin/--tmax/--tstep` into individual JDOS temperatures and reads FC2 without requiring FC3 by default. Use the sampler's temperature options rather than forwarding them through `--ph3-arg`.

For **Grüneisen**, CAT generates exactly one full-mesh job. Native Phono3Py ignores grid-point selection for this calculation; batching it would cause multiple jobs to overwrite `gruneisen.hdf5`. `N_GP` remains a positional argument but grid-point selections are not used in this mode.

| Option | Meaning |
|:-------|:--------|
| `--ts "T1 T2 ..."` | Select temperatures; range options are also available. |
| `--job-dir DIR`, `--tag TAG` | Set generated-script location and names. |
| `--partition`, `--time`, `--mem`, `--cpus`, `--nodes`, `--ntasks`, `--account`, `--qos` | Set Slurm resources. |
| `--module "MODULES..."`, `--no-module` | Select the module stack or omit module loading. |
| `--write-gamma`, `--read-phonon`, `--write-pp`, `--bi`, `--sigma`, `--bmfp`, `--isotope` | Common calculation options; see `--help` for the full list. |
| `--ph3-arg ARG` | Append a raw Phono3Py argument; repeat as needed. |
| `--dry-run` | Print planned jobs without writing or submitting them. |
| `--overwrite` | Replace existing generated scripts. |
| `--submit` | Submit generated jobs; the default is generation only. |

Jobs default to `generated_jobs`, partition `long`, and a 48-hour wall time. Memory defaults depend on mode. The inherited module default is `palma/2023a foss/2023a phono3py/2.7.0`; select the stack appropriate or use `--no-module`. Resource defaults can also be supplied through the environment variables listed by `--help`.

Generated jobs check the installed Phono3Py version when they run. Wigner and NAC arguments are adjusted by the script; its Phono3Py 4 Wigner route expects the external `phono3py-wte` plugin. Check your chosen runtime's supported options before executing a generated job, be careful about potential mismatches between your local and the HPC version of Phono3Py.

The version probe imports Phono3Py with the module's `python`, avoiding the native `--version` pipeline. Set `PHONO3PY_PYTHON` to that interpreter's path when necessary. It must match the `phono3py` executable used by the job. Commands, start/finish timestamps, stdout and stderr are appended to the sampling log and mirrored to the Slurm log. A failed calculation stops the job with a nonzero status, including when output is piped through `tee`. Existing generated jobs are not automatically changed by a CAT upgrade. For Phono3Py 4, CAT also adapts force-constant I/O: the removed `--fc2/--fc3` flags are omitted, while `--no-fc2/--no-fc3` use `--no-read-fc2/--no-read-fc3`. Phono3Py 4.4.0 has a further CLI check that incorrectly requires FC3 for JDOS. For that exact version, CAT adds `--read-gamma` to bypass the unrelated FC3 requirement. The native JDOS branch exits before any gamma-file I/O; no linewidth files are read. A real FC2-only JDOS integration test covers this workaround. Use the same runtime version when continuing a scientific dataset.

---

### `audit` — results, reruns, and collection

**Purpose:** Find missing or invalid grid-point outputs and prepare rerun or final collection jobs.

**Usage:** `cat-tools phono3py audit MODE [M1 M2 M3] [OPTIONS]`

```bash
cat-tools phono3py audit rta 21 21 21 --ts "100 200 300"
cat-tools phono3py audit rta --root /path/to/run --rerun 5 --no-collect
cat-tools phono3py audit jdos --match nac
```

Supported modes are `rta`, `wigner`, `lbte`, `pp`, `jdos`, `spf`, `ise`, and `rse`. Expected grid points come from `ir_grid_points.yaml` or `--ir-grid`; the mesh can be inferred from that YAML when omitted. Supply temperature and calculation options matching the original jobs.

The command writes `qmesh_eval_MODE.csv`, `missing_gps.txt`, `invalid_gps.txt`, and `rerun_gps.txt` under the scan root. `--prefix` changes the CSV prefix; `--rerun-file` changes the combined rerun list. Existing reports are replaced.

| Option | Meaning |
|:-------|:--------|
| `--root DIR`, `--ir-grid FILE`, `--match REGEX` | Select the calculation, expected grid, and matching output files. |
| `--recursive` | Include nested output directories; use a narrowly selected root to avoid mixing runs. |
| `--num-freq-points N` | For JDOS, require exactly N frequency rows in each file. |
| `--shallow` | Check HDF5 existence/size without inspecting its contents. |
| `--rerun N_GP` | Generate jobs for missing/invalid points using the bundled sampler. |
| `--sampler EXECUTABLE` | Override the bundled sampler. |
| `--no-collect` | Skip final collection-job generation. |
| `--force-collect` | Allow collection-job generation for an incomplete dataset. |
| `--submit` | Submit generated rerun/collection jobs. |
| `--dry-run` | Suppress job-script writes/execution/submission; audit reports are still written. |

RTA/Wigner, LBTE, and `pp` modes have collection steps; the other modes do not. Complete datasets generate collection jobs by default where applicable. An incomplete dataset returns status `2`, including with `--no-collect`; completed datasets return `0`. Invocation or processing errors also return a non-zero status.

JDOS checks three finite columns, a strictly increasing frequency axis, and coverage of the requested temperatures. With `--recursive`, duplicate filenames are reported as ambiguous. Provide the original `--ts` and `--num-freq-points` for a complete check. These checks establish file integrity, not which `BORN` tensors generated the results.

JDOS recovery preserves completed temperatures. For example, if the original multi-temperature job saved only 1 K:

```bash
cat-tools phono3py audit jdos 21 21 21 \
  --ts "1 100 300 600 1000" --num-freq-points 601 \
  --rerun 3 --job-dir jobs_jdos_recovery \
  --nac --pa "1 0 0 0 1 0 0 0 1" \
  --ph3-arg=phono3py_disp.yaml \
  --cpus 4 --mem 16G --time 24:00:00 \
  --module "palma/2023a foss/2023a phono3py/2.7.0"
```

This writes `rerun_gps_t*.txt` and generates jobs in `jobs_jdos_recovery/t*/` only for missing or invalid temperatures. The audit still returns `2` until those jobs finish. Add `--submit` to submit them; otherwise inspect the scripts and submit manually. The sampler receives the audit's module, resource, job-directory and raw-argument settings. Do not start recovery for grid points still being written by another job.

RTA/Wigner collection creates a fresh `collection/MODE_bmfp_LENGTHum/` directory (or uses `--tag` in place of MODE), retaining native `kappa-m*.hdf5` filenames inside it. Only audited per-GP checkpoints are linked; an old aggregate kappa file in the scan root cannot override them. Text inputs are copied, so writing a summary YAML cannot alter the source input. Existing collection directories are refused even when `--overwrite` replaces a generated job script.

`--bmfp` selects the boundary length **for collection**. The intrinsic `gamma` and isotope rates are reusable across lengths; their original boundary metadata need not match. For example, repeat `audit wigner --bmfp 0.1`, `--bmfp 0.5`, `--bmfp 1`, `--bmfp 2`, and `--bmfp 4` against the same corrected checkpoints, with identical mesh, temperatures, force constants and NAC settings.

---

## Slurm Log Commands

### `scan` — job diagnostics and performance

**Purpose:** Summarize states and diagnostics from slurm output files, independently of scheduler accounting.

Define this helper in your shell, adjusting the checkout path:

```bash
slurm_scan() {
    python3 /path/to/PhononTools/1_Tools/General/SlurmScan.py "$@"
}
```

The script can also be copied elsewhere and run independently; CAT is not req.

**Usage:** `slurm_scan [OPTIONS]` or `slurm_scan array ARRAY_ID ... [OPTIONS]`

```bash
slurm_scan --root /path/to/calculations
slurm_scan --all --recursive --details
slurm_scan --performance
slurm_scan --diag OOM --paths
slurm_scan --failed restart_dirs.txt
slurm_scan array 123456 --root /path/to/array_logs --details
```

Directory mode examines immediate subdirectories by default and treats the highest numeric job ID in each as the latest run. Clean completed jobs are hidden unless `--all` is supplied. `--history` includes older runs; `--recursive` includes deeper directories. `--state` and `--diag` filter results; `--performance` reports completed-run resource data.

`--failed [FILE]` writes directories whose latest run is not classified as completed, including unfinished jobs. Its default is `failed_subdirs.txt`; relative output filenames are resolved against the scan root and existing files are replaced. It does not resubmit jobs. Detected calculation failures do not themselves make the scanner exit non-zero.

---

## Preparation Commands

`cat-tools prepare` and `cat-prepare` use the same functions and argument parser.

| Command after `cat-tools prepare` | Purpose and output |
|:---------------------------------|:-------------------|
| `structure convert INPUT OUTPUT` | Convert structure formats; optionally standardize with `--standardize`. Writes the selected output path. |
| `structure perturb INPUT OUTPUT` | Random atomic displacements for relaxation; `--amplitude` in Å and optional `--seed`. Existing outputs require `--force`; input replacement is refused. |
| `structure restore SITES CELLS OUTPUT` | Rebuild a crystal from canonical tables; select a frame with `--frame`. |
| `structure ase-export SITES CELLS OUTPUT` | Write canonical frames as an ASE-readable trajectory. |
| `topology POSCAR --config seed.yaml -o config.yaml` | Derive periodic framework contacts and distances into an editable LOBSTER config. |
| `lobsterin --config config.yaml -o lobsterin` | Render LOBSTER input and `lobsterin.bonding.json`; `--no-manifest` omits the sidecar. |
| `bvpa --cif FILE --cube FILE --ion Na1+` | Execute Windows BVPA from WSL and capture outputs. |
| `mass-scan CONFIG -o SCAN --species S --range 16 96 --points 17` | Prepare independent FreqBall mass configurations and run scripts. |

```bash
cat-tools prepare structure convert POSCAR structure.cif
cat-tools prepare lobsterin --config config.yaml -o lobsterin
cat-tools prepare mass-scan freqball.yaml -o mass-scan --species S --values 16 32 64
```

**LOBSTER:** Supply basis functions and interaction requests explicitly. The topology helper can derive requests from selected framework atoms; see the configuration workflow below. Generated inputs and sidecars replace the selected output files. Preserve the topology sidecar beside the input for bonding analysis.


**BVPA:** Set `--softbv-root` or `SOFTBV_ROOT`. `--conduction-path FILE` also converts output into FreqBall's path format. This command actually runs the external program. See [the BVPA procedure](../2_Workflows/99_Misc/README.md). This utility is not intended for use on the  HPC but rather assumes use of a wsl with SoftBV installed on the windows side.

**Mass scans:** Supply a working FreqBall YAML referencing existing force data. Select `--values` or `--range` with `--points`. Each mass gets a configuration and `run.sh`; `--job-template` selects a job template. `--submit-command sbatch` writes `submit_all.sh` but does not submit it. Existing point configurations/scripts are refused unless `--force` is supplied. See [the mass-scan procedure](../2_Workflows/10_LatticeDynamics/README.md).

Reusable Python interfaces and the package layout are described in [CAT/README.md](CAT/README.md).

### Perturbing a structure before relaxation

```bash
cat-prepare structure perturb POSCAR POSCAR.perturbed --amplitude 0.02 --seed 42
```

Each mobile atom moves exactly `--amplitude` Å (default 0.02) in a random direction, uniform within its allowed Cartesian subspace. This is a displacement length, not a per-component standard deviation. VASP selective-dynamics flags constrain direct coordinates, including in skew cells; fully fixed atoms stay fixed. An entirely fixed structure is rejected. The cell, species order and site properties are retained in memory; output formats may not support all properties. `POSCAR` retains selective dynamics. Coordinates are not wrapped or symmetry-standardized.

Use a non-negative `--seed` to reproduce a realization; omitting it gives a fresh random realization. `--input-format` and `--output-format` override format detection. Install the `structure` extra. The command writes only the requested structure and does not submit jobs or modify `INCAR`. Existing outputs require `--force`, and the input cannot be overwritten. Inspect distances before running a relaxation: there is no minimum-distance rejection. This prepares a new relaxation geometry, not an MD restart with preserved velocity/predictor history. See the [relaxation workflow](../2_Workflows/00_Common/02_Relaxation/VASP/README.md).

## LOBSTER Configuration and Framework Detection

**Use:** Prepare `lobsterin` for your system, chemistry agnostic. Start with [the example seed](CAT/examples/lobsterin.seed.yaml), adapting its species, basis functions, and distance window to your calculation:

```bash
cat-prepare topology POSCAR --config seed.yaml -o config.yaml
cat-prepare lobsterin --config config.yaml -o lobsterin
```

The first command uses `POSCAR` coordinates and periodic boundary conditions to find selected bonds, derive distances, and write atom numbers and cell translations. Review the generated YAML, then run the second command to render it. You can also write a complete config by hand and run only `lobsterin`.

```yaml
system:
  structure: POSCAR
  StartE: -20
  EndeE: 15
  basisSet: pbeVaspFit2015
  basisFunctions:
    P: 3s 3p
    S: 3s 3p
  save: true
  write: false

Generator:
  - types: [P, S]
    distance: [1.8, 2.4]
    orbitalWise: true

Atom:
  - atoms: [1, 5]
    cell: [0, -1, 0]
    orbitalWise: true
```

This example requests:

```text
cohpGenerator from 1.8 to 2.4 type P type S orbitalWise cohpBetween atom 1 atom 5 cell 0 -1 0 orbitalWise
```

| Setting | Meaning |
| --- | --- |
| `system.structure` | `POSCAR`/`CONTCAR` path relative to the YAML; defaults to `POSCAR`. A positional structure argument overrides it. |
| `system.StartE`, `system.EndeE` | COHP energy window in eV; defaults to -20 and 15. `EndE` is also accepted. |
| `system.basisSet` | Required LOBSTER basis-set name. |
| `system.basisFunctions` | Required mapping with exactly one orbital string per element type in `POSCAR`, including atoms outside the framework. Choose these to match your potentials. |
| `system.save` | Emit `saveProjectionToFile` when true; default false. |
| `system.write` | Emit `writeBasisFunctions` when true; default false. |
| `system.cohp`, `coop`, `cobi`, `dos` | Optional booleans; false emits the corresponding skip keyword. Default true. |
| `Generator` | List of type-pair requests. Each entry has `types: [A, B]`, optional `distance: [min, max]` in angstrom, and optional `orbitalWise` (default false). |
| `Atom` | List of explicit requests. Each has `atoms: [i, j]`, optional `cell: [a, b, c]` (default zero), and optional `orbitalWise` (default false). |

Atom numbers are **1-based `POSCAR` indices**. `cell` translates the second atom relative to the first, in integer lattice vectors. A generated `Atom.distance` records the reconstructed distance in angstrom and is checked against `POSCAR` when rendering, so changed structures require regenerating the configuration.

To derive contacts automatically, add an explicit `framework` section to the seed:

```yaml
framework:
  centers: {types: [P]}       # Or: {atoms: [1, 2]}
  ligands: {types: [S]}       # Or: {atoms: [3, 4, 5, 6]}
  distance: [0, 2.4]
  # coordination: 4         # Optional; otherwise keep all bonds in the window.
  padding: 0.05
  orbitalWise: true
  interactions: [bonds, intra]
```

`types` and `atoms` select their union when both are supplied. Both center and ligand selections must be nonempty. The required distance window defines what counts as a center–ligand bond. Optional `coordination` keeps the nearest N ligands within it and fails if too few exist; near-tied extra neighbors produce a warning. This supports different coordination numbers without assuming tetrahedra.

`bonds` adds center–ligand requests; `intra` adds ligand–ligand pairs sharing a center, including pairs across cell boundaries. These are the default selections. The helper also identifies connected bonded components and distinguishes finite polyanions from periodic networks; this information is saved in the LOBSTER sidecar. Connectivity follows your supplied bond criteria, not an inferred oxidation state. It does not automatically request every pair in a connected network or between separate polyanions.

If a `Generator` entry omits `distance`, the helper derives its range from matching selected center–ligand bonds, expanding both ends by `padding` angstrom (default 0.05). Without matching framework bonds, provide the range explicitly. Type-based LOBSTER generators apply to **all pairs of those types in that range**, including sites outside your framework selection. For strictly site-specific requests, use `Generator: []`; topology then writes explicit `Atom` entries for selected contacts. Contacts already covered by a generator or an explicit request are not added again.

The renderer uses the listed `Generator` and `Atom` requests; it does not silently add contacts from `framework`. The old chemistry-specific defaults and generation flags have been replaced by this YAML interface. Analysis continues to use the shared topology code and the generated `lobsterin.bonding.json` sidecar.

### Inter-polyanion and Inter-chain Contacts

Add `inter` to the topology selections to request ligand–ligand contacts between **different connected components**. This works for finite polyanions, molecular copies, polymer chains, and extended networks. It does not assume that a component has a particular charge.

```yaml
framework:
  centers: {types: [P]}
  ligands: {types: [S]}
  distance: [0, 2.4]          # Defines the bonded framework itself.
  interactions: [bonds, intra, inter]
  orbitalWise: true
  inter:
    selection: nearest
    cutoff: 4.5             # Search radius for contacts BETWEEN components, in Å.
    tolerance: 0.03
    orbitalWise: true
```

Run the same two commands after editing the seed:

```bash
cat-prepare topology POSCAR --config seed.yaml -o config.yaml
cat-prepare lobsterin --config config.yaml -o lobsterin
```

| `framework.inter` setting | Meaning |
| --- | --- |
| `selection: nearest` | Default. Keep the shortest ligand–ligand contact for **each neighboring component pair**, plus contacts within `tolerance` of that minimum. |
| `selection: all` | Keep every inter-component ligand contact within `cutoff`. |
| `cutoff` | Required positive search radius in Å. Contacts beyond it are never included, even within the tie tolerance. |
| `tolerance` | Nonnegative distance tolerance in Å; default 0.03. Applies to `nearest`; zero keeps exact numerical ties. |
| `orbitalWise` | Boolean; inherits `framework.orbitalWise` unless supplied here. |

The framework's `distance` window builds the bond network; the separate inter `cutoff` searches between its components. Sharing a coordination center is not the criterion for membership: centers joined through shared ligands belong to one component, so their ligand pairs are excluded from `inter` even if they share no individual center. Only selected ligands connected to the framework participate; unattached atoms are not treated as additional polyanions.

For finite molecules, a periodic image is a separate component instance. For a polymer, translating along its bonded repeat direction stays within the same chain and is excluded. Translating to a neighboring chain is included. This distinction also handles oblique chains, sheets, and disconnected interwoven networks: the code uses the exact integer translation subgroup, not merely whether the structure is periodic. Distinct images may share atom numbers but have different `cell` vectors. A single fully connected network has no inter-component contacts with itself.

For periodic chains or sheets, equivalent component pairs related by bonded repeat translations are grouped together when finding the minimum. Reversed contacts and exact duplicate requests are removed; tied contacts remain separate atom requests. Use `selection: all` when you want longer inter-chain contacts as well as the nearest.

**Polymer example:** centers and ligands may overlap in element type, so a carbon backbone can be selected with `centers: {types: [C]}` and `ligands: {types: [C, H]}`. See [the polymer seed](CAT/examples/lobsterin.polymer.seed.yaml). Adapt the bond window to your system and include all species needed to represent its connectivity. Selecting only a fragment of a chain, or omitting a bridging species, can split one physical molecule into separate detected components.

The generated `Atom` rows contain the selected indices, periodic cells, measured distances, and orbital flags. For strict component-specific requests, keep `Generator: []`: a type-based generator can additionally request contacts within a component. Existing explicit requests are preserved; topology does not delete them. Generate a fresh config from the seed when changing selection criteria, then render that config. Rendering alone does not regenerate its contact list.

The sidecar records the inter-contact settings and each component's periodic dimension (0: finite, 1: chain, 2: sheet, 3: network) and translation basis. With these settings present, bonding analysis uses the same inter-component selection and distinguishes `X-X_same_component` contacts from `X-X_inter` where no shared coordination center exists. Older sidecars without inter-component settings retain their previous inter-polyhedron classification.
