# Setting Up a VASP Calculation

**Use:** Assemble checked VASP inputs and verify the execution setup before production.

**Files:** [INCAR.example](INCAR.example) · [KPOINTS.example](KPOINTS.example) · [job.sh.example](job.sh.example).

**Prerequisites:** Access to VASP and the potential datasets you will use, and a compatible execution environment. CAT's structure conversions need Python 3.12+ and its `structure` extra. Preparation can run without VASP; calculation requires your local installation or an allocation on the machine where VASP is available.

**Scope:** VASP 6; execution examples use CPU/MPI and Slurm where a scheduler is needed. Adapt them to your launcher and hardware. The example `INCAR` illustrates nonmagnetic bulk PBEsol; choose the physical model for your own material first.

**Navigation:** [Environment](#1-set-up-access-and-the-environment) ·[Inputs](#2-input-files) · [Structures](#3-prepare-the-structure) · [POTCAR](#4-assemble-and-inspect-potcar) · [INCAR](#5-choose-incar-settings-by-purpose) · [KPOINTS](#6-generate-an-initial-kpoints-mesh) · [Parallelization](#7-match-parallelization-to-the-node-layout) · [Submission](#8-write-jobsh-and-submit) · [Outputs](#9-read-the-outputs-and-judge-the-run)

---

## 1. Set Up Access and the Environment

### Make VASP available where you will calculate

Use your existing VASP installation or follow your site's access and environment instructions. A cluster is optional. On a module-based installation, for example:

```bash
module load YOUR_VASP_MODULE
command -v vasp_std
```

Replace `YOUR_VASP_MODULE` with the installed module; skip `module load` when your executable is already available. Repeat the required environment setup in your job script if using a scheduler. The [PALMA site notes](../../../../1_Tools/CAT/vasp-site-notes.md) retain this group's access and module examples separately.

### Make CAT available

If CAT is already available, continue to the next step. Otherwise, from your repository checkout (replace the path):

```bash
cd /path/to/PhononTools
uv sync --extra structure
source .venv/bin/activate
cat-tools --help
cat-prepare structure convert --help
```

Alternatively, use `python -m pip install -e '.[structure]'` in an existing Python 3.12+ environment. The full [Tools setup](../../../../1_Tools/README.md#setup) also covers analysis extras. Keep the Python preparation environment separate from the compiler/MPI choices required by your VASP executable.

### Locate the potential library

Use the potential release and variants required by your study. For reproducing an existing calculation, retain its datasets. CAT assembles local `POTCAR` files; it does not download them. The [potential overview](https://vasp.at/wiki/Available_pseudopotentials) describes releases and variants.

The CAT builder expects a root directory containing library directories such as `PAW_PBE_64`, then element/variant directories containing individual `POTCAR` files. Keep the library wherever it is accessible to you. After checking `POSCAR`, you can supply the location for just this command:

```bash
VASP_POTCAR_DIR=/path/to/potential-library cat-tools vasp potcar --sources --verbose
```

Replace the path with the root containing `PAW_PBE_64` (or another supported library), not an individual element folder. If already configured, simply run `cat-tools vasp potcar --sources --verbose`. Persisting `VASP_POTCAR_DIR` in your shell configuration is optional. Unlike a workflow convenience variable, this is an actual CAT setting; VASP reads only the assembled calculation-specific `POTCAR`. See the [builder reference](../../../../1_Tools/README.md#potcar--potential-selection) for supported library layouts and settings.

For convenient monitoring, install the optional [Bash helpers](../../../../0_bashrc.d/README.md#setup). These are independent of the VASP executable and the Python environment.

---

## 2. Input Files

Use a separate working directory for each calculation:

```text
calculation/
├── POSCAR          # Lattice, species, atom counts, coordinates, optional constraints.
├── POTCAR          # Concatenated PAW datasets in POSCAR species-block order.
├── INCAR           # Physical model, numerical settings, task, and output choices.
├── KPOINTS         # Brillouin-zone sampling.
└── job.sh          # Scheduler resources, software environment, and launch command.
```

VASP reads `INCAR`, not the template names `INCAR.static` or `INCAR.d`. This workflow uses explicit `KPOINTS` rather than the alternative `KSPACING` interface. Keep each calculation in its own directory, including separate benchmark trials; mixing calculation inputs can result in unexpected behavior by VASP or, more often than not, a failed calculation run.

When revising a setup, check the downstream dependencies:

| Change | Recheck |
|:-------|:--------|
| Atom/species order | `POTCAR` order, `MAGMOM` and constraints, and all atom-index-based analysis inputs. |
| Cell, supercell, or symmetry | K-point mesh, irreducible k-point count, `KPAR` balance, and restart compatibility. |
| Potential variant or functional | Valence/band requirements, numerical convergence, and consistency of reference energies. |
| MPI layout | Actual `NBANDS`, timings, and peak memory; nominally identical INCARs can acquire different band counts. |
| Restart strategy | Whether the calculation starts from scratch, reuses orbitals/charge, or fixes a charge density for a non-SCF task. |

---

## 3. Prepare the Structure

### Start from an experimental CIF

A publication's supporting information or an experimental database such as the [Crystallography Open Database](https://www.crystallography.net/cod/) can provide a starting CIF. A refined experimental structure is not automatically a zero-temperature DFT equilibrium structure, relaxation is always necessary.

Check occupancy and disorder before conversion. A partially occupied or mixed site needs an explicit ordered model, often in a supercell; exporting a file does not make that modeling choice for you. Inspect atom counts, composition, lattice, interatomic separations, and any missing species. If you leave partial occupancies during conversion they will be filled fully throwing off your stoichiometry.

With CAT:

```bash
cat-prepare structure convert input.cif POSCAR --output-format poscar
```

By default CAT does not request symmetry standardization. Use `--standardize primitive` or `--standardize conventional` only when that cell change is intended. Avoid symmetry refinement of an MD snapshot or a distorted structure, it will undo the distortion.

### VESTA conversion

Inspect the CIF in VESTA, then use **File → Export Data → VASP/`POSCAR`**. Reopen the export to verify the cell, occupancies, species blocks, and atom count; displayed periodic neighbors are not automatically the intended simulation supercell. [VESTA export documentation](https://jp-minerals.org/vesta/en/doc/VESTAch4.html)

### Write an existing Pymatgen Structure

```python
from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Poscar

structure = Structure.from_file("input.cif")
# If your script already has a Structure object, start from here
if not structure.is_ordered:
    raise ValueError("Choose an explicit ordered model before writing POSCAR")
Poscar(structure, sort_structure=True).write_file("POSCAR")
```

Here `sort_structure=True` groups/reorders sites; any existing atom indices must be updated. For a structure whose ordering already matters, use `sort_structure=False` and inspect the resulting species blocks. Regenerate `POTCAR` from the final `POSCAR` rather than retaining one from an earlier ordering. [Pymatgen structure I/O](https://pymatgen.org/pymatgen.io.vasp.html)

### Write an ASE Atoms object or an MD snapshot

For an existing ASE object, conversion is similarly direct:

```python
from ase.io import read, write

atoms = read("input.cif")
# Replace the preceding line with your existing Atoms object if appropriate.
write("POSCAR", atoms, format="vasp", direct=True, vasp5=True, sort=True)
```

To select one saved frame from an ASE trajectory:

```python
from ase.io import read, write

atoms = read("trajectory.traj", index=-1)  # Last saved frame.
if atoms.cell.rank != 3:
    raise ValueError("Supply the complete simulation cell before exporting to VASP")
write("POSCAR", atoms, format="vasp", direct=True, vasp5=True, sort=True)
```

The same approach can read another supported trajectory format. Preserve the cell of the selected frame, especially after variable-cell MD. A plain XYZ file often lacks this information; obtain it from the original simulation rather than guessing. LAMMPS dumps may also require an explicit atom-type-to-element mapping, details may follow.

`sort=True` changes the atom order. Check constraints after conversion and regenerate all atom-index-dependent inputs. For an MD continuation, velocities and thermostat state require additional handling; this example prepares a structural snapshot for a static calculation or subsequent relaxation. VASP uses a periodic cell even when an ASE object's `pbc` flags are false, so molecular/slab models also need an appropriate vacuum and sampling setup. [ASE's VASP writer](https://docs.ase-lib.org/_modules/ase/io/vasp.html)

### Start from a previously calculated structure

| Source | Useful starting point | What to retain/check |
|:-------|:----------------------|:---------------------|
| [Materials Project](https://docs.materialsproject.org/apps/explorer-apps/materials-explorer/tutorial) | Search a composition and export a calculated structure as CIF/`POSCAR`. | Material ID, selected polymorph, calculation method, and whether you exported the relaxed structure. |
| [OQMD](https://oqmd.org/documentation/structures) | DFT structures from experimental entries and computational prototypes. | Entry ID, structural origin, and relaxation/calculation metadata. |
| [AFLOW](https://aflow.org/documentation/) | Calculated materials and structural prototypes. | Entry identifier and whether the geometry is relaxed or an unrelaxed prototype. |

Downloaded energies and geometries may use different functionals, potentials, magnetic states, or numerical settings from yours. Use the geometry as a starting point and validate it with your own settings before making energy comparisons.

Start in the calculation directory you chose, then copy the accepted starting structure if it is elsewhere:

```bash
cp /path/to/starting/POSCAR POSCAR
```

Skip this copy when your checked `POSCAR` is already in place. Retain the original source and species order. Optional [material examples](../../../../1_Tools/CAT/examples/na3ps4/README.md) live in tools; they are not prerequisites for this workflow.

---

## 4. Assemble and Inspect `POTCAR`

After finalizing `POSCAR`, run:

```bash
cat-tools vasp potcar --verbose
```

Select the intended potential family/release, then a variant for each species. CAT assembles the datasets in `POSCAR` species-block order and checks the number of potential blocks before replacing `POTCAR`. It excludes GW-named variants from this interactive selection; a specialized GW workflow needs its own potential choices.

A quick inspection is enough here:

```bash
sed -n '6,7p' POSCAR
grep -E 'TITEL|VRHFIN|ZVAL|ENMAX' POTCAR
```

- `TITEL` identifies the dataset/variant: check the sequence against `POSCAR`.
- `ZVAL` is the valence-electron count represented by that dataset; including semicore   states changes the electronic problem and cost.
- `ENMAX` supplies a cutoff reference for that potential. Inspect the largest value  across the assembled datasets before setting `ENCUT`.

Use consistent potential choices across calculations you intend to compare. A variant such as `_pv` or `_sv` is a modeling choice, not an automatic improvement for every calculation. Keep the datasets intact; select the functional in `INCAR` rather than editing `POTCAR` headers. The detailed [POTCAR references](#references) explain variant selection and dataset contents. CAT also prints out the different ENMAX values for the chosen potentials for your convinience.

**Species-order mistakes can still produce a running calculation, be sure to check this when reusing POTCARs**

---

## 5. Choose `INCAR` Settings by Purpose

Keep one active assignment per tag. `OUTCAR` is the check on what VASP actually interpreted. Separate changes to the physical model, numerical accuracy, algorithm, and output policy when comparing calculations.

### Numerical accuracy and electronic convergence

| Tags | What you decide | Practical interpretation |
|:-----|:----------------|:-------------------------|
| `ENCUT` | Plane-wave cutoff, in eV. | Set it explicitly after inspecting `POTCAR`, then converge to the target observable. A larger basis costs time and memory. |
| `PREC` | Numerical grid/projector precision. | `Accurate` is a useful starting point for force/stress work. When using `LREAL` make sure to check ROPT explicitly. |
| `EDIFF` | Electronic self-consistency threshold, in eV. | `1E-6` is a good starting value. Tighten and test for demanding force calculations. |
| `ALGO` | Electronic minimization algorithm. | Use `Normal` as the reference here; compare alternative algorithms on a stable calculation, including iteration count rather than time per iteration alone. |
| `NELM` | Maximum electronic iterations per ionic step. | Reaching the limit is a reason to inspect convergence; increasing it does not repair persistent oscillations. |
| `ISMEAR`, `SIGMA` | Occupation scheme and broadening in eV. | Gaussian `ISMEAR=0` with a small width is a useful initial choice when the electronic character is not yet established. Metals need explicit smearing/k-point checks. Semiconductors may pref tetrahedron methods. |
| `NBANDS` | Number of electronic bands. | Inspect the actual value and occupations. Empty-state or projection workflows may need more bands than a basic ground-state calculation. More bands increase cost. |

`ISMEAR=-5` is a tetrahedron integration choice for suitable meshes and tasks, not a general replacement for every relaxation or Γ-only calculation. The [smearing guide](https://vasp.at/wiki/Smearing_technique) explains when it is appropriate. See also [ENCUT](https://vasp.at/wiki/ENCUT), [PREC](https://vasp.at/wiki/PREC), and [EDIFF](https://vasp.at/wiki/EDIFF) for the distinct roles of these settings.

### Physical model and optional refinements

| Tags | What to decide before copying them |
|:-----|:----------------------------------|
| `GGA` / `METAGGA` / `LHFCALC` | Which exchange-correlation approximation defines the calculation? The example below uses PBEsol (`GGA=PS`) with PBE potentials. Meta-`GGA` and hybrid setups need their own requirements and convergence checks. |
| `ISPIN`, `MAGMOM` | Is spin polarization needed, and which initial magnetic configurations will you test? `ISPIN=1` assumes a non-spin-polarized calculation; initial moments are not a proof of the final magnetic ground state. |
| `LDAU*`, `IVDW`, `LSORBIT` | DFT+U, dispersion, and spin-orbit coupling are physical-model choices. Add them for a stated reason and a compatible workflow, not because another `INCAR` contains them. |
| `LREAL` | Reciprocal-space projectors (`.FALSE.`) or optimized real-space projectors (`Auto`)? The latter can help large cells; validate its accuracy using ROPT against a reciprocal-space reference. |
| `ROPT` | Real-space projector optimization accuracy, with one entry per `POTCAR` species block. Relevant for `LREAL=Auto`; no ROPT tuning is needed with `LREAL=.FALSE.`. |
| `LASPH` | Whether nonspherical contributions inside PAW spheres are included. Retain an intentional, consistent choice, particularly for methods/materials where these contributions matter. |
| `ADDGRID` | An additional support grid can affect force noise. Test whether it helps your system rather than treating it as a universal accuracy switch. |
| `ISYM` | Which symmetry treatment is appropriate? Imposed distortions, magnetic configurations, and some downstream workflows need their own symmetry choices. |

Check whether inherited `NELMIN`, `NELMDL` and mixing overrides are needed. 

For difficult SCF cases, distinguish an iteration limit from charge/magnetic oscillations or unsuitable occupations before changing several controls together. Relevant tag references are collected [below](#references).

### Validate `LREAL` and ROPT together

For `LREAL=Auto`, the default ROPT is `-5E-4` with `PREC=Normal` and `-2.5E-4` with `PREC=Accurate`. Setting ROPT explicitly overrides that `PREC`-dependent choice. For an input with three species blocks, an illustrative trial is (adjust the number and order of entries to your `POTCAR`):

```ini
LREAL = Auto
PREC  = Accurate
ROPT  = -2.5E-4 -2.5E-4 -2.5E-4
```

In this small-magnitude precision regime, **smaller absolute values are tighter**: try `-1E-4` for each species next, then tighten further if needed. ROPT is an approximate projector-accuracy target, not a guaranteed total-energy error bound. These values are for `Auto`; avoid transferring settings from the older `On` scheme.

Keep geometry, potentials, `ENCUT`, k mesh, and electronic settings fixed. Compare energies, forces, stress, and timing against `LREAL=.FALSE.`, then retain the least expensive ROPT meeting your tolerances. Inspect `OUTCAR`'s real-space projector optimization report and revalidate after changing potentials or cutoff. Keep a consistent validated projection setup across related energy comparisons. [ROPT reference](https://vasp.at/wiki/ROPT)

### Calculation type and restarts

| Tags | Role |
|:-----|:-----|
| `NSW`, `IBRION` | `NSW=0`, `IBRION=-1` select the static example below. Relaxation requires an ionic algorithm and a positive step budget; MD is a separate setup. |
| `ISIF` | Which ionic/cell degrees of freedom can relax, and whether stress is evaluated. In ordinary relaxation, `2` keeps the cell fixed while moving ions; `3` also permits cell changes. |
| `EDIFFG` | Ionic stopping criterion. A negative value specifies a force threshold in eV/Å; a positive value specifies an energy-change criterion. It does not affect a static calculation. |
| `ISTART`, `ICHARG` | Where initial orbitals and charge density come from. The fresh-start example uses `0` and `2`. Reading an existing `WAVECAR` or `CHGCAR` requires an intentional, compatible restart setup. |
| `SYSTEM`, `NWRITE` | A human-readable run label and `OUTCAR` verbosity; `NWRITE=2` is suitable for this walkthrough, but consumes more memory on disk. |

For a restart, distinguish reusing a compatible orbital/charge guess from holding charge fixed: for example, `ICHARG=11` belongs to a non-self-consistent calculation, not a generic SCF restart. Keep starting-file choices identical across performance trials. See [ISTART](https://vasp.at/wiki/ISTART) and [ICHARG](https://vasp.at/wiki/ICHARG).

Preserve the previous run before restarting. [CAT's restart utility](../../../../1_Tools/README.md#restart--archive-and-prepare-a-restart) archives outputs and can promote `CONTCAR` to `POSCAR`.

### Static baseline example

For nonmagnetic bulk PBEsol, this provides an editable static baseline. **520 eV is illustrative:** replace it if necessary after checking the selected potentials, and establish the final cutoff in the convergence workflow. Other physical models and tasks require their own settings.

```ini
SYSTEM  = Initial_static
GGA     = PS

ENCUT   = 520
PREC    = Accurate
EDIFF   = 1E-6
ALGO    = Normal
NELM    = 100
ISMEAR  = 0
SIGMA   = 0.04
LREAL   = .FALSE.
LASPH   = .TRUE.
ADDGRID = .TRUE.

ISTART  = 0
ICHARG  = 2
NSW     = 0
IBRION  = -1
ISIF    = 2

NWRITE  = 2
LWAVE   = .FALSE.
LCHARG  = .FALSE.

KPAR    = 1
NCORE   = 1
```

`NCORE=1` is a small-cell starting point for the four-rank CPU/MPI example below, not a material parameter. Select it together with the resources requested in `job.sh`. The two output flags avoid large restart files for this setup check; enable them when a later calculation needs those files.

---

## 6. Generate an Initial `KPOINTS` Mesh

```bash
cat-tools vasp kpoints POSCAR 50
cat KPOINTS
```

CAT writes a Γ-centered mesh, replacing an existing `KPOINTS` file. The argument `50` is a reciprocal-density criterion, not a request for 50 *k* points. For an orthogonal cell it gives approximately `Ni ≈ 50 / Li`, with `Li` in Å; the implementation uses reciprocal-vector lengths for a general cell and rounds each dimension to at least 1.

Inspect the mesh produced for your actual cell. This is an initial sampling trial, not a converged prescription. The symmetry-reduced *k*-point count used for computation appears in `IBZKPT`/`OUTCAR` and can be much smaller than the full mesh product.

The helper assumes one positive `POSCAR` scale factor. Normalize other scale formats through a structure writer before using it. Inspect meshes for strongly skewed cells, slabs, or molecules; the helper does not choose a dimensionality-specific sampling policy. Converge the mesh for the observable and geometry of interest.

For large supercells, see the optional [Gamma-only calculation guide](gamma-only.md).

---

## 7. Match Parallelization to the Node Layout

### MPI ranks, threads, and VASP work groups

A node is a machine with physical cores, memory, and possibly multiple sockets or NUMA domains. An MPI rank is a process; OpenMP threads share work within that process. Slurm's allocated CPUs can refer to cores or hardware threads depending on site configuration. Check the compute-node layout with the cluster documentation or, inside an allocation, `lscpu` and `numactl --hardware` if available.

For a standard **CPU/MPI-only** calculation with one image and no additional algorithm-specific partitioning:

```text
R = nodes × MPI ranks per node
ranks per k-point group = R / KPAR
band-parallel groups = R / (KPAR × NCORE)
```

| Tag | Meaning | How to choose a first trial |
|:----|:--------|:---------------------------|
| `KPAR` | Number of *k*-point groups. | Must divide the total MPI ranks. Use the actual irreducible *k*-point count; balanced groups work best when that count is divisible by `KPAR`. Use 1 for a single *k* point. |
| `NCORE` | MPI ranks cooperating on one band's plane-wave/FFT work. | Must divide the ranks in a *k*-point group. Try factors of the physical cores in a socket/NUMA domain and benchmark. |
| `NPAR` | Number of band-parallel groups, the older alternative to `NCORE`. | For this setup, `NPAR = R / (KPAR × NCORE)`. Prefer setting `NCORE` and omit `NPAR`; setting both lets `NPAR` take precedence. |

`NCORE` is not a node count or an OpenMP-thread count. These relations and the mutual exclusion of `NCORE`/`NPAR` are described in the [NCORE reference](https://vasp.at/wiki/NCORE). Additional image or algorithm parallelization changes the available ranks.

**Worked resource example:** suppose an allocation has two nodes with 32 physical cores each and uses one MPI rank per core. With 64 ranks and 8 irreducible *k* points, `KPAR=2`, `NCORE=4` gives 32 ranks per *k*-point group and 8 band-parallel groups. A trial with `KPAR=4`, `NCORE=4` gives 16 ranks per *k*-point group and 4 band groups. These are valid configurations to compare, not a guarantee that either is fastest. Rank placement determines whether the groups actually align with node/NUMA boundaries.

For small cells, benchmark a modest allocation before increasing the node count. More ranks can add communication and increase the band count without reducing wall time. Inspect the actual `NBANDS` chosen by VASP, which may be rounded up for the parallel layout. [NBANDS and parallelization](https://vasp.at/wiki/NBANDS)

Hybrid MPI/OpenMP and GPU builds need different tuning: do not apply this MPI-only `NCORE` recipe unchanged. `NCORE` is reset to 1 in the threaded/GPU paths; use the build-specific thread/rank guidance. On a hybrid CPU job, `--cpus-per-task` reserves threads for each MPI rank and must agree with `OMP_NUM_THREADS`. [Parallelization guide](https://vasp.at/wiki/Optimizing_the_parallelization)

### Benchmark the documented 36-core layout

The supplied PALMA experience with `NCORE=18`, `NSIM=18` concerns large metallic systems (>250 atoms) on 36-core nodes with `VASP/5.4.4`-vtst. It is a candidate to measure, not a default for small cells, relaxations or response calculations. Use the [small/medium/large benchmarking exercise](benchmarking.md) to compare rank counts, `NCORE` and `NSIM` while tracking actual `NBANDS` and result agreement.

**Strongly recommended:** Record each calculation in the group's shared [benchmark log](benchmark-log.md), accessible to all colleagues for reading, contributions and discussion. Add one simple table row with input settings, resources, restart-file presence and measured timings/memory; include your name and note failures. Follow the [shared logging procedure](benchmarking.md#5-measure-repeat-and-choose) and link the shared table from your production calculation notes.

### Memory is a separate constraint

Memory demand changes with cell size, cutoff/FFT grids, bands, spin, and the method. Increasing `KPAR` can reduce communication but increases memory pressure because *k*-point parallelization does not distribute all associated storage. If memory is tight, test a smaller `KPAR`. Increasing `NCORE` can reduce replicated projector memory in CPU/MPI calculations, but the fastest setting still needs measurement. [VASP memory guidance](https://vasp.at/wiki/Not_enough_memory)

A node's memory must accommodate all ranks running there, not just one rank's reported usage. `OUTCAR`'s rank-0 memory breakdown is diagnostic, not the whole job's peak memory. Fewer ranks per node can provide more memory per rank if the allocation actually grants the corresponding memory; adding nodes is not a universal fix for replicated arrays. Request memory explicitly and inspect scheduler accounting. `MAXMEM` is only read by selected VASP routines and does not enforce a global memory cap; it is therefore omitted from these ordinary SCF and response examples. See [MAXMEM](https://vasp.at/wiki/MAXMEM).

---

## 8. Write `job.sh` and Submit

Use the execution method available where you calculate. For a local CPU/MPI installation, an example is:

```bash
mpirun -np RANKS vasp_std > vasp.out 2> vasp.err
```

Replace `RANKS` with the tested process count and use the launcher required by your build. Run from the directory containing the four VASP inputs. On a cluster, use allocated compute resources rather than a login node.

If using **Slurm**, use your reviewed job script or copy the linked [CPU/MPI example](job.sh.example):

```bash
cp /path/to/PhononTools/2_Workflows/00_Common/00_SettingUp/VASP/job.sh.example job.sh
```

Edit `job.sh` before submitting. Its PALMA module/account values and small pilot allocation are examples; they do not select appropriate resources for your system. Review these items for your installation:

| Item | Decision |
| --- | --- |
| Modules / executable / MPI launcher | Use a compatible installed stack; configure any required runtime library paths. |
| Nodes, tasks and threads | Match hardware, the VASP build and `INCAR` parallelization; test a representative workload. |
| Memory and wall time | Request measured needs plus headroom. Per-node and per-CPU memory requests have different meanings. |
| Partition, account and optional QoS | Use values valid for your allocation. |
| Exclusive allocation / requeue | Follow site policy and preserve attributable outputs on retries. |
| Output paths | Keep scheduler and VASP logs; avoid rerunning over an existing result. |

Slurm directives require literal values, not shell expansion. Runtime variables such as thread limits belong in the script body when the executable needs them; no exported project-directory variable is required. See the [sbatch reference](https://slurm.schedmd.com/sbatch.html).

After reviewing inputs and resources, submit from the calculation directory:

```bash
sbatch job.sh
squeue -u "$USER"
tail -f vasp.out
```

Use either the local launcher or submission route for a given run. Keep the chosen layout fixed during numerical convergence comparisons. If the optional [Bash helpers](../../../../0_bashrc.d/README.md) are installed, they provide extra monitoring; ordinary scheduler commands work without them.

---

## 9. Read the Outputs and Judge the Run

### What is written where?

| File | Look here for | Related choices |
|:-----|:--------------|:----------------|
| `vasp.out` | Startup, electronic progress, and warnings printed to standard output. | The filename comes from the redirect in `job.sh`, not an `INCAR` tag. |
| `slurm-*.out`, `slurm-*.err` | Environment/launcher output, scheduler errors, time-limit or OOM messages. | The `job.sh` log directives and resource limits. |
| `OUTCAR` | Interpreted settings, energies, forces, stress, magnetization, timings, and memory information. | `NWRITE` controls detail; inspect actual `ENCUT`, `NBANDS`, spin, and parallel layout. |
| `OSZICAR` | Compact electronic-iteration and ionic-step summaries. | `EDIFF`, `NELM`, `ALGO`; ionic progress also depends on `NSW`/`IBRION`. |
| `CONTCAR` | Latest written geometry. | Relaxation/MD and ionic constraints. Its existence does not prove that the run was fully successful. |
| `vasprun.xml` | Structured settings and results for parsers and analysis. | Run type and output choices; interrupted files may be incomplete, in principle this file contains all information about a calculation. |
| `IBZKPT` | Irreducible k points and weights. | `KPOINTS` and symmetry; useful when choosing `KPAR`. |
| `WAVECAR` | Orbitals for compatible restarts or downstream tools. | `LWAVE`; loading is controlled by `ISTART`. |
| `CHGCAR` | Charge density and PAW information for suitable restarts/analysis. | `LCHARG`; loading/charge treatment is controlled by `ICHARG`. |
| `DOSCAR`, `EIGENVAL` | DOS and band eigenvalue information. | Mesh, smearing, `NBANDS` and DOS settings determine whether these are suitable for the intended analysis. |
| `PROCAR` | Site/orbital projections. | Appropriate `LORBIT` settings; can be large. |
| `LOCPOT` | A selected local-potential representation. | `LVTOT`/`LVHAR`; distinguish the full local potential from the electrostatic contribution. |
| `AECCAR*` | All-electron/related charge-density outputs. | `LAECHG`, for workflows that explicitly need them. |
| `XDATCAR` | Saved sequence of positions in ionic/MD runs. | Ionic run settings and output cadence; not a replacement for all restart state. |

Choose output flags according to the next task. For example, the starter `INCAR` sets `LWAVE` and `LCHARG` false, so those conventional restart files are not expected. The conventional LOBSTER workflow needs a compatible `WAVECAR`; enabling `LORBIT` and getting `PROCAR` is not a substitute. Enabling more output does not by itself make a DOS, band structure, or charge analysis converged.

### Check correctness before timing

```bash
tail -n 20 OSZICAR
tail -n 40 OUTCAR
grep -E 'TITEL|ENCUT|NBANDS|NKPTS' OUTCAR
```

An electronic-iteration line is not an ionic step. `OSZICAR` summarizes changes and residuals during the electronic loop; reaching `NELM` repeatedly calls for inspection. `F` and `E0` represent different energy conventions associated with the occupation scheme. When comparing runs, use the same intended energy definition rather than mixing the free energy, energy without entropy, and extrapolated energy. [OSZICAR reference](https://vasp.at/wiki/OSZICAR)

Distinguish three separate questions: did Slurm finish the job without an error, did the electronic loop converge, and—if this was a relaxation—did the ionic criterion converge before its step limit? Inspect forces/stress and the final geometry as well. A completed job, a final timing block, or a `CONTCAR` alone does not answer all three. The static setup run establishes neither a relaxed structure nor converged production parameters.

### Use the performance indicators

```bash
grep -E 'LOOP:|LOOP\+:' OUTCAR | tail -n 12
grep -A8 'total amount of memory used' OUTCAR
tail -n 30 OUTCAR
sacct -j JOB_ID --format=JobID,State,ExitCode,Elapsed,AllocCPUS,TotalCPU,MaxRSS
```

Replace `JOB_ID` with the submitted number and inspect the calculation step as well as the job-level row; accounting availability depends on cluster configuration.

| Indicator | How to use it |
|:----------|:--------------|
| `LOOP` real time | Compare representative electronic iterations when benchmarking MPI layouts. Initialization steps can differ from later steps. |
| `LOOP+` real time | Includes the electronic cycle and additional work such as forces. Useful for the cost of a complete ionic step. |
| End-of-`OUTCAR` elapsed time | Overall VASP run time, including setup and output overhead. Compare with scheduler elapsed time to understand other job overhead. |
| Actual `NBANDS` and iteration count | Check whether changing parallelization rounded up the bands or changed the amount of iterative work; then timings are not a pure scaling comparison. |
| `OUTCAR` memory breakdown | Identify large allocations reported by the printing rank. Do not interpret it as a measurement of every rank's peak memory. |
| Slurm `MaxRSS` | Peak reported resident memory for a task in a job step, not the sum across all nodes. Use step-level accounting and site tools to assess memory pressure. |
| `TotalCPU / (Elapsed × AllocCPUS)` | A rough utilization measure when all times use the same units and the same step is compared. It is not a scientific-convergence or parallel-efficiency score. |

For a useful benchmark, keep the structure, electronic settings, and starting charge/orbitals identical. Change one resource/layout choice at a time; record wall time, actual bands, iteration count, and memory. Compare allocated core-hours as well as speed: doubling the allocation for a small speedup may cost more overall. [VASP timing/tuning guidance](https://vasp.at/wiki/Optimizing_the_parallelization) and [Slurm accounting](https://slurm.schedmd.com/sacct.html)

For multiple directories, `cat-tools vasp check --prefix PREFIX` summarizes convergence runs. The independent [Slurm scanner](../../../../1_Tools/General/README.md) provides diagnostic and performance views; it complements inspection of VASP's own convergence output.

You are ready to begin convergence testing once the chosen structure/potentials are consistent, the intended `INCAR` settings appear in `OUTCAR`, the first electronic calculation finishes without unresolved errors, and its resource use is understood. Use the tested settings in the next workflow.

---

## References

Use these for detailed tag semantics and version-specific behavior; the practicalsequence above is the entry point.

- **Access:** [VASP user registration](https://vasp.at/info/faq/vasp_access/) and [licensing FAQs](https://vasp.at/info/faq/).
- **Structures:** [VESTA export](https://jp-minerals.org/vesta/en/doc/VESTAch4.html), [Pymatgen VASP I/O](https://pymatgen.org/pymatgen.io.vasp.html), and [ASE VASP I/O](https://docs.ase-lib.org/_modules/ase/io/vasp.html).
- **`POTCAR`:** [VASP Portal downloads](https://www.vasp.at/vasp-portal/), [available datasets](https://vasp.at/wiki/Available_pseudopotentials), [Preparing a POTCAR](https://vasp.at/wiki/Prepare_a_POTCAR), [choosing potentials](https://vasp.at/wiki/Choosing_pseudopotentials), and [dataset contents](https://vasp.at/wiki/POTCAR).
- **`INCAR` basics:** [INCAR](https://vasp.at/wiki/INCAR), [ENCUT](https://vasp.at/wiki/ENCUT), [PREC](https://vasp.at/wiki/PREC), [EDIFF](https://vasp.at/wiki/EDIFF), [ALGO](https://vasp.at/wiki/ALGO), [NELM](https://vasp.at/wiki/NELM), and [smearing](https://vasp.at/wiki/Smearing_technique).
- **Model and numerical refinements:** [ISPIN](https://vasp.at/wiki/ISPIN), [MAGMOM](https://vasp.at/wiki/MAGMOM), [LREAL](https://vasp.at/wiki/LREAL), [ROPT](https://vasp.at/wiki/ROPT), [LASPH](https://vasp.at/wiki/LASPH), and [ADDGRID](https://vasp.at/wiki/ADDGRID).
- **Ions and restarts:** [IBRION](https://vasp.at/wiki/IBRION), [ISIF](https://vasp.at/wiki/ISIF), [EDIFFG](https://vasp.at/wiki/EDIFFG), [ISTART](https://vasp.at/wiki/ISTART), and [ICHARG](https://vasp.at/wiki/ICHARG).
- **Parallelization:** [NCORE](https://vasp.at/wiki/NCORE), [NPAR](https://vasp.at/wiki/NPAR), [KPAR](https://vasp.at/wiki/KPAR), [NBANDS](https://vasp.at/wiki/NBANDS), [tuning](https://vasp.at/wiki/Optimizing_the_parallelization), and [memory](https://vasp.at/wiki/Not_enough_memory).
- **Outputs:** [OUTCAR](https://vasp.at/wiki/OUTCAR), [OSZICAR](https://vasp.at/wiki/OSZICAR), [NWRITE](https://vasp.at/wiki/NWRITE), [LORBIT](https://vasp.at/wiki/LORBIT), [LOCPOT](https://vasp.at/wiki/LOCPOT), [LAECHG](https://vasp.at/wiki/LAECHG), and [restart/output choices](https://vasp.at/wiki/Restart_and_output_files_cheat_sheet).
- **Scheduler:** [Slurm sbatch](https://slurm.schedmd.com/sbatch.html) and [Slurm sacct](https://slurm.schedmd.com/sacct.html).
