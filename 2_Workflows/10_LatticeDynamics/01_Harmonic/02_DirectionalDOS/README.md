# FreqBall: directional spectra, conduction fatbands and mass scans

**Standalone scripts are stored in this directory.** Start with
[the standalone workflow](STANDALONE.md) to run `FreqBallz.py`, the family/mass/view
helpers and `FreqPlot.py` without installing CAT or CATPlot. Copy the complete
directory, including `freqball_helpers/`. The procedure below uses CAT entry points.

This is the **execution procedure**, from completed harmonic force data to plots and
cross-mass comparisons. The [FreqBall reference](../../../../3_Analysis/docs/freqball.md)
defines the algorithms, columns and compatibility rules. Complete the
[Phonopy force workflow](../01_PhonoPy/README.md) or [DFPT](../03_DFPT/README.md) first.
No new DFT calculations are needed for this procedure or its artificial-mass scans.

## 1. Prepare a working directory

Install CAT from the repository root (`uv sync`; add `--extra mpi` for MPI and
`--extra families --extra cloud` for family construction and 3D viewing).
Make the installed `cat-*` commands available in your shell. Install CATPlot following
[its instructions](../../../../4_Plotting/README.md) if it is not already available.

Copy this directory's `freqball.yaml` into a new analysis directory alongside your
primitive structure (`POSCAR.primitive`). Edit `system.path` to the completed force
calculation. The input must contain compatible Phonopy metadata and either FORCE_SETS,
FORCE_CONSTANTS, or force constants embedded in the YAML. For FORCE_CONSTANTS replace
`force_sets` with `force_constants: FORCE_CONSTANTS`; for embedded data use
`force_sets: null`. Add `born: BORN` only when that file exists and NAC is intended.

The structure used for families must have the **same lattice and sites as Phonopy's
primitive cell**. Atom indices are zero-based primitive indices, not supercell indices.
Use `atoms.elements: [Na]` or an explicit `atoms.indices` list, never both.

For an existing v6 configuration, normalize it first:

```bash
cat-freqball /path/to/v6.yaml --write-config freqball.yaml
```

This anchors input paths to the old YAML's directory. Inspect the generated output
paths, add the `modes` section from the example below, and use the rest of this
procedure. Direct execution of the old YAML also works. If the old run relied on a
different shell working directory, correct its relative paths before conversion.
Old v6 atom-wise HDF5 files must be recomputed as tensors.

## 2. Prepare each ion's conduction and tetrahedral directions

If `direction_families.yaml` from v6 already matches the primitive cell, reuse it.
Otherwise, obtain BVPA output using the [SoftBV procedure](../../../99_Misc/README.md)
and run:

```bash
cat-freqball-bvpa /path/to/BVPA_output -o conduction_path.yaml
cat-freqball-families POSCAR.primitive \
  --conduction-path conduction_path.yaml --mobile Na --tetra-center P \
  -o direction_families.yaml
```

Check the reported site matching. The file contains site-to-saddle `BVS`, connected
mobile-site `NaNa`, and mobile-site-to-tetrahedral-center `TetCenter` axes. Change
`Na` and `P` to the species in your material. Do not substitute a conventional cell
or a supercell with different coordinates. Without families, frequency balls remain
available but conduction/tetrahedral classification is undefined.

## 3. Set the numerical and analysis parameters

Edit the copied `freqball.yaml`:

- Set the mesh, frequency window, Gaussian width and pitch. Include negative
  frequencies when inspecting instabilities. A window that clips shifted modes
  will bias moments and mass comparisons.
- Set `analysis.directions.count` and the low-frequency cutoff. The sphere and each
  atom's exact local reference axes are both included.
- Set explicit reciprocal-fractional `modes.paths` and matching labels for your
  crystal. The supplied path is only an example. Omitting `paths` requests Phonopy's
  automatic path (and requires Seekpath in the environment).
- For mode meshes, retain `mesh_symmetry: false` and `time_reversal: false` when
  inspecting local directional character across the full Brillouin zone.

## 4. Compute and analyze

For a small serial check:

```bash
cat-freqball freqball.yaml
```

For production MPI work, run inside your allocation:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mpirun -n 8 cat-freqball freqball.yaml
```

The full run writes `freqball.h5`, then the tables and `manifest.yaml` under
`analysis/`. Rank zero writes files; MPI-enabled h5py is unnecessary. Tensor q-points,
mode-mesh chunks, band segments and per-atom directional analysis are distributed.
Band segments stay intact for branch tracking and NAC at Gamma. More ranks than
segments cannot accelerate a short band path further. Family-spectrum aggregation
and final table writing occur on rank zero.

To separate the stages:

```bash
mpirun -n 8 cat-freqball freqball.yaml --action compute
mpirun -n 8 cat-freqball freqball.yaml --action analyze
```

`analyze` reuses the tensor for directions and spectra. Configured fatbands still
require the original force inputs and an eigenvector calculation. Set `modes: []`
when changing only the ball sampling or descriptor cutoffs. Changing masses, mesh,
frequency window, broadening or NAC requires a new `compute`.

Before increasing resources, check that selected sites and active masses are correct
in the HDF5 and that `analysis/freqball_directions.csv` contains both sampled and exact
reference rows. Converge the mesh and broadening using the moments and family spectra.

## 5. Inspect the balls and plot the conduction fatband

```bash
cat-freqball-view freqball.h5 --families direction_families.yaml \
  --metric average_frequency_THz --atom 0
catplot render fatband.plot.yaml
```

Copy this directory's `fatband.plot.yaml` into the working directory first. It reads
`analysis/freqball_band_modes.csv` and writes `plots/draft/freqball-fatband.pdf` (CATPlot's default draft stage). Gray lines
show the dispersion; point size is scaled by mobile-ion participation, and color is
the conduction-minus-tetrahedral fraction on a fixed [-1, 1] scale. Edit labels/ticks
for your chosen path. The weights measure harmonic vibrational character along the
structural axes; they do not calculate ionic conductivity or hopping rates.

For spectra, copy the
[family-spectra recipe](../../../../3_Analysis/recipes/freqball/family-spectra.plot.yaml),
set its data path to `analysis/freqball_family_spectra.csv` and output directory to
`plots`, then render it with `catplot render`.

## 6. Generate an artificial-mass scan

Keep geometry, force constants, Born data, q-mesh and spectral settings fixed. Select
the species whose mass changes (it need not be the mobile species):

```bash
cat-prepare mass-scan freqball.yaml -o mass-scan \
  --species S --range 16 96 --points 17
```

Use `--values 16 32 64` instead for an explicit grid. Setup writes a config and executable
`run.sh` per mass, plus an authoritative job table and manifest. It neither runs jobs
nor submits them. Each point inherits the base workflow, including fatbands unless
`modes: []` was selected. Other existing element mass overrides remain in force.

## 7. Execute scan points

For bounded local parallel execution:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cat-freqball-mass-scan run mass-scan --workers 4
```

Each point has a `run.log`; `freqball_mass_status.csv` records exit codes. A failed job
makes the command fail, and collection rejects recorded failures. Inspect the logs,
fix the inputs and rerun. The runner executes every listed job; it does not silently
skip existing results. Do not run this orchestrator under MPI.

For a cluster, regenerate the scan in a **different output directory** with a template:

```bash
cat-prepare mass-scan freqball.yaml -o mass-scan-cluster \
  --species S --range 16 96 --points 17 \
  --job-template mass-job.sh --submit-command sbatch
bash mass-scan-cluster/submit_all.sh
```

Copy and edit `mass-job.sh` from this directory first (partition, memory, walltime,
modules/environment and MPI launcher). The generated script contains the exact
`cat-freqball` command in place of `{{FREQBALL_COMMAND}}`. Cluster submission parallelizes
mass points; each allocation can also use MPI within its point. Budget total ranks
and memory accordingly. Use the local runner only for scripts intended for local
execution, not scripts relying on scheduler allocations.

## 8. Collect, compare and render

After all points finish:

```bash
cat-freqball-mass-scan collect mass-scan -o mass-analysis --workers 4
cat-freqball-mass-scan render mass-scan -o mass-balls \
  --metric average_frequency_THz --equivalence-group 0 --range 1 8
```

Use your actual equivalence group and a common range suitable for all masses. Headless
3D rendering needs a working VTK display backend. The collection step verifies masses,
structure, mesh, frequency grid, broadening and recorded force-input identity. Missing
points fail by default; `--allow-missing` permits deliberately incomplete collections.

The outputs are `freqball_mass_directions`, `freqball_mass_family_spectra`,
`freqball_mass_spectra` (total and cone-selected spectra), and `freqball_mass_summary`.
Use the [mass plot recipes](../../../../3_Analysis/recipes/freqball/), editing their
input paths to `mass-analysis/` and output directories to `plots/` before rendering.
Use the same frequency window, cutoff and cone thresholds for every comparison.
