# Standalone FreqBalls workflow

**Use:** Run FreqBall analysis and mass scans with the standalone scripts, instead of the full CAT installation, might be handy on restrictive environments

**Files:** [Script bundle](.) · [requirements.txt](requirements.txt) · [freqball.yaml](freqball.yaml) · [mass-job.sh](mass-job.sh).

Copy **this whole directory**, including `freqball_helpers/`, to your workstation or cluster. All Python source needed by these commands is stored here as ordinary files. No CAT installation, repository checkout, CATPlot, or repository-specific PYTHONPATH is required. The top-level scripts call the adjacent helper modules; copying a single launcher without its helper directory is insufficient.

## One command for the workflow

`python FreqBallz.py --help` exposes `run`, `compute`, `analyze`, `bvpa`, `families`, `mass-scan`, `view`, and `plot`. Replace `cat-freqball` in the [main workflow](README.md) with `python FreqBallz.py`. For example:

```bash
python FreqBallz.py run freqball.yaml
python FreqBallz.py plot fatband analysis/freqball_band_modes.csv -o fatband.pdf
python FreqBallz.py mass-scan setup freqball.yaml -o mass-scan --species S --values 24 32 48
```

Top-level YAML `projections` supplies ordinary element/site, Cartesian, lattice-axis and plane-normal spectra, descriptors and fatbands through this same workflow. Collection includes the corresponding named-projection mass tables. The [reference](../../../../3_Analysis/docs/freqball.md) defines normalization, coordinate conventions and output columns.

## Files and dependencies

| File | Purpose |
| --- | --- |
| `FreqBallz.py` | Unified calculation, family, mass-scan, ball and plotting commands |
| `BVPA_Parser.py` | Convert existing BVPA output into conduction-path YAML |
| `FreqSiteFinder.py` | Construct local conduction and tetrahedral reference axes |
| `FreqMassScan.py` | Set up, run, collect and render mass scans |
| `FreqView.py` | Interactive balls, screenshots, GIFs and mesh export |
| `FreqPlot.py` | Standalone fatband, spectral and mass-summary plots |
| `freqball_helpers/*.py` | Complete local calculation, I/O, family, mass-scan and rendering implementation |
| `freqball.yaml` | Editable input configuration |
| `mass-job.sh` | Editable SLURM job template |

Use Python 3.12 or newer. From this directory:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Install optional features only when needed:

```bash
python -m pip install pymatgen  # Construct families from a structure
python -m pip install pyvista  # View/export 3D balls
python -m pip install mpi4py   # Requires a compatible MPI runtime/launcher
python -m pip install seekpath # Automatic band paths; explicit paths need no Seekpath
python -m pip install pyarrow  # Parquet; CSV needs no extra package
```

## 1. Connect completed harmonic force data

Edit `freqball.yaml`: set `system.path`, the Phonopy input filename, and the mobile species or zero-based primitive atom indices. Use either `force_sets: FORCE_SETS`, `force_constants: FORCE_CONSTANTS` (remove `force_sets`), or `force_sets: null` for force constants embedded in the Phonopy YAML. Enable `born: BORN` only for compatible Born data. Set the mesh, spectral window and broadening, then replace the example band path/labels with ones appropriate to your crystal.

Paths are relative to the configuration file, except force-input filenames, which are relative to `system.path`. Preserve the same primitive cell and site ordering when creating the structural direction families.

Existing v6 configurations can be executed directly or converted:

```bash
python FreqBallz.py /path/to/old-v6.yaml --write-config freqball.yaml
```

Conversion anchors external inputs to the old YAML directory. Correct old paths that previously relied on a different working directory, and inspect the generated output paths. Old atom-wise v6 HDF5 files must be recomputed; they are not tensor files.

## 2. Prepare structural directions

Reuse an existing compatible `direction_families.yaml`, or parse existing BVPA results and build the families from the same primitive structure used by Phonopy:

```bash
python BVPA_Parser.py /path/to/BVPA_output -o conduction_path.yaml
python FreqSiteFinder.py POSCAR.primitive --conduction-path conduction_path.yaml \
  --mobile MOBILE_ELEMENT --tetra-center CENTER_ELEMENT -o direction_families.yaml
```

Replace `MOBILE_ELEMENT` and `CENTER_ELEMENT` with your mobile species and tetrahedral-center element symbols. `NaNa` is an output-family label, not a requirement to use sodium. This creates BVS site-to-saddle, NaNa connected-site, and TetCenter reference axes. Without a family file, remove the `families` section: ordinary balls remain available, but conduction/tetrahedral classification and its colored fatband are undefined. Running SoftBV/BVPA itself is an external preparation step.

## 3. Compute and analyze

```bash
python FreqBallz.py freqball.yaml
```

For MPI inside an appropriate allocation:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
mpirun -n 8 python FreqBallz.py freqball.yaml
```

The example writes `freqball.h5` and tables under `analysis/`. The tensor's q-points, mode-mesh chunks, band segments and per-atom directional analysis use the same MPI implementation as the packaged version. Rank zero writes output; MPI-enabled h5py is unnecessary. Band segments stay intact for branch tracking and NAC at Gamma.

To reuse a tensor after changing the sampled directions or descriptor cutoffs:

```bash
python FreqBallz.py freqball.yaml --action analyze
```

Set `modes: []` to avoid recalculating fatbands during reanalysis. Changing masses, mesh, broadening or frequency window requires a new computation. A one-off mass change is supported with `--mass ELEMENT MASS --output trial-mass.h5` (replace both tokens); use the scan helper below when multiple runs must retain separate analysis outputs.

## 4. Plot the results without CATPlot

```bash
python FreqPlot.py fatband analysis/freqball_band_modes.csv -o plots/fatband.pdf
python FreqPlot.py spectra analysis/freqball_family_spectra.csv \
  --equivalence-group 0 -o plots/family-spectra.pdf
python FreqView.py freqball.h5 --families direction_families.yaml \
  --metric average_frequency_THz --atom 0
```

The standalone fatband uses gray dispersion lines, point area proportional to mobile participation with a fixed scale, and color for the conduction-minus-tetrahedral fraction on [-1, 1]. Optional `--ticks` and `--labels` set band-path labels; use the actual path distances in the mode table. PNG and SVG also work via the output suffix. These weights describe harmonic vibrational character, not ionic conductivity or migration rates. The existing `.plot.yaml` recipe is optional for CATPlot users.

3D screenshots/GIFs require a working PyVista/VTK display backend on headless systems.

## 5. Set up and run mass scans

```bash
python FreqMassScan.py setup freqball.yaml -o mass-scan \
  --species ELEMENT --values MASS_1 MASS_2 MASS_3
python FreqMassScan.py run mass-scan --workers 4
python FreqMassScan.py collect mass-scan -o mass-analysis --workers 4
```

Setup only writes files. Each generated job invokes the **current Python interpreter and this bundle's absolute `FreqBallz.py` path**, with no `cat-*` commands. Keep the bundle and environment in place while jobs run, or regenerate after moving them. For an alternative cluster launcher/interpreter, pass `--command` explicitly. Replace `ELEMENT` with a species in your structure and the `MASS_*` tokens with positive masses in atomic mass units. An evenly spaced scan can instead use `--range MIN_MASS MAX_MASS --points COUNT`, with the bounds/count replaced. Other configured element masses are retained.

Each point stores its own config, HDF5, analysis and `run.log`. The status table records failures; collection rejects failed jobs, inconsistent force inputs/settings or wrong masses. The local runner executes every listed job and uses bounded concurrency. Do not launch the scan orchestrator under MPI; MPI belongs inside individual jobs.

For SLURM, first edit the resources/environment in `mass-job.sh`, then:

```bash
python FreqMassScan.py setup freqball.yaml -o mass-scan-cluster \
  --species ELEMENT --values MASS_1 MASS_2 MASS_3 \
  --job-template mass-job.sh --submit-command sbatch
bash mass-scan-cluster/submit_all.sh
```

The template's `{{FREQBALL_COMMAND}}` is replaced by the standalone command. Account for concurrent mass jobs × MPI ranks × numerical-library threads when sizing resources.

After collection:

```bash
python FreqPlot.py mass-summary mass-analysis/freqball_mass_summary.csv \
  --equivalence-group 0 -o plots/mass-summary.pdf
python FreqPlot.py spectra mass-analysis/freqball_mass_spectra.csv \
  --equivalence-group 0 -o plots/mass-spectra.pdf
python FreqMassScan.py render mass-scan -o mass-balls \
  --equivalence-group 0 --metric average_frequency_THz --range 1 8
```

Choose a common frequency window, cutoff, cone thresholds and rendering range for the whole scan. Use your actual equivalence group and suitable scalar limits.

## Maintenance within the repository

The helper files are checked-in copies of the tested implementation with local imports. They do not locate or import the repository at runtime. To update them after changing the shared implementation, run from the repository root:

```bash
python 3_Analysis/scripts/sync_freqball_standalone.py
python 3_Analysis/scripts/sync_freqball_standalone.py --check
```

The standalone regression test copies this folder outside the repository and executes its commands with CAT imports unavailable. Shared-source synchronization is checked separately so corrections do not silently diverge between distributions.
