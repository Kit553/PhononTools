# FreqBallTools

For use without CAT, the [standalone distribution and workflow](../../2_Workflows/10_LatticeDynamics/01_Harmonic/02_DirectionalDOS/STANDALONE.md)
contain the scripts and all local helper source files.

FreqBallTools owns directional phonon spectra, atom-local structural directions,
frequency-ball geometry, and directional mode fatbands. It replaces the separate
calculators, analyzers, CSV splitters, site finders, and ball viewers in the original
folder.

## The canonical calculation

The raw result is exactly one `freqball.h5`. For each selected primitive atom it stores
six real spectral-tensor components on one frequency grid:

\[
D_{i,\alpha\beta}(f) = \frac{1}{N_q}\sum_{q\nu}
\operatorname{Re}\left[e_{i\alpha}(q\nu)e^*_{i\beta}(q\nu)\right]
g_\sigma(f-\omega_{q\nu}).
\]

The directional PDOS is the inexpensive contraction

\[
P_i(\hat n,f)=\hat n^\mathsf{T}D_i(f)\hat n.
\]

Consequently, 1000 sampled directions no longer invoke Phonopy 1000 times and no
direction-dependent spectra or cutoff-dependent descriptors are duplicated in HDF5.
The six stored components are `xx`, `yy`, `zz`, `yz`, `xz`, and `xy`; off-diagonal
components acquire the factor two during contraction.

The same HDF5 also embeds the complete primitive structure once, while the `atoms`
group identifies the tensor-selected sites. This makes later PhonoVis composition
portable without requiring another POSCAR merely to place a ball in structural
context; it is provenance and geometry, not a duplicate spectral representation.

The tensor calculation deliberately uses the full mesh without reciprocal-space
symmetry reduction. A direction in the laboratory Cartesian frame is not invariant
under the omitted rotations, so an irreducible mesh is not sufficient. Q-points are
processed in chunks; memory holds one eigenvector chunk, its Gaussian frequency
weights, and the final tensor, rather than the full mesh of eigenvectors. Optional MPI
splits q-points between ranks and reduces the tensor on rank zero; parallel HDF5 is not
required.

## Command interface

The operational procedure, with input preparation, execution order, checks, plot
commands and local/cluster mass scans, lives in the separate
[FreqBall workflow](../../2_Workflows/10_LatticeDynamics/01_Harmonic/02_DirectionalDOS/README.md).

The complete example is [`examples/freqball.yaml`](../examples/freqball.yaml).

```bash
uv sync
cat-freqball examples/freqball.yaml
```

Use `--action compute` to write only the tensor or `--action analyze` to derive new
views from an existing tensor without repeating the force-constant calculation.
`--format csv|parquet` selects one output format for the run.

For MPI calculation:

```bash
uv sync --extra mpi
mpirun -n 8 cat-freqball examples/freqball.yaml --action compute
cat-freqball examples/freqball.yaml --action analyze
```

## Structural direction families

The optional family preparation is a two-command path:

```bash
uv sync --extra families
cat-freqball-bvpa BVPA_output -o conduction_path.yaml
cat-freqball-families POSCAR \
  --conduction-path conduction_path.yaml \
  --mobile Na --tetra-center P \
  -o direction_families.yaml
```

When BVPA itself still needs to be run, StructureTools owns the WSL/Windows launcher
and can create the same handoff YAML in one command:

```bash
cat-prepare bvpa --cif structure.cif --cube bv.cube --ion Na1+ \
  -o BVPA --conduction-path conduction_path.yaml
```

This boundary is intentional: StructureTools launches the external structure/grid
calculation, while FreqBallTools owns the scientific path parser and BVS family model.

The first command parses explicit BVPA path connectivity into one canonical YAML and,
when present, checks `BVPA_gui.csv` and `BVPA_events.txt`. The second periodically
matches minima to mobile-ion structure sites and constructs:

- `BVS`: mobile site to the explicit path saddle;
- `NaNa`: mobile site to the connected mobile site image;
- `TetCenter`: mobile site to the nearest requested tetrahedral centers.

Parallel and antiparallel vectors are one axis and are deduplicated with retained
source names. The builder assigns crystallographic equivalence groups with Pymatgen.
Use the same primitive cell as the Phonopy calculation so site matching is unambiguous.

## Outputs

| File | Grain | Purpose |
|---|---|---|
| `freqball.h5` | atom × 6 tensor components × frequency | Only raw FreqBall result |
| `freqball_directions` | selected atom × sampled/exact direction | All moments, alignments, and direction scores |
| `freqball_family_spectra` | equivalence group × structural family × frequency | Exact BVS, NaNa, TetCenter, and isotropic spectra |
| `freqball_*_modes` | q-point × band | Mesh, path, or arbitrary-q fatband data |
| `manifest.yaml` | one run | Keys, units, columns, settings, and source identity |

These are distinct concepts, not competing layouts of the same data. There is no
wide-spectrum copy, column guide, extrema file, ZIP of per-direction tables, or
pre-rendered plot data. Select extrema by sorting/filtering `freqball_directions`; plot
all other views from the same table.

Directional descriptors currently include total spectral weight, average frequency,
central and raw second moments, standard deviation, a configurable low-frequency
weight/fraction/average, and inverse-square-weighted measures. Those settings are
recorded in the manifest and can be changed with `--action analyze` without changing
the tensor.

## Fatbands and Phono3py intersection

Every optional mode artifact starts with the standard PhonoPyTools mode schema and
adds:

- `mobile_participation`;
- absolute conduction and tetrahedral direction weights;
- the same weights divided by mobile participation;
- `direction_score`, the conduction fraction minus tetrahedral fraction;
- a weight column for every named subtype such as `BVS`, `NaNa`, or `TetCenter`.

For an atom and family, its directional weight is the maximum
`|direction · eigenvector|²` over equivalent reference axes. Contributions are then
summed over the selected atoms. This avoids multiplying a mode merely because the
family YAML lists symmetry-equivalent copies.

The q-vector and zero-based band columns are deliberately identical to the Phono3py
mesh identity. Phono3PyTools can therefore join `conduction_direction_weight` or any
other FreqBall mode column directly, without importing FreqBall internals.

## Artificial-mass scans

Mass scans reuse the ordinary FreqBall workflow. The setup command anchors external
input paths to the base YAML, changes only the requested element mass, and writes one
independent configuration and executable job per point:

```bash
cat-prepare mass-scan examples/freqball.yaml -o mass-scan \
  --species S --range 16 96 --points 17 --job-template job.sh \
  --submit-command sbatch
```

The optional job template must contain `{{FREQBALL_COMMAND}}`; resource directives,
modules, and the launcher stay entirely under user control. Other available tokens are
`CONFIG`, `WORKDIR`, `POINT_DIR`, `OUTPUT`, `SPECIES`, `MASS`, `INDEX`, and `JOB_NAME`.
Without a template, each generated `run.sh` is a directly executable local job. When
`--submit-command` is present, setup also writes (but never executes) `submit_all.sh`.
Setup itself never submits or runs a calculation.

After the jobs finish, one command derives the cross-mass tables:

```bash
cat-freqball-mass-scan collect mass-scan -o freqball-mass-analysis --format parquet
```

`freqball_mass_directions` retains the complete directional descriptors;
`freqball_mass_family_spectra` contains exact BVS/NaNa/TetCenter spectra;
`freqball_mass_spectra` contains the compact mobile-total and configurable
conduction/tetrahedral cone selections; and `freqball_mass_summary` contains their
moments versus assigned mass. These are four different grains, not reordered copies.
Incomplete scans fail by default; `--allow-missing` makes an intentional partial
collection.

The CATPlot recipes `mass-summary`, `mass-spectra`, and `mass-family-spectra` consume
those tables. Comparable 3D balls can be exported in a batch without a second renderer:

```bash
cat-freqball-mass-scan render mass-scan -o mass-balls \
  --metric average_frequency_THz --range 1 8 --equivalence-group 0
```

That command calls PhonoVis's common PyVista backend. An explicit shared `--range` is
recommended when comparing masses so color and radial deformation retain the same
meaning in every image.

## Visualization and CATPlot

FreqBallTools derives the scientific surface directly from `freqball.h5`, then hands
that surface and the structural-family arrows to PhonoVis's shared PyVista backend:

```bash
uv sync --extra cloud
cat-freqball-view freqball.h5 --families direction_families.yaml \
  --metric average_frequency_THz
```

The `cat-freqball-view` domain command remains the convenient entry point. The shared
renderer supports interactive rotation, atom or equivalence-group selection, shared scalar
ranges, structural-family arrows, publishable screenshots, rotating GIFs, and VTK mesh
export. Useful switches include `--atom`, `--equivalence-group`, `--range`,
`--deformation`, `--screenshot`, `--gif`, and `--export-mesh`.
Interactive rendering uses the local display. Headless screenshot/GIF jobs require a
working EGL, OSMesa, or Xvfb backend for VTK; ordinary desktop PyVista needs no extra
setup.

The 2D plotting layer remains CATPlot. Ready-to-edit recipes are in
[`recipes/freqball`](../recipes/freqball):

- `family-spectra.plot.yaml` for exact structural-direction spectra;
- `direction-map.plot.yaml` for descriptor/alignment relations;
- `fatband.plot.yaml` for a dispersion whose point area is mobile participation and
  whose color is the conduction–tetrahedral direction score.
- `mass-summary.plot.yaml`, `mass-spectra.plot.yaml`, and
  `mass-family-spectra.plot.yaml` for artificial-mass results.

FreqBallTools contains no PyVista calls. This keeps camera, background, scalar-bar,
screenshot, and GIF behavior consistent with the Phono3Py cloud and other PhonoVis
views. The viewer and CATPlot consume canonical data; neither writes another
scientific data representation.

## v6 compatibility and migration contract

`cat-freqball` accepts `FreqBallz_v6.py` configuration files directly. Both dialects
use paths relative to the configuration file; v6 itself used the process working
directory for force inputs/output, so those paths may need correction. To write an
editable canonical configuration with anchored external input paths:

```bash
cat-freqball old.yaml --write-config freqball.yaml
cat-freqball freqball.yaml --mass S 48 --output sulfur48.h5
```

`--mass SYMBOL U` is repeatable for distinct elements, overrides configured element
masses, and cannot be mixed with a primitive mass list. All masses must be finite
and positive. `--output` changes the tensor path only. `--write-config` writes a
configuration and exits without computing; output paths then belong to the new
configuration's directory. Mass-scan setup also accepts v6 YAML.

| v6 setting | Canonical setting |
| --- | --- |
| `phonopy.phonopy_yaml`, `force_sets`, `born`, `masses` | Same keys under `system` |
| `phonopy.mesh` | `spectrum.mesh` |
| `ions.indices` | `atoms.indices` (zero-based primitive atoms) |
| `sampling.freq_min`, `freq_max` | `spectrum.frequency_range_THz` |
| `sampling.freq_pitch`, `sigma` | `spectrum.pitch_THz`, `sigma_THz` |
| `sampling.n_directions` | `analysis.directions.count` |
| `selection.low_freq_cutoff` | `analysis.directions.low_frequency_cutoff_THz` |
| `direction_families_file` | `families.file` |
| `direction_family_match_tolerance_A` | `families.match_tolerance_A` |
| `output.file`, `spectra_dtype` | `output.tensor`, `spectra_dtype` |

Normalized v6 configurations use `gamma_center: false` to preserve the original
Phonopy default on even meshes; the CAT example explicitly uses `true`. Frequency
grid endpoints follow Phonopy's DOS grid convention. Gaussian spectra and the seven
original moments retain v6 definitions, including the nonnegative low-frequency
window. Sampling accepts two or more directions, although triangulated 3D balls
need a nondegenerate sphere with more points. Floating-point storage may be float32
or float64 (default float64 for CAT, float32 for converted v6); accumulation is float64.

The direction table retains v6's best conduction/tetrahedral reference indices,
angles in degrees, and maximum alignment strength, alongside the signed score.
`avg_freq`, `second_central_moment`, `std_freq`, `low_freq_weight`,
`low_freq_fraction`, and `low_freq_avg` become `average_frequency_THz`,
`second_central_moment_THz2`, `frequency_std_THz`, `low_frequency_weight`,
`low_frequency_fraction`, and `low_frequency_average_THz`. `total_weight` is unchanged.
Reference indices index the atom's references within that category, in YAML order.
Absent families produce NaN alignments/angles and reference index -1.

The **output format changes**: v6's `metadata/schema_version=2` atom-wise PDOS file
is not the tensor schema identified by the root `freqball_schema` attribute. Recompute
from the force inputs; do not rename an old file and treat it as a tensor. Raw spectra
for arbitrary directions remain available through `directional_spectra`. Full
effective configuration and input hashes are recorded for provenance. HDF5 publication
is atomic, so readers see a completed result rather than a partially written tensor.

## Parallel execution and reproducibility

MPI is optional. A serial call and a multi-rank call use the same scientific path:

| Stage | Work partition | Final output |
| --- | --- | --- |
| Spectral tensor | Strided q-points, bounded `spectrum.qpoint_chunk` | Tensor sum reduced to rank zero |
| Direction descriptors | Selected atoms | Ordered table gathered to rank zero |
| Mesh/arbitrary-q fatbands | `modes.qpoint_chunk` chunks (default 64) | Projected tables gathered, not full eigenvectors |
| Band fatbands | Complete path segments | Ordered segments with cumulative distances |
| Mass calculations | Independent generated jobs; local `run --workers N` or scheduler | One tensor and analysis directory per mass |
| Mass collection | `collect --workers N` processes mass points | Four combined tables and manifest |

Band segments preserve Phonopy's branch connection and NAC handling at Gamma.
Family-spectrum aggregation and final I/O use rank zero; 3D batch rendering is
serial. Avoid nested oversubscription: total processes are concurrent mass jobs
multiplied by ranks per job, with BLAS/OpenMP threads per rank controlled separately.
MPI stage errors are propagated collectively before continuing to another stage.

Local scan execution writes per-point `run.log` and a scan-level exit-code table;
failed recorded jobs block collection. Setup remains generation-only. Collection
checks the job table's assigned masses against HDF5, unchanged masses of other
species, force-input hashes when available, structure, mesh, grid and broadening.
`--allow-missing` permits absent outputs, not failed recorded jobs or inconsistent
results. No automatic resume/skip policy is inferred from an existing HDF5 file.

The default directional mode mesh now uses **no symmetry or time-reversal reduction**.
Explicitly requesting reduction gives representative q-points with their weights;
the local directional weights are not rotated onto omitted q-points. For q-tables,
set `kind: qpoints` and either `qpoints: [[...], ...]` or `from_table` plus optional
`where` filters. Set `nac_direction` when an explicit Gamma limit is required.

`--action analyze` recomputes configured mode outputs from force data. It verifies
selected atoms, primitive structure/masses and available input hashes against the
tensor. Set `modes: []` to reanalyze a portable tensor without the original force data.

## Interpretation of the conduction fatband

For mode eigenvectors in Phonopy's normalization, mobile participation is
`sum_selected_atoms,sum_xyz |e|²`. A family weight uses each atom's maximum squared
projection over its reference axes, then sums over atoms. Dividing by mobile
participation gives a conditional directional fraction; zero mobile participation
has an undefined fraction (NaN). Opposite or repeated axes do not multiply weights.
The conduction and tetrahedral families can overlap, so their fractions need not
sum to one. The signed difference ranges from -1 to +1 where defined.

The recipe draws the dispersion in gray, scales scatter area by mobile participation
within the plotted data range, and uses the signed difference for color. This is a
harmonic mode-character diagnostic related to the chosen conduction geometry, not
an ionic conductivity or migration-rate calculation. Mass scans rediagonalize the
mass-modified dynamical matrix; they do not simply rescale all frequencies by one
square-root factor.

## Verification

From the installed repository environment:

```bash
python -m unittest discover -s tests -p 'test_freqball*.py'
python -m unittest discover -s tests -p test_recipes.py
```

Regression coverage includes the v6 mesh/projection call and seven moments after a
mass change, oblique primitive coordinates, exact local references, float32 storage,
NAC-aware band segments against native Phonopy, serial versus concurrent collective
orchestration (including idle ranks and rank-local errors), and local parallel mass
jobs followed by serial/parallel collection. The collective test simulates MPI;
a deployment still needs an actual `mpirun` smoke test with its MPI/mpi4py stack.
