# FreqBallTools

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
processed in chunks and only the small final tensor remains in memory. Optional MPI
splits q-points between ranks and reduces the tensor on rank zero; parallel HDF5 is not
required.

## Workflow

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
