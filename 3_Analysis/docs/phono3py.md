# Phono3PyTools

The Phono3py collection separates faithful data extraction from optional reduction and
plotting. The normal workflow has no plotting state. CATPlot receives tables that are
already shaped for the requested visual, while `mesh_modes` remains the lossless
analysis table from which new views can be derived.

## Normal workflow

```bash
cat-phono3py path/to/workflow.yaml
```

Paths are resolved relative to the workflow file. Automatic kappa discovery selects
complete `kappa-m*.hdf5` files and excludes per-grid-point `-gN` files. Explicit paths
or globs always take precedence. Only requested temperatures are sliced from HDF5.

| Artifact | Row identity | Purpose |
|---|---|---|
| `mesh_modes` | dataset, temperature, irreducible q-index, band | Canonical modal properties and projection joins |
| `thermal_conductivity` | dataset, temperature, channel | Total, population, coherence, RTA, or LBTE tensors when present |
| named profile | groups, frequency bin | Weighted mean, median, and interquartile range |
| named density | groups, frequency bin, value bin | Complete weighted 2D histogram for heatmaps |
| named spectrum | groups, frequency bin | Signed frequency sum, density, cumulative sum, and fraction |
| named accumulation | groups, coordinate bin | Direct cumulative modal property versus frequency, MFP, or another coordinate |
| `jdos_mesh` | temperature, channel, frequency | Weighted mesh JDOS with separate difference and sum channels |
| `gruneisen_modes` | mesh q-point/band or path segment/point/band | Mesh or true-path Grüneisen data |
| `grid_path_modes` | segment, point, band | Existing mesh modes selected on an exact commensurate path |
| `mode_pairs` | temperature, q-point, band i, band j | Pair-resolved Wigner coherence |
| `interactions` | dataset, triplet, band i, band j, band k | Explicitly selected full phonon-phonon interaction tensor |

The modal lifetime is derived from the Phono3py linewidth as

\[
\tau = \frac{1}{4\pi\Gamma}.
\]

Modal conductivity from a complete mesh file is divided by the sum of irreducible-grid
weights. A per-grid-point file is divided by the full mesh size. Both normalizations
are recorded in the artifact manifest.

For generalized regular grids, the full size is the absolute determinant of the grid
matrix. The exporter does not invent simple `q × mesh` addresses for such a grid;
q-vectors stay authoritative and the grid matrix is retained in the manifest.

## Mesh representations

Use the smallest transformation matching the scientific question:

- `profile` for a frequency-resolved central trend and spread, such as mean Γ at each
  temperature;
- `density` for a dense two-variable cloud where overplotting would hide structure;
- `spectrum` for additive quantities such as cumulative thermal conductivity;
- raw `mesh_modes` scatter for mode character, Na participation, or another meaningful
  color channel.

Profiles and densities use irreducible-grid weights by default. Empty density bins are
kept so heatmap geometry is complete. A spectrum never clips negative contributions.
Bins can be an integer count, explicit edge list, or a mapping with `start`, `stop`,
`width` or `count`, and `scale`.

The bundled recipes cover κ(T) by boundary MFP, cumulative κ by frequency, Γ profiles
by temperature, Γ density, Γ and Grüneisen colored by Na participation, Grüneisen
colored by an external conduction-direction weight, averaged PP, and both JDOS channels.

For direct accumulation from `mesh_modes`, use an `accumulation` view. It does not
shell out to `phono3py-kaccum` or parse its human-oriented stdout:

```yaml
mesh:
  views:
    - {name: kappa_accum_frequency, kind: accumulation,
       value: mode_kappa_iso_W_mK, coordinate: frequency_THz, bins: 160}
    - {name: kappa_accum_mfp, kind: accumulation,
       value: mode_kappa_iso_W_mK, coordinate: mean_free_path_A,
       bins: {start: 0.1, stop: 1.0e6, count: 120, scale: log}}
```

The two outputs have the same tidy accumulation columns, but remain separate artifacts
because their coordinate units and physical questions differ. The modal conductivity
stored in `mesh_modes` is normalized once by the full grid size; accumulation then
applies the irreducible q-point multiplicity once. Signed modal contributions are
preserved. Set `weight: null` only for a property that is already multiplicity-weighted.

## Harmonic and FreqBall intersections

The `mesh.harmonic` section loads `phono3py.yaml` with `fc2.hdf5` through PhonoPyTools
and evaluates eigenvectors at the exact q-points stored in the anharmonic file. Element,
site, group, direction, and participation-ratio columns are joined by periodic q-vector
and zero-based band index.

`mesh.joins` accepts CSV or Parquet tables using the same `qx`, `qy`, `qz`, and `band`
identity. This is the FreqBall interface. Joins require complete coverage by default,
include `system` and temperature when available, reject duplicate external keys, and
understand periodically equivalent coordinates such as `-0.5` and `+0.5`. A deliberate
partial join requires `require_complete: false`.

## Mesh paths and true paths

`grid_path` selects points already present on the Phono3py mesh and fails clearly when
a requested point is not commensurate. An odd 21-point grid, for example, cannot
contain a reduced coordinate of exactly `1/2`. Use an even grid or a true path.

Grüneisen band YAML is a true arbitrary path and retains segment, point, trace, and
distance information. Other interaction properties remain mesh-derived unless they
were calculated separately at the requested grid addresses.

An optional Lorentzian view can be generated from the exact path modes and their
Phono3py linewidths:

```yaml
lorentzian:
  name: lorentzian_band
  frequency_range_THz: [0, 15]
  frequency_points: 800
  minimum_gamma_THz: 1.0e-4
```

This is deliberately labeled an approximate visualization. It sums independent
diagonal Lorentzians and therefore cannot describe mode mixing or the off-diagonal
spectral function. Every run emits a warning, explicitly counts adjacent modes whose
half-widths overlap, and records that count and the linewidth convention in the
manifest. The `lorentzian-band.plot.yaml` recipe renders the resulting regular grid.

## Experimental full-BZ cloud

```bash
cat-phono3py-cloud mesh_modes.csv phono3py.yaml \
  --value gamma_THz --dataset kappa-m212121 --temperature 300 \
  --frequency-range 0 8 --reduce max
```

The command constructs the actual first Brillouin zone as the reciprocal-space
Wigner-Seitz cell. For scalar quantities it expands irreducible q-points into symmetry
stars, folds them into that cell, and removes coincident copies. Bands at the same
q-point are reduced by `max` by default; `mean`, `sum`, and `none` are explicit options.
The visibility cutoff slider means:

- `0`: display every retained point;
- values between `0` and `1`: hide progressively weaker points;
- `1`: display only points at the maximum normalized strength.

`--scale log` applies the cutoff in log-strength space. `--absolute` ranks signed data
by magnitude while retaining the signed value for color. `--export cloud.csv` writes
the expanded points for inspection or another renderer.

Cloud preparation remains in Phono3PyTools, while its renderer-neutral layers are
displayed by PhonoVis's common PyVista backend. Symmetry expansion assumes the chosen
value is a scalar invariant. Disable it for a
direction-dependent projection unless its transformation law is handled before export:

```bash
cat-phono3py-cloud mesh_modes.csv phono3py.yaml \
  --value conduction_direction_weight --no-symmetry-expansion
```

The viewer is an optional dependency because VTK is large. Install it with
`pip install 'cat-tooling[cloud]'` or `uv sync --extra cloud`.
