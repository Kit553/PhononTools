# PhonoVis

PhonoVis owns 3D rendering, eigenvector inspection, and visual export. PhonoPyTools exports the band, DOS and participation tables consumed by CATPlot. PhonoVis calculates the same aligned mode columns in memory when an interactive dispersion needs colors, so no second set of near-identical CSV files is created.

## Commands

```bash
uv sync --extra visualization

# Click a dispersion point and inspect its animated eigenvector beside the DOS.
cat-phonovis explore examples/phonovis.yaml

# Inspect one exact arbitrary q-point and zero-based band in the rotatable viewer.
cat-phonovis explore examples/phonovis.yaml --q 0.5 0 0 --band 7

# Use the lighter Matplotlib viewer for the explicit mode.
cat-phonovis explore examples/phonovis.yaml --q 0.5 0 0 --band 7 --viewer matplotlib

# Export the q-point and band in the YAML, or override them on the command line.
cat-phonovis export examples/phonovis.yaml
cat-phonovis export examples/phonovis.yaml --q 0.5 0 0 --band 7 --amplitude 0.8
```

`explore.kind: dashboard` opens the clickable view. `explore.kind: mode` instead opens the q-point in the `mode` section, which is convenient for a saved inspection setup. The dashboard uses Matplotlib because it can place dispersion, projected DOS, and a rotatable 3D mode in one reliable window. The explicit viewer defaults to PyVista and has live phase and amplitude sliders.

Dashboard controls are intentionally small:

| Input | Action |
|---|---|
| click on dispersion | nearest sampled q-point and band |
| space | pause/resume phase animation |
| left/right | previous/next sampled q-point |
| up/down | previous/next zero-based band |
| `+` / `-` | increase/decrease amplitude |
| `L` | labels on/off |
| `A` | arrows on/off |

The dispersion `color_by` value is any numeric column from the ordinary PhonoPyTools mode schema or from a joined canonical mode table. Common local choices are `Na_participation`, `participation_ratio`, and `group_velocity_x`. A FreqBall path table can add `direction_score`; a commensurate Phono3Py grid path can add `gamma_THz` or `gruneisen`. Configure joins under `explore.dispersion`:

```yaml
color_by: direction_score
color_map: coolwarm
color_range: [-1, 1]
joins:
  - path: ../freqball-analysis/freqball_band_modes.csv
    columns: [direction_score]
    require_complete: true
```

Path joins validate segment, point index, periodic q-vector, zero-based band, and system. Mesh joins use periodic q-vector and band. Duplicate keys and incomplete coverage fail instead of silently generating a partly coloured dispersion.

## Shared 3D boundary

Scientific packages create small renderer-neutral layers: surfaces and family vectors from FreqBallTools, or a BZ boundary and point cloud from Phono3PyTools. PhonoVis alone converts those layers to PyVista meshes and owns backgrounds, camera presets, axes, scalar bars, screenshots, rotating GIFs, and interactive widgets. The domain commands remain `cat-freqball-view` and `cat-phono3py-cloud`; this keeps the interface close to the scientific task without duplicating rendering code.

StructureTools supplies a `StaticDisplacement` model after it maps two structures. PhonoVis reuses its atom, bond, polyhedron, cell, arrow, and camera machinery but replaces the phonon phase control with a physical displacement-scale slider. At scale one, arrows and displaced atoms represent the actual mapped structural difference.

## Mode and supercell convention

One frozen `Mode` dataclass contains the exact q-point, frequency, complex eigenvector, primitive structure, symbols, and masses. An atom in translated cell \(T\) moves as

\[
u_{jT}(\phi) = A\,\Re\left[\frac{e_j}{\sqrt{m_j}}
\exp\left(i(2\pi q\cdot T + \phi)\right)\right].
\]

This keeps finite-q phase relations correct in repeated cells. `mass_weighting: false` shows the raw dynamical-matrix eigenvector instead. `normalization` may be `max`, `unit`, or `none`. `repeat_from_qpoint: true` derives the smallest axis-aligned repeat commensurate with q; a manual repeat is often clearer during exploration.

Bonds and polyhedra are derived once from Pymatgen periodic neighbors. Their vertices retain atom indices and lattice-image translations, so coordination objects crossing a cell boundary remain complete and deform with the animated atoms. `boundary_images: true` displays the required image atoms. The same topology feeds Matplotlib, PyVista, GIF, and VMD coordinate generation.

## Export contract

- PDF and SVG are true Matplotlib vector graphics by default. PNG is available when requested. `export.backend: pyvista` instead uses PyVista sphere/cylinder rendering; its GL2PS PDF/SVG exporter remains vector when `raster: false`.
- GIF is one periodic mode cycle. The Matplotlib backend is deterministic and works headlessly; the optional PyVista backend produces sphere/cylinder rendering.
- VMD export is exactly one multi-frame XYZ trajectory plus one short `.vmd` loader. Run it with `vmd -e exported-mode.vmd`.
- Export names encode system, q-point, and zero-based band unless `export.stem` is set.
- One `manifest.yaml` records the input hashes, exact frequency, scene settings, and chosen visual files. It is provenance, not a duplicate numerical data format.
- No mode CSV is written. Use the canonical `band_modes` or `qpoint_modes` artifact from PhonoPyTools when numerical data are needed by CATPlot.

Static vector figures show equilibrium atom positions with displacement arrows; GIF and VMD frames move the atoms through the real part of the complex mode. Increase `scene.repeat` and `scene.amplitude_A` for presentation views without changing the underlying eigenvector.


## Recorded relaxation trajectories

PhonoVis can also display actual ionic geometries, independently of phonon-mode YAML or force constants:

```bash
cat-tools vasp relaxation --out relaxation-view
cat-phonovis trajectory relaxation-view/relaxation.extxyz
MPLBACKEND=Agg cat-phonovis trajectory relaxation-view/relaxation.extxyz \
    --gif relaxation-view/relaxation.gif --fps 6
```

Use CAT's `structure` and `visualization` extras. The exporter obtains species from a matching `POSCAR` beside `OUTCAR` (`--poscar` overrides it). It writes extended XYZ with cell vectors and forces, plus a plain XYZ and `.vmd` loader for `vmd -e relaxation-view/relaxation.vmd`. The [VMD XYZ reader](https://www.ks.uiuc.edu/Research/vmd/plugins/molfile/xyzplugin.html) accepts the coordinate trajectory; the loader supplies per-frame cell lengths/angles separately. Full lattice orientation remains in extxyz.

The Matplotlib viewer/GIF draws element-colored atoms and the cell with fixed bounds, labeling the recorded step, volume and maximum force. No sinusoidal motion, amplification, interpolation or physical timestamps are added. Periodic crossings remain as recorded. `--stride N` on the exporter reduces frame count while keeping the last complete block. Outputs are replaced on rerun; an incomplete final force block is omitted with a warning. This is a snapshot viewer, not a live `OUTCAR` watcher. A one-frame input yields a still GIF.

For split relaxations, `cat-tools vasp relaxation --gif` automatically joins a local `OUTCAR.sN` archive chain followed by `OUTCAR`, or numbered `relax_N` directories when there is no local chain. Sorting is numeric. PhonoVis labels segments and local steps, while the overall step index remains continuous. The extxyz and CSV retain source paths. Duplicate restart-boundary geometries are preserved. Each directory needs a matching `POSCAR`, or supply a shared `--poscar PATH`.
