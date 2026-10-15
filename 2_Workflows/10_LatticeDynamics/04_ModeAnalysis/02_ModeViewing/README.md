# Inspect and export phonon modes

**Use:** Explore a dispersion interactively, inspect an exact eigenmode, and export figures or animations with the same geometry and mode convention.

**Files:** [phonovis.yaml](phonovis.yaml) · [PhonoVis reference](../../../../3_Analysis/docs/phonovis.md) · [Phonopy analysis example](../../../../3_Analysis/examples/phonopy.yaml).

## 1. Connect harmonic inputs

Copy the YAML to your work directory, set `system.path`, the native Phonopy YAML and either `FORCE_SETS` or force constants. Remove `born` if NAC is not intended; otherwise supply matching response data. Install the `visualization` extra. All YAML paths are relative to the YAML, with force filenames under `system.path`.

Set a material-appropriate reciprocal-fractional path explicitly, or install the optional Seekpath dependency for automatic paths. Choose an initial zero-based band and a repeat large enough to show the periodic pattern. Configure bonds and polyhedra with actual species and verified cutoffs rather than copying another material's coordination assumptions.

## 2. Inspect the motion

```bash
cat-phonovis explore phonovis.yaml
cat-phonovis explore phonovis.yaml --q 0.5 0 0 --band 3
cat-phonovis explore phonovis.yaml --q 0.5 0 0 --band 3 --viewer matplotlib
```

The dashboard couples a clickable dispersion, projected DOS and animated structure. Use space to pause, arrow keys to change point/band, and +/- to change display amplitude. Explicit mode viewing also offers the PyVista phase/amplitude controls. A larger display amplitude does not change the calculated frequency or PES.

`repeat_from_qpoint: true` derives a commensurate diagonal repeat. A manual repeat can show a local section but should be commensurate when presenting a periodic supercell. The mass-weighted complex eigenvector and Bloch phase determine motion; `normalization: max` is a visual scale, not a mass-normalized PES coordinate.

## 3. Add directional/scattering information

Generate [FreqBall](../../01_Harmonic/02_DirectionalDOS/README.md) tables with the identical path, primitive convention, masses, NAC and band-connection settings. Under `explore.dispersion`, set `color_by: direction_score` and add:

```yaml
joins:
  - path: freqball-analysis/freqball_band_modes.csv
    columns: [direction_score]
    require_complete: true
```

Ordinary `<projection>_participation` columns also work. Path joins validate segment, point, q-vector and band identity; do not interpret identical indices as evidence of mode equivalence between different physical calculations. Phono3Py linewidths can be joined on a compatible exact mesh path. They are not available at arbitrary interpolated band points merely because the harmonic dispersion is smooth.

## 4. Export

```bash
cat-phonovis export phonovis.yaml --q 0.5 0 0 --band 3
```

The supplied config uses the Matplotlib backend for headless PDF/SVG and GIF, plus a VMD XYZ trajectory and loader. Open the generated `.vmd` with `vmd -e FILE.vmd`. The manifest records input hashes, frequency and scene settings. Native mode tables come from `cat-phonopy`; PhonoVis exports visual files without duplicating them. PyVista screenshots on cluster nodes require a working VTK rendering backend.

At exact Gamma in polar materials, an explicit directional NAC limit may differ from a path endpoint. The current exact-mode viewer does not expose a separate NAC-direction CLI flag. Inspect the path endpoint when that limit matters and do not silently equate it to a direction-free Gamma view.
