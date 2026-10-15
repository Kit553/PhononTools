# Scientific analysis

Install the suite from the [repository root](../README.md). The command examples below assume your working directory is `3_Analysis`; adapt the example YAML input paths to your data. Preparation commands are documented in [1_Tools](../1_Tools/README.md).

The analysis code lives in [src/cat_tools](src/cat_tools/), with shared data formats and geometry imported from `cat_shared`.

`cat-tooling` consolidates scientific post-processing around stable, tidy tables that CATPlot can consume directly. Scientific modules calculate and export data; CATPlot provides declarative plotting. FreqBall also offers a small Matplotlib plotting command for its canonical tables, and PhonoVis handles interactive 3D rendering. A run writes one representation of each requested data concept plus a `manifest.yaml`; it does not produce competing wide, long, and plot-specific copies.

The shared artifact writer validates the complete batch before touching an existing result, serializes each table once in a sibling staging directory, and commits the manifest last. Each manifest records the primary key, units, metadata, row/column shape, and SHA-256 of every table. A later run removes only files owned by the previous manifest, so changing CSV to Parquet cannot leave a stale parallel representation and unrelated notes or figures remain untouched.

## PhonoPyTools

`cat_tools.phonopy` covers:

- mesh, path, and arbitrary-q phonons;
- element, site, custom-group, and Cartesian-direction mode projections;
- group velocities and participation ratios;
- total/projected DOS and full per-system averaged descriptors from both DOS and exact mesh modes;
- harmonic thermodynamic properties and thermal displacements;
- coherent one-phonon neutron dynamic structure factors;
- near-Gamma acoustic-velocity summaries.

Run a workflow with:

```bash
cat-phonopy examples/phonopy.yaml
```

Paths in the configuration are relative to the configuration file. Every run writes one table per data concept and one `manifest.yaml`; it never writes parallel wide and long representations of the same result.

Band and mesh indices are zero-based, matching the Phonopy Python API. Eigenvector projection columns use the suffix `_participation`. Group velocities remain in Phonopy's native length-times-THz unit and the unit is recorded in the manifest. See the [harmonic analysis reference](docs/phonopy.md) for projection conventions, dependencies between calculation sections, and the boundary with FreqBall.

## Phono3PyTools

`cat_tools.phono3py` reads completed Phono3py outputs without wrapping the Phono3py object or loading every large tensor into memory. Its canonical `mesh_modes` table has one row per dataset, selected temperature, irreducible q-point, and band. It contains the available linewidth, lifetime, group velocity, mean-free-path, heat-capacity, interaction-strength, modal-conductivity, Grüneisen, and joined projection columns.

```bash
uv sync
cat-phono3py examples/phono3py.yaml
```

The example workflow demonstrates weighted profiles, 2D mesh densities, direct kappa accumulation versus frequency and mean free path, approximate Lorentzian path maps, the harmonic PhonoPyTools bridge, JDOS, exact mesh paths, Wigner mode-pair exports, full interaction tensors, and external FreqBall column joins. See [`docs/phono3py.md`](docs/phono3py.md) for the schemas and design choices. Ready-to-use CATPlot configurations are in [`recipes/phono3py`](recipes/phono3py).

The experimental full-first-Brillouin-zone viewer is intentionally a separate command:

```bash
uv sync --extra cloud
cat-phono3py-cloud phono3py-analysis/mesh_modes.csv phono3py.yaml \
  --value gamma_THz --temperature 300
```

Phono3PyTools expands an irreducible scalar mesh by crystal symmetry and folds it into the true reciprocal Wigner-Seitz cell. It then hands a renderer-neutral point layer to PhonoVis, which opens the rotatable PyVista cloud. Its slider runs from all points at `0` to only the strongest points at `1`. Use `--no-symmetry-expansion` for direction-dependent values such as a conduction-direction weight unless the input already contains the correctly transformed full mesh.

## FreqBallTools

`cat_tools.freqball` replaces the old repeated direction-by-direction projected-DOS calculation with one exact six-component spectral tensor. From that single `freqball.h5` it derives frequency balls, one tidy directional-descriptor table, exact structural-family spectra, and named Cartesian/lattice/plane-normal projections. Mesh/path fatband columns use the same projection definitions with force-data eigenvectors; they are not reconstructed from the frequency-only tensor. BVPA parsing and BVS/NaNa/tetrahedral family construction are included. FreqBallTools creates the ball surface and family vectors; PhonoVis owns their 3D rendering.

```bash
cat-freqball run examples/freqball.yaml
cat-freqball plot fatband freqball-analysis/freqball_band_modes.csv -o fatband.pdf
cat-freqball view freqball.h5 --families direction_families.yaml
cat-freqball mass-scan setup examples/freqball.yaml -o mass-scan \
  --species S --range 16 96 --points 17
cat-freqball mass-scan collect mass-scan -o freqball-mass-analysis
```

Follow the separate [FreqBall workflow](../2_Workflows/10_LatticeDynamics/01_Harmonic/02_DirectionalDOS/README.md) for execution. See [`docs/freqball.md`](docs/freqball.md) for the schema, family-preparation commands, MPI calculation, fatband definitions, and the CATPlot recipes in [`recipes/freqball`](recipes/freqball).

## PhonoVis

`cat_tools.phonovis` is the dedicated 3D/exploration owner and eigenvector visual exporter. Its dashboard places a clickable, optionally participation-colored dispersion beside a projected DOS and an animated 3D mode. An explicit q-point mode opens in a rotatable PyVista view with phase and amplitude sliders.

```bash
uv sync --extra visualization
cat-phonovis explore examples/phonovis.yaml
cat-phonovis explore examples/phonovis.yaml --q 0.5 0 0 --band 7
cat-phonovis export examples/phonovis.yaml
```

The export workflow writes selected PDF/SVG/PNG figures, a GIF, and one VMD-compatible XYZ trajectory plus loader script, with one manifest for provenance. Periodic bonds and complete coordination polyhedra share one Pymatgen-derived topology across every renderer. PhonoVis does not fork the scientific data model: CATPlot continues to consume the canonical PhonoPyTools tables. See [`docs/phonovis.md`](docs/phonovis.md) for controls and the mode convention.

## BondingTools

`cat_tools.bonding` discovers completed LOBSTER runs and exports the data available in those files: ICOHP/ICOBI/ICOOP values and curves, interaction summaries, charges, DOS, electronic FATBAND dispersions, and run-quality diagnostics.

```bash
uv sync --extra bonding
cat-bonding examples/bonding.yaml
```

An optional LobsterPy benchmark independently constructs its `Analysis` object and compares element-pair means against direct extraction. The comparison is visibly marked diagnostic-only and never feeds final plotting. See the comprehensive [bonding README](src/cat_tools/bonding/README.md), including installation and plotting with [CATPlot](../4_Plotting/README.md), and the recipes in [`recipes/bonding`](recipes/bonding).

## StructureTools

`cat_tools.structure` is the canonical periodic-structure interchange layer. It reads and writes Pymatgen-supported crystal formats, converts structures to and from one cell table plus one site/species table, exports symmetry-aware LaTeX crystallographic tables, and compares reordered structures without confusing cell strain with internal atomic motion.

```bash
uv sync --extra structure
cat-structure normalize structure.cif -o structure-analysis --latex structure.tex
cat-structure compare POSCAR CONTCAR -o structure-comparison
cat-structure trajectory run.traj -o structure-trajectory --format parquet
```

ASE-readable trajectories use the same tables as static structures, including optional time, energy, temperature, velocity, and force columns. StructureTools prepares mapped displacements, while PhonoVis remains the sole 3D renderer. See [`docs/structure.md`](docs/structure.md) and the CATPlot recipes in [`recipes/structure`](recipes/structure).

The WSL-to-Windows SoftBV/BVPA launcher is also a StructureTools subcommand. It can immediately hand the completed calculation to FreqBallTools' canonical parser:

```bash
cat-prepare bvpa --cif structure.cif --cube bv.cube --ion Na1+ \
  -o BVPA --conduction-path conduction_path.yaml
```

## ThermalTools and ImpedanceTools

ThermalTools and ImpedanceTools provide measurement parsing and model fitting; their interfaces are under development. ThermalTools parses original PPMS thermal-transport `.dat` files and optionally fits the two-channel propagon/diffuson model against a PhonoPyTools DOS. ImpedanceTools discovers RelaxIS 3 exports by content, supports one- and two-process circuits, and emits spectra, parameters, diagnostic fit quality, first/total resistance, and Arrhenius fits.

```bash
cat-thermal examples/thermal.yaml
cat-impedance examples/impedance.yaml
```

See [`docs/thermal.md`](docs/thermal.md), [`docs/impedance.md`](docs/impedance.md), and their CATPlot recipes. DRT and further thermal models can extend these packages later without adding one-off plot scripts.

## Diffraction

TOPAS refinement parsing and CATPlot export live in `cat_tools.diffraction`:

```bash
cat-diffraction refinement.out -o refinement.csv --profile profile.txt --reflections reflections.txt
```

See [Diffraction documentation](docs/diffraction.md). The original `cat_tools.structure.refinement` import remains available for compatibility.

## Package layout

Build configuration, dependency locking, and tests live at the repository root. See [the original migration notes](CONSOL_MIGRATION.md) for source history.
