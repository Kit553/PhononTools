# Harmonic phonon analysis

`cat-phonopy` reads completed harmonic force data and exports mesh/path/q-point tables, projected DOS, thermal properties and descriptors. For a single workflow covering directional spectra, plane normals, conduction fatbands, balls and mass scans, use [FreqBall](freqball.md); both tools use the common mode-table schema.

## Inputs and execution

Start from [the example YAML](../examples/phonopy.yaml), edit paths and numerical settings, and run from the installed CAT environment:

```bash
cat-phonopy phonopy.yaml
cat-phonopy phonopy.yaml --out analysis-check --format parquet
```

The workflow YAML is a CAT configuration, distinct from the native Phonopy YAML named in `system.phonopy_yaml`. Paths resolve relative to the CAT configuration; force-input filenames resolve relative to `system.path`. Use compatible primitive metadata and force sets/constants. Enable `system.born` only for the corresponding structure, with the intended non-analytical correction.

Calculation sections set to `null` are disabled. `mesh`, `band` and `qpoints` export mode tables; `dos`, harmonic thermal calculations and neutron calculations require a `mesh` section. Set `mesh.export: false` when the mesh is needed internally but its mode table is not wanted. `band.paths` contains reciprocal-fractional vertices; omitting it requests the automatic path and its optional Seekpath dependency.

## Projections and modes

Each mode/DOS section has its own `projections` mapping. Element projections default to enabled. Optional `sites` and `symmetry_sites` add primitive-site or symmetry-orbit groups. `directions` maps names to Cartesian three-vectors; `custom` combines a name, atom/element selection and optional Cartesian direction. Example:

```yaml
band:
  paths: [[[0, 0, 0], [0.5, 0, 0]]]
  labels: [G, X]
  npoints: 101
  projections:
    elements: true
    custom:
      - name: Na_x
        elements: [Na]
        direction: [1, 0, 0]
```

The path is illustrative; choose one appropriate to the actual primitive cell. PhonoPyTools directions are Cartesian. The richer `{fractional: [...]}` and `{miller: [...]}` direction syntax belongs to FreqBall's top-level `projections`. Do not transfer those mappings into this interface.

Rows identify zero-based q-point/band indices, reduced q-coordinates and frequency in THz. Band rows also contain segment/point indices, cumulative distance and a trace key for plotting disconnected branches separately. Projection columns end in `_participation`. Default `weighting: eigenvector` measures fractions of the full eigenvector norm; mode-table `weighting: displacement` instead divides atomic weights by masses before normalization. The participation ratio is computed with the same selected weighting convention. DOS has its own projection normalization.

Group velocities use Phonopy's native length-times-THz units. For VASP inputs with length in Å, multiply by 100 for m/s. `velocity_directions` projects them onto named Cartesian directions. Symmetry-reduced mesh representatives carry weights, but directional values are not automatically rotated onto omitted q-points. Use a full mesh when interpreting fixed Cartesian or atom-local directions.

## Derived tables and provenance

The example demonstrates Gaussian DOS, harmonic thermal properties and thermal displacement matrices. `descriptors.averaged` combines DOS moments with available mesh-mode information; `group_velocity` and `sound_velocity` select velocity summaries. The latter uses near-Gamma sampling and needs a stated native-unit to m/s conversion. Configure frequency cutoffs explicitly when inverse-frequency measures or unstable/near-zero modes matter.

The optional neutron calculation requires explicit scattering lengths, q-points, temperature and basis choice; CAT does not infer isotope-specific scattering data. Thermal-matrix CIF exports use the configured temperatures.

Each run writes one CSV or Parquet table per requested concept plus `manifest.yaml` with keys, units, settings and hashes. A rerun replaces files owned by the previous manifest; use a different output directory to preserve an earlier variant. Figures consume these tables via [Phonopy recipes](../recipes/phonopy/) or the interactive [PhonoVis workflow](phonovis.md).

Validation includes `test_phonopy_integration.py`, `test_modes.py`, `test_projections.py`, `test_descriptors.py` and `test_thermal_neutron.py`.
