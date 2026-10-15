# StructureTools

StructureTools is the suite's periodic-structure interchange layer. It uses the real Pymatgen `Structure` as the in-memory object; it does not wrap it in another structure API. Dataclasses are reserved for concepts Pymatgen does not represent by itself, such as a trajectory frame with time, energy, forces, and velocities.

```bash
uv sync --extra structure
```

## Canonical representation

Two tables are sufficient for both a static crystal and an MD trajectory:

| Artifact | Grain | Contents |
|---|---|---|
| `structure_cells` | system × frame | full lattice, PBC, time, energy, temperature, charge, volume, source |
| `structure_sites` | system × frame × site × species | fractional/Cartesian coordinates, occupancy, labels, properties, velocities, forces |

The species component in the site key records each species and occupancy on disordered sites. The lattice occurs once per frame rather than once per atom. Together the tables round-trip labels, occupancies, oxidation-state species strings, structure and site properties, variable cells, forces, and velocities. CSV and Parquet are the only table representations; Parquet is the sensible choice for long trajectories.

`structure_summary` adds cell metrics, density, atom density, formula, and detected space group. `symmetry_sites` records zero-based site indices, Wyckoff symbols, and equivalence groups. These are derived reports, not competing structure formats.

## Structure conversion and tables

Pymatgen handles CIF, `POSCAR`/`CONTCAR`, CSSR, JSON and its other supported periodic formats. Conversion does not symmetrize unless requested explicitly:

```bash
cat-prepare structure convert POSCAR structure.cif
cat-prepare structure convert structure.cif POSCAR --standardize primitive

cat-structure normalize structure.cif -o structure-analysis \
  --system Na3PS4 --latex tables/Na3PS4-structure.tex

cat-prepare structure restore structure-analysis/structure_sites.csv \
  structure-analysis/structure_cells.csv restored.cif
```

`--standardize primitive|conventional|refined` is intentionally visible because standardization can change cell choice, atom count, and site order. Plain conversion preserves the input structure.

The Python boundary is equally direct:

```python
from cat_tools.structure import read_structure, structure_artifacts, static_frame

structure = read_structure("structure.cif")       # Pymatgen Structure
cells, sites = structure_artifacts([static_frame(structure, system="sample")])
```

## LaTeX crystallographic tables

The LaTeX exporter writes a compact `booktabs` table containing symmetry-unique sites, Wyckoff labels, fractional coordinates, occupancy, multiplicity, cell parameters, and space group. Use `--all-sites` when every expanded site is required.

```bash
cat-structure latex structure.cif structure.tex

# The same export directly from the canonical machine tables:
cat-structure latex structure-analysis/structure_sites.csv structure.tex \
  --cells structure-analysis/structure_cells.csv --frame 0
```

LaTeX is a presentation export, not a machine-readable authority. Reading arbitrary hand-edited LaTeX back into a structure would be ambiguous; table-to-structure uses the canonical site and cell tables instead.

## Structure comparison

```bash
cat-structure compare POSCAR-reference CONTCAR \
  -o structure-comparison --system Na3PS4
```

The matcher tolerates atom reordering and periodic origin shifts. It reports cell length, angle, volume, space-group, and composition changes separately from mapped internal atomic displacements. The outputs are:

- `structure_comparison`: one diagnostic row with RMS/mean/maximum displacement and concise threshold-based major differences;
- `site_displacements`: the reference-to-other site map and Cartesian displacement vector for every mapped site.

The default warning thresholds are 0.1 Å internal displacement, 1% lattice-length change, 1° angle change, and 3% volume change. These thresholds classify the report; they do not alter mapping or coordinates.

For interactive inspection, StructureTools creates a static-displacement model and hands it to PhonoVis:

```bash
uv sync --extra visualization --extra structure
cat-structure compare POSCAR-reference CONTCAR -o structure-comparison \
  --view --scene examples/structure-scene.yaml
```

PhonoVis owns atom, bond, polyhedron, cell, arrow, camera, and slider rendering. StructureTools contains no plotting backend. The displacement-scale slider multiplies the physical mapped displacement vectors; scale `1` is the actual structural change.

## ASE and future MD work

ASE supplies the open trajectory boundary, so the same command can normalize `.traj`, extended XYZ, and LAMMPS formats supported by the local ASE installation:

```bash
cat-structure trajectory run.traj -o structure-trajectory \
  --system Na3PS4-md --format parquet

cat-prepare structure ase-export structure-trajectory/structure_sites.parquet \
  structure-trajectory/structure_cells.parquet recovered.traj
```

ASE velocities are converted to Å/fs; energies and forces are retained when the input format carries a calculator result. StructureTools currently requires a fully periodic 3D ASE cell because its canonical object is a Pymatgen `Structure`, not a `Molecule`. That restriction is explicit and can later be relaxed by adding a molecular record without weakening the crystal schema.

## SoftBV/BVPA launcher

The platform-specific WSL wrapper belongs here because it launches a structure/grid calculation and manages its files. Parsing the resulting path topology and turning it into atom-local BVS/NaNa families remains FreqBallTools' responsibility.

```bash
cat-prepare bvpa \
  --cif structure.cif --cube bv.cube --ion Na1+ \
  --softbv-root /mnt/c/path/to/SoftBV \
  -o BVPA --conduction-path conduction_path.yaml
```

`SOFTBV_ROOT` can replace `--softbv-root`. BVPA 1.32 must run in the real SoftBV GUI root, so the wrapper temporarily protects pre-existing `BVPA_*` files there, moves the new outputs to the requested directory, and restores the old files. Standard output, standard error, and a small run-provenance YAML are kept with the result. Supplying `--conduction-path` calls FreqBallTools' parser directly; the resulting YAML is the input to `cat-freqball-families`.

Ready-to-edit CATPlot recipes under `recipes/structure` show site-displacement bars and MD cell volume versus time. Later MD analysis—MSD, van Hove functions, conductivity, and correlation measures—can consume these canonical frames without owning another `POSCAR`/XYZ parser.
