# Mode Grüneisen parameters from three harmonic volumes

**Use:** Calculate mode-resolved volume derivatives on a mesh or band path.

**Files:** Supply three harmonic calculation directories, each with its own primitive `POSCAR`, matching `FORCE_SETS` or `FORCE_CONSTANTS`, and optional compatible `BORN`. Follow [volume preparation](../README.md) first.

## 1. Prepare matched volumes

Use a central volume and small expansion/contraction about it. Start with a trial ±1–2% volume interval and converge the derivative against a second interval; these are trial values, not universal accuracy settings. At each volume relax the allowed internal coordinates and run the same supercell/displacement construction.

For the examples below, `v0`, `vplus`, and `vminus` are your result directories. Put the actual relaxed unit/primitive structure used for each force calculation at `POSCAR` inside that directory. All three calculations must have the same atom order, mass assignment, supercell matrix and primitive mapping. Do not use the same equilibrium `POSCAR` for all three force datasets.

## 2. Calculate a mesh or path

Run from their parent directory. The native argument order is **central, expanded, contracted**, independent of directory names:

```bash
phonopy-gruneisen v0 vplus vminus --dim="2 2 2" \
  --pa="1 0 0 0 1 0 0 0 1" --mesh="21 21 21" --gc
```

This writes `gruneisen_mesh.yaml` in the locked Phonopy 4.4 environment. The example assumes the supplied cells are already primitive and the forces used a 2×2×2 supercell. Replace both settings with those actually used. Add `--readfc` for text `FORCE_CONSTANTS`; `--readfc --hdf5` reads `force_constants.hdf5`. Add `--nac` only with a valid `BORN` in each directory. The native Grüneisen CLI reads these explicitly named files; a CAT YAML is not accepted as its configuration.

Run the band calculation in a different output directory to preserve the mesh file:

```bash
mkdir band
cd band
phonopy-gruneisen ../v0 ../vplus ../vminus --dim="2 2 2" \
  --pa="1 0 0 0 1 0 0 0 1" \
  --band="0 0 0  0.5 0 0  0.5 0.5 0" --band_points=101 -p -s
```

Replace this illustrative reciprocal-fractional path for your crystal. Keep copies of the command, versions and input hashes with the native YAML and plot. The band calculation writes `gruneisen_band.yaml` in this environment. Check `phonopy-gruneisen --help` for the installed version's plotting options.

## 3. Interpret and converge

The mode parameter is −d ln(ω)/d ln(V). Acoustic zero modes and soft modes near zero can have ill-conditioned derivatives; inspect them separately. Near crossings or degeneracies, do not subtract frequency rows by sorted band index yourself. Phonopy's Grüneisen machinery works with dynamical matrices and eigenvectors. Compare the central dispersion, chosen volume interval and primitive convention before relating these values to another dataset.

A signed Grüneisen parameter is not a linewidth, lifetime or conductivity. QHA captures volume dependence of harmonic modes; [FC3 Grüneisen/scattering](../../02_Anharmonic/01_Phono3Py/README.md) uses a different calculation route. Agreement is a useful convergence check only when structures, strain/relaxation convention and NAC treatment match.

Use the native plotting path for this output. CAT's current Grüneisen reader is primarily tested on Phono3Py output; do not assume every nested band-YAML layout from `phonopy-gruneisen` is interchangeable with it.

Source: [Phonopy Grüneisen documentation](https://phonopy.github.io/phonopy/gruneisen.html). This replaces the outline-only legacy `07_Grüneisen` workflow with the full procedure.
