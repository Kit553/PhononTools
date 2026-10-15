# Quasi-harmonic workflows

**Use:** Obtain volume-dependent phonons, mode Grüneisen parameters, and equilibrium thermal expansion from completed harmonic calculations at several volumes.

**Files:** [Three-volume Grüneisen procedure](01_Grüneisen/README.md) · [QHA thermal properties](02_ThermalProperties/README.md) · [Harmonic force workflow](../01_Harmonic/01_PhonoPy/VASP/README.md).

| Question | Route | Required calculations |
| --- | --- | --- |
| How does each phonon respond to volume? | Three-volume Phonopy Grüneisen | Matched harmonic force constants at V0, V+ and V− |
| How do V, expansion, bulk modulus and Cp change with temperature? | Phonopy QHA | Static energies and harmonic thermal properties over a volume series |
| What are the lifetimes or thermal conductivity? | [Phono3Py](../02_Anharmonic/01_Phono3Py/README.md) | Third-order force constants and scattering calculations |

CAT supplies cell rescaling, harmonic analysis and plotting of its own tables. `phonopy-gruneisen` and `phonopy-qha` are upstream executables installed with Phonopy; there is no `cat-qha` command. Use their native outputs for the fitting steps. The procedures below are checked against the locked Phonopy 4.4 CLI, while the linked upstream documentation may describe a newer version.

## Shared preparation

Start from a converged primitive structure. Keep the exchange-correlation model, potentials, cutoff, electronic convergence, smearing strategy, magnetic state, primitive basis and force-displacement settings consistent across the series. Converge the electronic k mesh and phonon q mesh separately.

For isotropic volume scaling, from a new working directory containing `POSCAR`:

```bash
cat-tools vasp rescale POSCAR --iso -2 --dir
cat-tools vasp rescale POSCAR --iso 0 --dir
cat-tools vasp rescale POSCAR --iso 2 --dir
```

This creates `iso_m2pct/POSCAR`, `iso_p0pct/POSCAR`, and `iso_p2pct/POSCAR`. `--iso 2` changes volume by 2%, not each vector length by 2%. Relax internal coordinates at each fixed cell (`ISIF=2` with an ionic optimizer), then run the [static reference](../../00_Common/03_StaticReference/README.md) and harmonic force workflow at that cell. A full variable-cell relaxation would erase the scan.

Fixed cell shape under isotropic scaling is a constrained QHA approximation. Materials with substantial anisotropic strain need a separately converged fixed-volume shape-relaxation or strain-dependent free-energy treatment. A c/a scan alone does not implement anisotropic QHA. Track phase changes or internal rearrangements rather than fitting across unrelated structures.

A volume series needs compatible site ordering and primitive transformations. For NAC, calculate consistent Born/dielectric data for each volume, or explicitly state and assess an approximation using fixed response tensors. Never silently mix `BORN` files from different structures. Imaginary branches require investigation; discarding negative modes does not establish a valid harmonic free energy.
