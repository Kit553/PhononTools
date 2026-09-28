# Baseline Harmonic Workflows

Start from the accepted structure and numerical settings in
[Common](../../00_Common/README.md).

| Route | Use | Handoff |
| --- | --- | --- |
| [Phonopy finite displacements](01_PhonoPy/README.md) | Generate supercells, calculate VASP forces and fit harmonic force constants | Displacement metadata and FORCE_SETS, with standard bands, DOS and thermal outputs |
| [VASP DFPT](03_DFPT/README.md) | Obtain harmonic force constants from a perfect supercell | Matching metadata and FORCE_CONSTANTS; the same post-processing |
| [Γ-point phonons](03_DFPT/VASP/README.md) | Independently compare primitive-cell Γ frequencies | Comparison with finite-displacement frequencies |

The finite-displacement route migrates `OLD_Workflows/02_PhonoPy`, with unchanged
source files under its VASP `Legacy` folder. DFPT builds on the already migrated
`09_DFPT-SinglePoints` material; its original files now reside in the harmonic DFPT archive.

Converge electronic k-point sampling, displacement amplitude, supercell size and
phonon q-mesh separately. A dense post-processing mesh cannot repair short-ranged
force constants from an insufficient supercell. Use consistent cell/atom mapping, masses, NAC and symmetrization settings
when comparing results.

Directional DOS/FreqBall, fatbands and mode mapping are outside this baseline pass.
Their existing material is retained for later work.
