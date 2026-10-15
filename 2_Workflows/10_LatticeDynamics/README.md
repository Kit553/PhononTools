# Lattice dynamics

**Use:** Follow harmonic, anharmonic, quasi-harmonic and mode-analysis calculations from prepared structures to interpretable results.

**Files:** [Common preparation](../00_Common/README.md) · [legacy coverage](../COVERAGE.md) · [analysis reference](../../3_Analysis/README.md).

| Task | Maintained procedure | Main result |
| --- | --- | --- |
| Harmonic phonons | [Phonopy/VASP](01_Harmonic/01_PhonoPy/VASP/README.md) | Frequencies, modes, DOS, thermal properties, displacement tensors |
| Independent phonon check | [DFPT](01_Harmonic/03_DFPT/VASP/README.md) | Hessian/commensurate phonons with matched conventions |
| Directional spectra and character | [FreqBall](01_Harmonic/02_DirectionalDOS/README.md) | Balls, projections, conduction fatbands, descriptors and mass scans |
| Three-phonon scattering and transport | [Phono3Py](02_Anharmonic/01_Phono3Py/README.md) | FC3, Γ/lifetimes, RTA/LBTE/Wigner κ, JDOS, PP, self energies |
| Corrected-input / incomplete GP recovery | [Audit and collection](02_Anharmonic/01_Phono3Py/recovery.md) | Temperature-complete JDOS and separate boundary-MFP collections |
| Volume-dependent modes | [Three-volume Grüneisen](03_QuasiHarmonic/01_Grüneisen/README.md) | Mode Grüneisen mesh/path |
| Thermal expansion and Cp | [QHA](03_QuasiHarmonic/02_ThermalProperties/README.md) | Equilibrium V(T), expansion, bulk modulus and Cp |
| Eigenvector inspection | [PhonoVis](04_ModeAnalysis/02_ModeViewing/README.md) | Interactive modes, figures, GIF and VMD |
| Frozen-mode energy landscape | [Mode mapping](04_ModeAnalysis/01_ModeMapping/README.md) | Validated 1D/2D PES from static DFT |

Begin with converged structure, potentials, electronic settings and a clear primitive cell. Born/NAC, masses and force data must refer to that same model. A shared table key is not proof that two independently calculated modes are physically identical. Each guide distinguishes installed CAT features, upstream commands and external calculations. Native outputs are retained when CAT has no matching extractor.
