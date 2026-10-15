# Maintained workflow coverage

This inventory compares scientific tasks with `OLD_Workflows`, including topics whose old guide was only an outline. A documented native route is distinguished from a CAT wrapper. Local synthetic checks do not establish production DFT accuracy or replace an actual cluster/software validation.

| Legacy topic | Maintained route | Reproduced scope and improvements |
| --- | --- | --- |
| Root `INCAR` notes | [Setup](00_Common/00_SettingUp/VASP/README.md) | Electronic/model settings, parallelization and explicit benchmark procedure; corrected spellings and qualified defaults |
| `00_Convergence` | [Convergence](00_Common/01_Convergence/VASP/README.md) | Cutoff/k sampling, precision, force/stress checks and recorded acceptance criteria |
| `01_Relaxation` | [Relaxation](00_Common/02_Relaxation/VASP/README.md), [static reference](00_Common/03_StaticReference/README.md) | Input preparation, restarts, geometry comparison, trajectory inspection and final static verification |
| `02_PhonoPy` | [Harmonic forces](10_LatticeDynamics/01_Harmonic/01_PhonoPy/VASP/README.md), [analysis reference](../3_Analysis/docs/phonopy.md) | Displacements, force collection, bands/DOS/projections, velocities, thermal properties, displacement tensors and mode data |
| `03_Phono3Py` | [Anharmonic workflow](10_LatticeDynamics/02_Anharmonic/01_Phono3Py/README.md) | Completes the old outline: FC2/FC3, distributed sampling, Γ/lifetimes, RTA/LBTE/Wigner, JDOS, Grüneisen, PP and native self-energy spectra |
| `04_LO-TO` | [Born response](00_Common/04_DielectricBorn/VASP/README.md) | Dielectric/Born calculation, extraction and matching primitive/NAC conventions |
| `05_Lobster` | [LOBSTER](20_ElectronicStructure/01_BondingAnalysis/01_LOBSTER/README.md) | Wavefunction/basis preparation, periodic contacts, quality checks, populations, DOS, charges, FATBAND and explicit advanced-native-output routes |
| `06_DirectionalDOS` | [FreqBall](10_LatticeDynamics/01_Harmonic/02_DirectionalDOS/README.md), [standalone](10_LatticeDynamics/01_Harmonic/02_DirectionalDOS/STANDALONE.md), [BVPA](99_Misc/README.md) | Direction families, moments, curves and balls through the compact tensor; adds shared projections, fatbands and mass scans |
| `07_Grüneisen` | [Three-volume Grüneisen](10_LatticeDynamics/03_QuasiHarmonic/01_Grüneisen/README.md) | Completes the old outline; distinguishes volume derivatives from FC3 Grüneisen |
| `08_ModeMapping` | [Frozen-mode PES](10_LatticeDynamics/04_ModeAnalysis/01_ModeMapping/README.md), [mode viewing](10_LatticeDynamics/04_ModeAnalysis/02_ModeViewing/README.md) | Modern 1D/2D generation, commensurability and coordinate checks, static energy collection, tidy CSV and figures |
| `09_DFPT-SinglePoints` | [DFPT](10_LatticeDynamics/01_Harmonic/03_DFPT/VASP/README.md) | Hessian extraction and phonon comparison with correct cell/supercell interpretation |

## Additional procedures now connected

- [QHA thermal fitting](10_LatticeDynamics/03_QuasiHarmonic/02_ThermalProperties/README.md): ordered E(V)/thermal data, normalization, volume/temperature coverage and fit checks.
- [Phono3Py recovery](10_LatticeDynamics/02_Anharmonic/01_Phono3Py/recovery.md): corrected `BORN` dependencies, missing JDOS temperatures, multiple boundary lengths and array-task/log interpretation.
- [Electronic bands/DOS](20_ElectronicStructure/02_BandsAndDOS/README.md): separate integration mesh and band path, native/Pymatgen analysis and the boundary with CAT's LOBSTER-based electronic exporters.
- [Ordering](00_Common/05_Ordering/README.md): the existing generation, calculation, collection and hull workflow remains part of common preparation.

## Intentional compatibility differences and remaining boundaries

| Area | Current contract |
| --- | --- |
| FreqBall files/configuration | The native tensor replaces atom-wise historical HDF5; old results are not converted by renaming. Documented legacy YAML normalization and numerical comparisons are covered by tests. The archived v7 `--check`/auto-config-generator interface is not claimed to exist in CAT. Use maintained templates and the documented input checks. |
| ModeMap indices and coordinates | The maintained helper uses zero-based bands and a recorded mass-normalized real field on the mapping cell. JMS CLI bands are one-based and its historical Q normalization must be checked before combining data. |
| Native QHA | Installed Phonopy performs Grüneisen/QHA. CAT has no general QHA fitting command or guaranteed importer for every native Grüneisen YAML layout. Native data/plots remain part of the workflow. |
| Native LOBSTER extras | CAT exports the documented populations, DOS, FATBAND, charges and diagnostics. Advanced distributions, real-space matrices, fragment energies and related outputs remain native/external analyses; the workflow identifies the native outputs and the tools needed to analyze them. |
| Anharmonic models | Three-phonon perturbative transport and frozen-mode PES are supported routes. Self-consistent finite-temperature renormalization, general higher-order solvers and ALAMODE integration are not implemented CAT workflows. |
| Execution evidence | Examples distinguish locked local interfaces from legacy cluster commands. Licensed VASP/LOBSTER calculations and cluster MPI/Slurm validation require the actual installation/data. |

Historical files are retained as references. Corrections include the actual JMS postprocessor name, one-based JMS band input, explicit distinction between complex mode phase and PES dimensionality, basis-derived `NBANDS`, and the separation of native spectral functions from qualitative linewidth-based Lorentzian plots.

## Local verification (2026-10-06)

- The full test suite ran 190 tests successfully, with one optional CATPlot test skipped in the locked review environment. Both recipe tests also passed in the separate environment containing CATPlot.
- Seven frozen-mode tests check mass normalization against harmonic curvature, finite-q image phases and commensurability, independent 2D coordinates, figures, reference energies and rejection of changed structures or unsuitable static runs. Energy-reader tests use mocked VASP results; they are not DFT calculations.
- Synthetic harmonic datasets exercised Phonopy 4.4 thermal-file generation, mesh/band Grüneisen commands with a native PDF, and a seven-volume QHA fit. The mesh and band filenames are `gruneisen_mesh.yaml` and `gruneisen_band.yaml`, respectively, in this environment.
- A synthetic Si/O structure exercised the supplied LOBSTER seed and Pymatgen `INCAR` generation, including the basis-derived band count. LOBSTER itself was not executed.
- Workflow links, new Python lint and source-archive inclusion were checked. Licensed VASP/LOBSTER jobs, production-material convergence and native cluster Phono3Py 2.7/MPI/Slurm execution remain installation-specific validation.
