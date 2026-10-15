# Anharmonic lattice dynamics

**Use:** Calculate third-order force constants, scattering spectra and lattice thermal conductivity, and connect them to harmonic mode character.

**Files:** [Phono3Py procedure](01_Phono3Py/README.md) · [Recovery and boundary-length studies](01_Phono3Py/recovery.md) · [Frozen-mode PES analysis](../04_ModeAnalysis/01_ModeMapping/README.md).

Phono3Py provides perturbative three-phonon scattering, RTA/LBTE transport and version-dependent Wigner transport. CAT prepares grid-point jobs, audits completion, generates recovery/collection jobs, and extracts canonical analysis tables. A frozen-mode PES probes selected anharmonic directions but does not by itself supply a renormalized finite-temperature spectrum or conductivity.

The maintained route below covers the scope listed in legacy `03_Phono3Py`, whose calculation section was unfinished. ALAMODE, self-consistent phonons and higher-order force-constant models remain external methods; CAT currently has no general adapter for them. An unstable harmonic reference is a reason to establish the appropriate physical treatment before interpreting perturbative lifetimes or conductivity.
