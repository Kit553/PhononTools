## Convergence:

Prior to any investigation, the convergence of your calculations should be tested to not only make sure that the reported values are stable but also to optimize the timing of your calculations. Phonon Workflows are especially susceptible to convergence problems.
Typical Convergence issues discussed here:
  - Basis Set Completeness via the `ENCUT` parameter
  - Brillouin Zone Integration via *k*-meshes (electronic PAW) and *q*-meshes (phonon wavevector)
  - Electronic Convergence and Ionic Force Convergence via `EDIFF` and `EDIFFG`
  - Finite Size effects and Supercell Construction

## Files and helpers:
