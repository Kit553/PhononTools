# Dielectric Tensor and Born Effective Charges

**Use:** Obtain dielectric and Born tensors for non-analytical phonon corrections.

**Files:** [INCAR.dielectric](VASP/INCAR.dielectric) · [submit.sh](VASP/submit.sh).

- [VASP procedure](VASP/README.md): response calculation, convergence and `BORN` export.
- [Original LO–TO workflow](VASP/Legacy/README.md): preserved source material.

Use the same accepted structure and physical model as the phonon calculation. This stage supplies response tensors; the lattice-dynamics workflow applies them to force constants and evaluates the directional limit near Γ.
