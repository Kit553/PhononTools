# Dielectric Tensor and Born Effective Charges

Obtain the electronic dielectric tensor and Born charges for non-analytical phonon
corrections in polar insulators.

- [VASP procedure](VASP/README.md): response calculation, convergence and BORN export.
- [Original LO–TO workflow](VASP/Legacy/README.md): preserved source material.

Use the same accepted structure and physical model as the phonon calculation.
This stage supplies response tensors; the lattice-dynamics workflow applies them
to force constants and evaluates the directional limit near Γ.
