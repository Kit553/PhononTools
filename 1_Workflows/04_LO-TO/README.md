## LO-TO Splitting
**Use:** Obtain the dielectric tensor $\epsilon$~$\infty$~ and the dielectric charges *Z**.  
         These can then be used to write the BORN files for Non-analytic term correction or elucidate LO-TO splitting in the material

## Workflow
**Folder Setup:**
- `POSCAR`  &rarr; fully relaxed structure of your material; no supercell
- `POTCAR`  &rarr; corresponding potential file
- `KPOINTS` &rarr; *k*-point sampling grid; it is recommended to use a very fine odd numbered sampling grid e.g. double the KPOINTS per direction
- `INCAR`   &rarr; job parameters; be sure that the run is A) Static and B) `LEPSILON` is set to .TRUE. to compute the dielectric tensor
- `job.sh`  &rarr; Job starter; see the uploaded file for recommended timings on PALMA-II

---

#### 1. Perform a static run
As VASP has to calculate the dielectric tensor in this run, it can be far longer than a simple static SCF calculation. While an odd numbered grid is
recommended to capture the behavior of the material around $\Gamma$ correctly, even numbered grids can be used to prevent visual weirdness around the
zone center.

---

#### 2. Writing the BORN file
After this run has concluded the BORN file can be written with PhonoPy using the created `OUTCAR` or `vasprun.xml`. 

```python
  phonopy-vasp-born > BORN
```

While it is recommended to write out the symmetrized BORN file, the raw BORN file can also be written out using the `--nost` tag, this is only useful
when trying to debug or for low-symmetry/defective structures.

```python
  phonopy-vasp-born --nost > BORN_NoSymm
```

---

## Theory:
\[1\] [Gonze, X.; Lee, C. Dynamical matrices, Born effective charges, dielectric permittivity ten-
sors, and interatomic force constants from density-functional perturbation theory. Phys.
Rev. B 1997, 55, 10355–10368.](https://doi.org/10.1103/PhysRevB.55.10355)

\[2\] [Gonze, X.; Vigneron, J.-P. Density-functional approach to nonlinear-response coefficients
of solids. Phys. Rev. B 1989, 39, 13120–13128.](https://doi.org/10.1103/PhysRevB.39.13120)
