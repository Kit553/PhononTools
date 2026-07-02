## Phonon Frequencies from DFPT
**Use:** Obtain the frequencies of all bands at $\Gamma$ from DFPT to validate frozen phonon approach.

## Workflow:
**Folder Setup:**
- POSCAR &rarr; relaxed primitive cell of your system
- POTCAR &rarr; corresponding potentials
- KPOINTS &rarr; choose the same settings that give good electronic convergence
- INCAR &rarr; make sure you set the correct `IBROIN = 8` setting for a DFPT run
- job.sh &arr; job starter for PALMA-II; see file for rough memory estimates and timing

### 1. Run SCF Calculation:

Make sure to modify your `INCAR` accordingly:

```bash
  SYSTEM = gamma_dfpt_phonons
  NWRITE = 3
  ALGO = Normal

  ENCUT = 520
  EDIFF = 1E-8
  ISMEAR = -5
  PREC = Accurate
  LASPH = .TRUE.

  IBRION = 8       
  NSW = 1          
``` 

If non-analytic term correction is desired, additionally set `LEPSILON = .TRUE.` but be aware of the memory and timing requirements this opens up!

`IBRION = 7` is recently depreceated in VASP so only pertubation theory with symmetry is considered stable, make sure to consider this in further comparisons with `PhonoPy`!

In general `NSW` should always be set to 1 when doing DFPT phonon calculations, `NSW = 0` is unsafe as the static behavior might override the request for the Hessian; higher values waste time since the structure should already be at a stable minimum.

---

It is possible to obtain the 'full' phonon dispersion from DFPT, however, this requires calculations on increasingly large supercells and subsequent Fourier interpolation making the comparison to frozen phonon cumbersome; it is recommended to only look for a comparison around $\Gamma$ so that the primitive cell of your system is sufficient. If you desire to get zone boundary points, look up the smallest commensurate supercell, perform the calculations the same then adjust the `--qpoints` settings of the last step; be careful!

---

### 2. Convert DFPT Hessian to `PhonoPy` readable input:

After the calculation is done, your `OUTCAR` should contain the section:

```bash
Eigenvectors and eigenvalues of the dynamical matrix
```

Using the `phonopy_disp.yaml` from your normal `PhonoPy` run this output can be converted into `PhonoPy`'s format with the command:

```bash
phonopy --fc vasprun.xml
```

For this to work without issues, rename the original `POSCAR` to something else to prevent `PhonoPy` from taking it as the structure input, we want it to use the `phonopy_disp.yaml`

The frequencies at $\Gamma$ can then be checked with:

```bash
phonopy-load --readfc --qpoints="0 0 0" -p
```

Make sure that `--qpoints` is correctly specified to only be $\Gamma$ for primitive cell DFPT.
