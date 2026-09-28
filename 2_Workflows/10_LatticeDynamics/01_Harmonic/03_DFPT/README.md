# Harmonic Force Constants from VASP DFPT

Use DFPT as an alternative source of harmonic force constants, or as an independent
finite-displacement check. For a primitive-cell Γ-only comparison, follow the
[Γ-point phonon procedure](VASP/README.md).
Its [legacy archive](VASP/Legacy/README.md) preserves
`OLD_Workflows/09_DFPT-SinglePoints`.

For dispersions/DOS, this route calculates the perfect **supercell** and obtains
real-space force constants for interpolation. A primitive-cell Hessian alone is
not a converged full-zone dispersion. The calculation can require substantially
more memory/time than an ordinary SCF run.

## 1. Generate the matching supercell

Use the accepted common baseline and the same physical model as the phonons you
want to compare. In a new working directory, supply `POSCAR-unitcell` and run:

```bash
phonopy-init -d --dim 2 2 2 --pa auto -c POSCAR-unitcell
mkdir dfpt
cp SPOSCAR dfpt/POSCAR
```

The dimensions are an example to converge. Keep `phonopy_disp.yaml`, SPOSCAR and
the reported primitive mapping together. The generated displaced POSCARs are not
used for DFPT. Unlike finite-displacement collection, there is one perfect-supercell
calculation and no reference-force subtraction.

## 2. Prepare and submit the response calculation

Use the harmonic [INCAR.dfpt](VASP/INCAR.dfpt) and
[submission script](VASP/submit.sh) as references:

```bash
dfpt_templates="$PHONONTOOLS_ROOT/2_Workflows/10_LatticeDynamics/01_Harmonic/03_DFPT/VASP"
cp "$dfpt_templates/INCAR.dfpt" dfpt/INCAR
cp "$dfpt_templates/submit.sh" dfpt/submit.sh
```

Restore all converged electronic/model settings from your accepted baseline, retain
`IBRION=8`, `NSW=1`, and supply the matching POTCAR and a reviewed **supercell**
electronic KPOINTS mesh. Verify that the VASP build supports the selected model in
DFPT. Review rank layout and response-specific restrictions; do not copy the
finite-displacement SCF NCORE tuning automatically.

Benchmark the actual response workload and record it in the group's shared
[performance log](../../../00_Common/00_SettingUp/VASP/benchmark-log.md). The
36-rank, 32G, two-hour template is a pilot starting point; resize it for this
supercell before submitting.

```bash
cd dfpt
sbatch submit.sh
cd ..
```

## 3. Extract and analyse

Check SCF/response convergence, normal termination and a complete Hessian in
`dfpt/vasprun.xml`. Verify the calculated POSCAR still matches SPOSCAR. From the
generation root, extract:

```bash
phonopy-init --fc dfpt/vasprun.xml
```

Keep the resulting FORCE_CONSTANTS with this run's `phonopy_disp.yaml`; never borrow
metadata from a differently ordered or sized supercell. Follow the
[Phonopy VASP-DFPT procedure](https://phonopy.github.io/phonopy/vasp-dfpt.html).

For example, a total DOS calculation in an isolated analysis directory:

```bash
workflow="$PHONONTOOLS_ROOT/2_Workflows/10_LatticeDynamics/01_Harmonic/01_PhonoPy/VASP"
mkdir dos
cp phonopy_disp.yaml FORCE_CONSTANTS dos/
cp "$workflow/dos.conf" dos/settings.conf
cd dos
phonopy --config settings.conf --nonac phonopy_disp.yaml > phonopy.log 2> phonopy.err
cd ..
```

Use the other [standard post-processing configurations](../01_PhonoPy/VASP/README.md#4-standard-post-processing)
for bands, site DOS, thermal properties and displacement tensors. Optional NAC uses
the [common dielectric/Born workflow](../../../00_Common/04_DielectricBorn/VASP/README.md)
and matching mapping; use the same settings on both sides of any comparison.

Converge the supercell, electronic accuracy and q-mesh. Compare frequencies and
degenerate-mode subspaces at equivalent q-points using the same masses, units,
NAC and force-constant symmetrization conventions.
