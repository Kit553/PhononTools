# VASP Γ-Point Phonons with DFPT

**Use:** Calculate primitive-cell Γ-point phonons using DFPT and, when available, compare them with a finite-displacement calculation. Start from the accepted [static reference](../../../../00_Common/03_StaticReference/VASP/README.md).

**Files:** [INCAR.dfpt](INCAR.dfpt) · [submit.sh](submit.sh) · [original inputs and notes](Legacy/README.md).

The example uses nonmagnetic bulk PBEsol. Confirm that your physical model is supported by VASP's DFPT implementation before using it.

## Prepare and run

Use a relaxed primitive cell with the same geometry, species order, potentials and physical model as the finite-displacement reference. A conventional-cell calculation also contains modes folded to its Γ point and is not the same comparison.

From the parent of `static_01`, provided its structure is that primitive cell:

```bash
mkdir dfpt_01
cp static_01/{POSCAR,INCAR,POTCAR,KPOINTS,submit.sh} dfpt_01/
```

Edit the copied `INCAR`, preserving your converged electronic settings:

```text
SYSTEM = gamma-dfpt
IBRION = 8
NSW = 1
LEPSILON = .FALSE.
```

Remove relaxation-only `EDIFFG` and `IWAVPR`; use [INCAR.dfpt](INCAR.dfpt) to review the remaining baseline. `IBRION=8` uses symmetry in DFPT; it computes the response rather than relaxing the starting structure. Tight electronic convergence matters especially for the Hessian! See [VASP's IBRION reference](https://vasp.at/wiki/IBRION).

Review the response run's resources and executable compatibility. Then:

```bash
cd dfpt_01
sbatch submit.sh
```

## Inspect the VASP result

After completion, check electronic/response convergence, the dynamical-matrix eigenvalues/eigenvectors in `OUTCAR`, and a complete Hessian in `vasprun.xml`. Retain imaginary modes and the three acoustic modes in the comparison. Investigate small acoustic residuals using stricter numerical settings before interpreting them.

Use the same masses, units and non-analytical-correction convention on both sides. For the initial comparison, disable NAC in both calculations. If subsequently including [Born/dielectric response](../../../../00_Common/04_DielectricBorn/VASP/README.md), use matching tensors and the same limit direction near Γ.

## Optional Phonopy conversion

Use the Phonopy 4 command split installed with this repository: `phonopy-init` prepares files and `phonopy` performs post-processing. Return to the parent of `dfpt_01` and create a fresh analysis directory with metadata from this calculation's own primitive cell:

```bash
mkdir dfpt_analysis
cd dfpt_analysis
cp ../dfpt_01/POSCAR POSCAR-unitcell
phonopy-init -d --dim 1 1 1 --pa P -c POSCAR-unitcell
phonopy-init --fc ../dfpt_01/vasprun.xml
phonopy --qpoints 0 0 0 --pa P --nonac phonopy_disp.yaml
```

Phonopy reads `FORCE_CONSTANTS` automatically in this directory. The generated displacement POSCARs are metadata-setup byproducts and need no VASP runs for this DFPT route. Inspect `SPOSCAR` against the DFPT `POSCAR` before using the extracted force constants; atom order and cell must agree. Keep `phonopy_disp.yaml` and `FORCE_CONSTANTS` together for post-processing. The [Phonopy VASP-DFPT guide](https://phonopy.github.io/phonopy/vasp-dfpt.html) documents Hessian extraction and the relation between the calculated cell and metadata.

Do not borrow metadata from an unrelated finite-displacement supercell or rename away inputs to let automatic discovery guess the mapping. This identity-supercell example supports the Γ comparison only, not a converged dispersion across the zone.

## Compare the phonons

Compare frequencies in common units, inspect eigenvectors when modes cross, and
treat degenerate modes as subspaces rather than requiring identical vector labels.
Check Phonopy's force-constant symmetrization/acoustic-sum-rule settings: processed
frequencies can differ from the raw VASP eigenvalues because of these operations.

If agreement fails,
check geometry, masses, mapping, numerical accuracy and the finite-displacement
supercell/amplitude. Agreement at Γ validates that comparison; a full stability
claim needs the relevant wave-vector and supercell checks.
