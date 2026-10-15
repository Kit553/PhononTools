# Ordering calculation acceptance record

Copy this into the calculation project and fill it while following the workflow. An empty field is an unperformed check, not an accepted default.

## Study and structures

- Date / operator / repository revision:
- Chemical substitution space and endpoint compositions:
- Parent sources, input symmetry and site selections:
- Atomic versus coupled-motif ensembles:
- Pressure / cell constraints / magnetic state / functional / other model terms:
- Expected energy scale (state eV/atom or meV/f.u. and conversion):
- Target uncertainty in relative and mixing energies (eV/atom):
- Force criterion (eV/angstrom) and residual-stress criterion (GPa):

## Numerical convergence pilots

Record actual meshes, energy definition, and stricter references. For static comparisons hold geometry fixed. Compare energies per atom; `EDIFF` is per cell.

| Candidate ID / geometry | Tested parameter values | Stricter reference | Change in relative/mixing energy (eV/atom) | Force-vector/stress-tensor differences | Decision |
| --- | --- | --- | --- | --- | --- |
| | | | | | |

- Potential release / per-element variants / source records and checksums:
- Common accepted `ENCUT` / `PREC` / `LREAL` / `LASPH` / `ADDGRID`:
- Sampling policy / centering / actual meshes for representative cells:
- `ISMEAR` / `SIGMA` / smearing sensitivity / energy definition:
- `EDIFF` / attained SCF convergence:
- Joint stricter cutoff + mesh + SCF validation:
- `EDIFFG` / residual forces and stress / restart stability:
- Revalidation after cell relaxation / changed cell sizes:
- Final physical/numerical protocol and shared energy_model label:

## Configuration and SQS convergence

| Batch / composition / phase | Cell matrix or enumeration sizes | Method / limits / seeds | Unique / selected / calculated | Correlation residuals (if SQS) | Change in hull / observable | Decision |
| --- | --- | --- | --- | --- | --- | --- |
| | | | | | | |

- SQS cluster cutoffs / residuals by shell / search effort / ensemble:
- Cell-size and cell-shape tests:
- Additional composition points and motif patterns:
- Selection expansion / uncomputed candidates / coverage limitations:
- Final symmetries / distinct starts that relaxed into the same structure:
- Symmetry-breaking checks and any separate phonon-stability evidence:

## Hull acceptance

- Included collections (one accepted record per candidate ID):
- Chosen endpoint IDs and evidence they are lowest among tested endpoint phases:
- Hull revision and newly added candidates:
- Changes in vertices / decompositions / mixing energies from previous revision:
- Near-ties unresolved within measured numerical uncertainty:
- Missing competing compounds / entropy terms / other scope limitations:
- Accepted for which physical claim, or outstanding work:
