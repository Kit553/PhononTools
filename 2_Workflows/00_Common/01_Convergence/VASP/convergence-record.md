# Convergence record: <system / phase / purpose>

Status: draft / provisional / validated / superseded Date: YYYY-MM-DD Author/contact: Related records / replacement:

## Scope and reproducibility

- Intended property and conditions (pressure, temperature, etc.):
- Structure source, phase, composition, atom count, cell dimensions/volume:
- Structure file or portable relative path; `POSCAR` SHA256:
- VASP version/build and executable; CAT version or repository commit:
- Functional, dispersion, spin/magnetic initialization, SOC, +U, charge:
- PAW release and variants in `POSCAR` order; `POTCAR` SHA256; max ENMAX (eV):
- Input/output archive location (and how another user can access it):
- Geometry used for static tests, including any displacement or strain:
- Restart policy; extra input files; physical model assumptions:

Keep licensed potential files out of a publicly shared record. Record their identifiers and hashes so input can be reconstructed.

## Acceptance criteria, chosen before testing

| Quantity | Allowed difference from reference | Reference / comparison method |
| --- | --- | --- |
| Energy per atom (state which energy) | <eV/atom> | |
| Largest atom-wise force-vector difference | <eV/angstrom> | |
| Largest stress-component difference | <kbar> | |
| Relaxed geometry / residual stress | <criterion and units> | |
| Target property | <criterion and units> | |

## Selected settings

These are starting suggestions for other systems; state below what was actually validated.

```text
ENCUT = <eV>
PREC = <value>
EDIFF = <eV>
ISMEAR = <value>
SIGMA = <eV>
NBANDS = <actual OUTCAR value; or default plus actual value>
LREAL = <value>
ROPT = <per-species values, defaults with PREC, or not applicable>
ADDGRID = <value>
LASPH = <value>
IBRION / NSW / ISIF / EDIFFG = <relaxation settings, if applicable>
```

- `KPOINTS`: Gamma/Monkhorst-Pack, actual mesh, shift, CAT density if used:
- Actual irreducible k-point count (NKPTS), NELECT, symmetry settings:
- Other relevant `INCAR` tags / path to complete final `INCAR`:
- Supported structure/chemistry/volume range and known limitations:

## Machine and execution layout

- Cluster, partition, CPU/GPU model, memory per node:
- Nodes, MPI ranks per node, threads per rank, requested memory:
- Modules / compiler / MPI / executable:
- `NCORE`, `KPAR` (or `NPAR` if used), actual `NBANDS`:
- Shared benchmark-table link (strongly recommended):
- Cold start or restart; filesystem / output settings affecting timings:

## Scan results

Copy rows as needed. Keep raw `cat-tools vasp check --verbose` output alongside this record. Use `NA` for unavailable measurements; never substitute zero. Differences below must be calculated against the named reference, not copied from CAT's absolute MaxF/MaxStress columns.

| Run / relative path | Axis and value / actual mesh | Reference | Delta E/atom (eV) | Max delta F (eV/angstrom) | Max delta stress (kbar) | SCF / ionic convergence | Decision |
| --- | --- | --- | --- | --- | --- | --- | --- |
| <run> | | | | | | | |

| Run | Job ID / state | SCF / ionic steps | Wall time (s) | Total CPU (s) | Allocated core-hours | MaxRSS/task (MiB) | Repeats / timing spread |
| --- | --- | --- | --- | --- | --- | --- | --- |
| <run> | | | | | | | |

Allocated core-hours = wall seconds × allocated CPU count / 3600; this differs from measured CPU time and may differ from the site's billing policy. MaxRSS/task is a task-level peak, not total node memory. Record total time and step counts so a change in iteration count is visible.

## Final validation and reuse

- Selected versus stricter settings on the same final geometry:
- Reference values and measured differences (including target property):
- Final relaxation force, stress, cell and restart checks:
- Coupled checks (k mesh / `SIGMA`, `ENCUT` / cell, `LREAL` / ROPT):
- Remaining property-specific convergence tests:
- Suggested first tests for the next user, and cases needing a fresh scan:
- Timing estimate applies to which system size, build and hardware:
