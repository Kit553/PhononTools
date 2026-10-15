# Common Workflows

**Use:** Establish structures, numerical settings and reference calculations for later workflows.

**Files:** [Setup inputs](00_SettingUp/VASP/) · [convergence inputs](01_Convergence/VASP/) · [relaxation inputs](02_Relaxation/VASP/).

Use this collection to establish the structure, numerical settings and reference outputs shared by the lattice-dynamics and electronic-structure workflows. The maintained calculation examples currently use VASP with Slurm and CAT.

| Stage | Procedure | Result to carry forward |
| --- | --- | --- |
| 00 | [Setting up](00_SettingUp/README.md) | Checked structure, matching potentials and working execution environment |
| 01 | [Convergence](01_Convergence/README.md) | Numerical settings and a record of their accuracy |
| 02 | [Relaxation](02_Relaxation/README.md) | Accepted geometry, constraints and restart checks |
| 03 | [Static reference](03_StaticReference/README.md) | Fixed-geometry energy, forces, stress and selected electronic outputs |
| 04 | [Dielectric and Born tensors](04_DielectricBorn/README.md) | Electronic dielectric tensor and Born charges for polar systems |
| 05 | [Ordering](05_Ordering/README.md) | Local configuration generation, convergence checks, cluster calculations and a mixing hull |

The numbering groups topics; ordering/disorder decisions belong before convergence and relaxation. Dielectric response is an optional branch selected by the intended property. Phonon calculations belong to the [harmonic collection](../10_LatticeDynamics/01_Harmonic/README.md).

## Working conventions

Make CAT available as described in [setup](00_SettingUp/VASP/README.md), or use an existing installation. Work in directories of your choice and copy only the needed templates; no project-path environment variable is required. See the [path and command conventions](../README.md#how-to-follow-a-workflow). Keep reference files in the repository distinct from run outputs. Example cutoffs, meshes, thresholds and Slurm resources require validation for your system. Preserve the chosen functional, spin/charge model and potential variants when adapting a template.

Use the settings from your [convergence record](01_Convergence/VASP/convergence-record.md) throughout the calculation sequence.

## Migration map

| Original material in `OLD_Workflows` | Maintained destination | Original files |
| --- | --- | --- |
| Root `INCAR` notes | [VASP setup](00_SettingUp/VASP/README.md) | [Original notes](../../OLD_Workflows/README.md) |
| `00_Convergence` | [VASP convergence](01_Convergence/VASP/README.md) | [Legacy](01_Convergence/VASP/Legacy/README.md) |
| `01_Relaxation` | [VASP relaxation](02_Relaxation/VASP/README.md) and [static reference](03_StaticReference/VASP/README.md) | [Legacy](02_Relaxation/VASP/Legacy/README.md) |
| `04_LO-TO` | [VASP dielectric/Born](04_DielectricBorn/VASP/README.md) | [Legacy](04_DielectricBorn/VASP/Legacy/README.md) |

`Legacy` retains the old inputs and documentation for provenance; its settings and commands are not the maintained procedure. The originals remain in `OLD_Workflows` while migration continues. Specialized workflows remain separate under [Workflows](../README.md).
