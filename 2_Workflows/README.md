# Workflows

**Use:** Choose and follow a calculation workflow using your own inputs and working directory.

**Files:** [Preparation tools](../1_Tools/README.md) · [analysis commands](../3_Analysis/README.md) · templates linked in each workflow.

Workflows connect the reusable preparation commands in `1_Tools` with calculations and the analysis commands in `3_Analysis`.

| Area | Procedure |
| --- | --- |
| [00_Common](00_Common/README.md) | Setup, convergence, relaxation, static references and dielectric/Born response |
| [10_LatticeDynamics](10_LatticeDynamics/README.md) | Harmonic/anharmonic phonons, QHA, FreqBall, mode viewing and frozen-mode PES |
| [20_ElectronicStructure](20_ElectronicStructure/README.md) | VASP electronic bands/DOS, LOBSTER preparation, bonding and orbital fatbands |
| [99_Misc](99_Misc/README.md) | SoftBV/BVPA to directional phonon analysis |

## How to follow a workflow

Each guide opens with a short **Use:** statement and a **Files:** line linking the relevant templates or naming the inputs you supply. The numbered steps follow below.

Start in a calculation or analysis directory of your choice with the required tools available. Each procedure states its inputs and where to run its commands. No particular home/scratch layout or exported project variable is required.

- Paths such as `/path/to/accepted/POSCAR` and uppercase tokens such as `ELEMENT` are placeholders: replace them before running the command. Quoted paths can contain spaces. Angle brackets should be avoided because shells interpret them as redirection.
- Short names such as `relax_01`, `pool-r01`, and `settings.yaml` are example names. Keep your own naming scheme and update subsequent references consistently. Names required by the program, such as VASP's `INCAR`, are identified as such.
- Copy only the linked templates you need into your working directory and edit them with your preferred editor. Material examples live in the tools examples; their species, model settings and numerical parameters are not general recommended defaults.
- Choose local execution or your scheduler. `sbatch` examples apply to Slurm; module names, resources and launch commands must match your installation. A preparation command does not require a cluster unless explicitly stated.
- Installation instructions describe one supported environment. An existing compatible installation works too. Environment variables used by a program (for example thread limits or a potential-library setting) are distinct from optional shortcuts for project paths.

Historical files under `Legacy` retain their original assumptions for provenance; use the maintained procedure linked beside them for current work.

## Coverage of the previous workflows

The [coverage matrix](COVERAGE.md) maps every `OLD_Workflows` topic to a maintained procedure, identifies corrections to old commands, and distinguishes CAT features from upstream/native output handling. Reproducing the scientific task is the minimum scope; retaining an obsolete filename or output layout is not required. Existing legacy data and scripts remain untouched for provenance.
