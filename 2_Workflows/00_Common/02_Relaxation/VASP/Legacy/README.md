# Original relaxation workflow

These are unchanged copies from the old Relaxation workflow: `INCAR.relax`, `INCAR.static`, `KPOINTS`, `job.sh`, and [README.original.md](README.original.md). Use the [maintained guide](../README.md) for new calculations.

To reproduce the old input setup in a separate working directory:

```bash
legacy="$PHONONTOOLS_ROOT/2_Workflows/00_Common/02_Relaxation/VASP/Legacy"
cp "$legacy/INCAR.relax" INCAR
cp "$legacy/KPOINTS" KPOINTS
cp "$legacy/job.sh" job.sh
```

Supply your own matching `POSCAR`/`POTCAR`. Review all inputs and adapt the job's modules, account, resources and email before `sbatch job.sh`. After convergence, use a new directory containing the final `CONTCAR` as `POSCAR` and the archived `INCAR.static` as `INCAR`, reviewing its settings against the relaxation.

The original guide references external `apgu` and `varc` helpers. Their maintained routes are `cat-tools vasp potcar` and `cat-tools vasp restart`; consult the new guide before using restart.