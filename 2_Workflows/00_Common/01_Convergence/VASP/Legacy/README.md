# Original convergence workflow

These files are unchanged copies from the old convergence workflow. They are retained for reproducing older runs; use the [maintained CAT workflow](../README.md) for new work. [README.original.md](README.original.md) preserves the original discussion, including its numerical advice.

From a separate calculation directory containing your checked `POSCAR`, `POTCAR`, `INCAR`, and site-adapted `submit.sh`:

```bash
legacy="$PHONONTOOLS_ROOT/2_Workflows/00_Common/01_Convergence/VASP/Legacy"
bash "$legacy/KPOINTS.sh" POSCAR 50
bash "$legacy/ConvGen.sh" --encut --range 1.0 2.5 0.25 --no-submit
bash "$legacy/ConvCheck.sh" --prefix ENCUT_
```

The checker becomes useful once the generated jobs have run. The generator submits through `sbatch` by default: retain `--no-submit` for inspection. Submit generated directories using the procedure in the maintained guide; rerunning the generator over existing directories fails. Bash, standard Unix utilities, and Slurm are needed for generation/submission; structural comparisons with `ConvCheck.sh --relax` also need Python 3.

Use `bash "$legacy/ConvGen.sh" --help` and `bash "$legacy/ConvCheck.sh" --help` for the original options.
