# VASP Static Reference

**Use:** Calculate reference energies, forces and stress on an accepted geometry.

**Files:** [INCAR.static](INCAR.static) · [submit.sh](submit.sh).

**Prerequisites:** Complete [convergence](../../01_Convergence/VASP/README.md) and [relaxation](../../02_Relaxation/VASP/README.md) first.

The `INCAR` is a nonmagnetic bulk PBEsol starting example, with the same baseline as the relaxation template.

## 1. Prepare the static inputs

Start in a directory of your choice. Replace `/path/to/accepted-relaxation` with the run whose geometry you accepted; `static_01` is an example output name:

```bash
mkdir static_01
cp /path/to/accepted-relaxation/CONTCAR static_01/POSCAR
cp /path/to/accepted-relaxation/{INCAR,POTCAR,KPOINTS,submit.sh} static_01/
```

Edit the copied `INCAR`, retaining the converged electronic settings:

```text
SYSTEM = static-reference
IBRION = -1                 # Disable relaxation
NSW = 0                     # Set the number of ionic steps to zero
ISIF = 2
```

Remove relaxation-only `EDIFFG` and `IWAVPR`. Review `ISTART` and `ICHARG`: the provided template starts from scratch and the commands above do not copy electronic restart files. Keep the accepted geometry unchanged during this run. See VASP's [IBRION reference](https://vasp.at/wiki/IBRION).

If you do not already have a reviewed scheduler script, copy [submit.sh](submit.sh) and adapt its PALMA module stack, account, resources and time limit. Set output requirements before submission:

| Output | Decision |
| --- | --- |
| `OUTCAR`, `vasprun.xml`, `OSZICAR` | Retain for convergence, energies, forces |
| `WAVECAR` | Set `LWAVE=.TRUE.` if the downstream workflow requires wavefunctions, can be helpful for timings |
| `CHGCAR` | Set `LCHARG=.TRUE.` if charge analysis or a compatible restart needs it, rarely needed but cheap |

The template disables both large optional files. Specialized workflows can require additional settings or their own SCF run; consult them before assuming this reference supplies every input.

## 2. Execute the static calculation

Use your local launcher or scheduler from `static_01`. For Slurm:

```bash
cd static_01
sbatch submit.sh
```

## 3. Check the result

After completion, from the parent directory:

```bash
cat-tools vasp check --prefix static_ --verbose
```

Check normal termination and electronic convergence, then compare final forces and stress against the relaxation's acceptance criteria. If you changed sampling or occupations and forces change materially, return to relaxation with appropriate settings. A scheduler `COMPLETED` state alone does not establish convergence!

Carry this accepted reference into [lattice dynamics](../../../10_LatticeDynamics/README.md), [electronic structure](../../../20_ElectronicStructure/README.md), or [hull analysis](../../05_Ordering/README.md) as needed.
