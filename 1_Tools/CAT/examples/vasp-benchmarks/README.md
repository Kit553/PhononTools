# Optional Si/Al benchmark examples

These concrete Si/Al and PALMA examples illustrate benchmark construction. For a procedure using your own structures and hardware, follow the [generic benchmark workflow](../../../../2_Workflows/00_Common/00_SettingUp/VASP/benchmarking.md). Run them on allocated compute nodes with licensed potentials and CAT's `structure` extra installed.

Start with the small case; prepare larger trials only when you need them and run into problems with your usual resource allocations

**Strong recommendation:** The shared `benchmarkLog.md` is currently placed at `scratch/tmp/peckert/shared/sharedEdit/`. It records benchmarking test that may be able to help you decide the resource allocations needed for your jobs.  Consult existing entries before testing, record unsuccessful trials as well as improvements, and link the shared table from your production workflow. Personal notes alone make it difficult for colleagues to compare settings or build on your measurements. See [the shared logging procedure](#5-share-results-and-improve-settings-together).

## 1. Understand the PALMA recommendation

The supplied HPC guidance reports good experience with `NCORE=18`, `NSIM=18` on 36-physical-core normal-queue nodes for **more than 250 metal atoms**, using `palma/2019a`, `intel/2019a`, `VASP/5.4.4-vtst`. Treat that combination as a benchmark candidate for comparable work. Record the module stack actually used (!The newer VASP 6.2.0 version is placed on a different stack!).

| Choice | How this collection uses it |
| --- | --- |
| Exclusive nodes | PALMA launchers use `--exclusive`, following the supplied site guidance. This reserves whole nodes even with fewer MPI ranks. |
| Physical cores | CPU/MPI launchers use one thread per rank and `--hint=nomultithread`. The documented full-node trial has 36 ranks; check the actual node type. |
| Stack | `ulimit -s unlimited` before launching VASP. |
| Requeue | Keep `--no-requeue`, including for benchmark runs, so an automatic restart does not overwrite the run being measured. |
| Memory and time | Request measured needs plus headroom, memory can be especially annoying with larg(er) systems. |
| `NCORE` | Start small and test divisors of ranks per k-point group. Keep `NPAR` absent. |
| `NSIM` | Benchmark blocking separately after rank/`NCORE` tuning; 18 is an optional large-system trial. |
| `MAXMEM` | Selected-routine memory hint, not an overall cap or a scheduler memory request. Omitted from the ordinary SCF/relaxation/DFPT baselines. |

PALMA-II's `5139` estimate is approximately `185000 / 36` MB per task. Use the actual allocation, ranks and units if a method that reads `MAXMEM` needs it, allowing additional overhead. VASP documents this tag mainly for selected memory-intensive routines such as GW/RPA. It cannot make an oversized SCF fit into a small slurm allocation. See [MAXMEM](https://vasp.at/wiki/MAXMEM).

An exclusive allocation does not automatically grant all node memory! High memory requests can restrict eligible nodes, check node availability at the [HPC confluence](https://palma.uni-muenster.de/documentation/hardware/partitions/); a short wall time can help slurm's backfilling. For tiny workloads, compare fewer ranks within a single exclusive allocation, or use a site-approved smaller node/partition. Fewer active ranks do not imply lower whole-node charges. See [Slurm allocation options](https://slurm.schedmd.com/sbatch.html).

## 2. Create three explicit example structures

From a fresh working directory, with CAT's Python environment active:

```bash
mkdir vasp-benchmarks
cd vasp-benchmarks
python - <<'PY'

from pathlib import Path
from ase.build import bulk
from ase.io import write

cases = {
    "small-si": bulk("Si", "diamond", a=5.43),
    "medium-si": bulk("Si", "diamond", a=5.43, cubic=True).repeat((2, 2, 2)),
    "large-al": bulk("Al", "fcc", a=4.05, cubic=True).repeat((4, 4, 4)),
}
for name, atoms in cases.items():
    path = Path(name) / "seed"
    path.mkdir(parents=True, exist_ok=False)
    write(path / "POSCAR", atoms, format="vasp", direct=True, vasp5=True)
    print(name, len(atoms), "atoms")
PY
```

This gives a 2-atom primitive Si cell, a 64-atom conventional-cell Si supercell, and a 256-atom Al supercell. The lattice constants are explicit example values, not relaxed results. ASE's [bulk builder](https://docs.ase-lib.org/ase/build/build.html) defines the primitive/conventional choices.

| Case | Initial Gamma-centered electronic mesh | First allocation | Initial `KPAR`, `NCORE` |
| --- | --- | --- | --- |
| Small Si, 2 atoms | `6 6 6` | 1 node, 4 ranks, 8G, 15 min | 1, 1 |
| Medium Si, 64 atoms | `2 2 2` | 1 node, 18 ranks, 16G, 30 min | 1, 3 |
| Large Al, 256 atoms | `2 2 2` | 1 node, 36 ranks, 32G, 30 min | 1, 6 |

For each `seed` directory, assemble `POTCAR` with `cat-tools vasp potcar --sources --verbose`. Keep the same potential release/variant throughout a case. Create an explicit `KPOINTS` file, replacing the mesh line from the table:

```text
Benchmark mesh
0
Gamma
6 6 6
0 0 0
```

Copy `INCAR.example` and `job.sh.example` from `2_Workflows/00_Common/00_SettingUp/VASP` in the repository to the seed as `INCAR` and `submit.sh`. Set the following electronic/ionic choices explicitly:

```text
ISTART = 0
ICHARG = 2
IBRION = -1
NSW = 0
ISIF = 2
ALGO = Normal
LREAL = .FALSE.
EDIFF = 1E-6
NELM = 100
NSIM = 4
LWAVE = .FALSE.
LCHARG = .FALSE.
```

Keep the chosen functional and a fixed explicit `ENCUT` suitable for the selected potentials; inspect their `ENMAX` values. For this timing exercise use `ISMEAR=0`, `SIGMA=0.05` for Si and `ISMEAR=0`, `SIGMA=0.10` for Al. These are starting occupations for the exercise, not production convergence advice. Set `KPAR`/`NCORE` and scheduler resources from the table, and review the modules. Do not submit from the seed directory.

## 3. Change one performance setting at a time

Start with a pilot copied from the small seed, from `vasp-benchmarks`:

```bash
mkdir small-si/r4-nc1-repeat1
cp small-si/seed/{POSCAR,POTCAR,INCAR,KPOINTS,submit.sh} small-si/r4-nc1-repeat1/
cd small-si/r4-nc1-repeat1
sha256sum POSCAR POTCAR INCAR KPOINTS submit.sh > input.sha256
jobid=$(sbatch --parsable submit.sh) && printf '%s\n' "${jobid%%;*}" > jobid.txt
```

Suggested **small sets of trials**, not an all-combinations submission sweep:

| Case | First comparison, `KPAR=1`, `NSIM=4` | Follow-up |
| --- | --- | --- |
| Small Si | 4 ranks: `NCORE` 1 versus 2; optionally compare 1 or 2 ranks | Stop increasing ranks if startup/communication dominates. |
| Medium Si | 18 ranks: `NCORE` 3 versus 6; then compare 36 ranks using the better setting | Test `KPAR` only after inspecting `IBZKPT`. |
| Large Al | 36 ranks: `NCORE` 6, 9, 18 | At the best `NCORE`, compare `NSIM` 4, 8, 18. |

`NSIM` is a band-blocking parameter, especially relevant to RMM-DIIS, with some use in Davidson/force calculations; its payoff depends on the algorithm. Keep `ALGO=Normal` fixed for this exercise rather than changing `ALGO` to chase the `NSIM` recommendation. See [NSIM](https://vasp.at/wiki/NSIM).

For CPU/MPI trials require `R % KPAR == 0` and `(R / KPAR) % NCORE == 0`. Choose `KPAR` compatible with the actual irreducible *k*-point count. Keep it at 1 for one *k* point. One node per irreducible *k* point is a possible layout to test for sufficiently large work, not a node-request formula for small cells. See [NCORE](https://vasp.at/wiki/NCORE).

Inspect actual `NBANDS` in every `out_VASP`: VASP can round it for the parallel layout. If it differs, use a common adequately large band count compatible with all selected layouts and repeat the comparison. Do not go lower than the highest occupied band.

The optional large-metal trial can use `NCORE=18`, `NSIM=18`, 36 ranks and `KPAR=1`. Only if the measured memory requirement warrants the site's high-memory example, change the script to `#SBATCH --mem=185G`; keep the wall time tied to your measurement. Do not add `MAXMEM=5139` as a supposed global memory safeguard. Closely monitor how the high `NCORE` affects memory.

## 4. Measure and choose

For a completed trial (replace the job ID):

```bash
grep -E 'LOOP:|LOOP\+:|Elapsed time|Total CPU time|NBANDS' OUTCAR
sacct -j JOB_ID --format=JobID,State,Submit,Start,Elapsed,AllocCPUS,TotalCPU,MaxRSS,ReqMem
```

Retain job-step rows: MaxRSS is generally a task-level peak, not the sum of node memory. Use available site node-memory metrics and headroom for sizing; a single rank's `OUTCAR` memory estimate is also not a total. See [Slurm accounting fields](https://slurm.schedmd.com/sacct.html).

Record total elapsed time, median steady electronic-iteration time (excluding startup), SCF iteration count, actual `NBANDS`, energy/force agreement and memory. Repeat the promising configurations three times and compare medians and spread. For extremely short runs, serial repeats in separate directories within one exclusive allocation avoid repeated queue/startup overhead; never reuse a completed directory without accounting for its restart files!

Enter the input parameters, resources and measured results in the [shared benchmark table](../../../../2_Workflows/00_Common/00_SettingUp/VASP/benchmark-log.md), one row per calculation.

For exclusive jobs, use **allocated/charged cores**, not only active ranks, in `core-hours = allocated cores × elapsed seconds / 3600`. Record submission-to-start delay separately from runtime; it depends on queue conditions. Choose a layout that balances throughput, cost and waiting time, rather than maximizing reserved cores.

Finally, test the chosen layout on the actual production task: several representative ionic steps for relaxation, or a response calculation for dielectric/DFPT work. Then carry the measured resources into the workflow's submission script and freeze them during numerical convergence scans. See VASP's [parallelization guide](https://vasp.at/wiki/Optimizing_the_parallelization).

## 5. Share results and improve settings together

**Please contribute:** Add one row per calculation to the shared [benchmark table](../../../../2_Workflows/00_Common/00_SettingUp/VASP/benchmark-log.md), currently at `/scratch/peckert/shared/sharedEdit/benchmark.md`. Keep it accessible to everyone in the group so colleagues can compare settings and contribute their own timings or share restart files.

Record the system and input parameters, allocated resources, queue time, run time, CPU efficiency and peak memory. Mark whether `WAVECAR` and `CHGCAR` were present at startup with **Yes/No**. Include your name so colleagues can ask about the setup. Use the notes column for failures or relevant differences such as restart settings.
