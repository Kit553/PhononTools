# Benchmark VASP for your workload

**Use:** Develop an idea about the resource usage of VASP to raise efficiency.

**Files:** [Benchmark log](benchmark-log.md) · [INCAR.example](INCAR.example) · [job.sh.example](job.sh.example).

## 1. Choose a representative calculation

Use your own structure and physical model, with a provisional numerical setup appropriate to the intended production calculation. Benchmarking chooses execution resources; [convergence testing](../../01_Convergence/VASP/README.md) establishes scientific accuracy. Repeat relevant performance checks if convergence changes the workload substantially.

Start in a directory of your choice. Keep a seed input set and separate trial outputs; the example below calls them `seed` and `trial-01`. Other names and existing layouts work too. Local execution and scheduled execution are both suitable as long as you can record the hardware and resources actually used.

For a study spanning sizes or methods, choose a small pilot and a representative expensive case. Optional [Si/Al case studies](../../../../1_Tools/CAT/examples/vasp-benchmarks/README.md) and [site notes](../../../../1_Tools/CAT/vasp-site-notes.md) can be found in tools.

## 2. Prepare a fixed input set

Put your `POSCAR`, matching `POTCAR`, explicit `KPOINTS` and reviewed `INCAR` in the seed directory; if using a scheduler, include your job script. Use the same geometry, potential release, physical model, cutoff, mesh and electronic convergence criteria across trials. For an SCF timing comparison, use static inputs (`IBRION=-1`, `NSW=0`).

Decide whether each trial starts from scratch or uses identical compatible restart files, and record that decision. Do not compare one cold start with another warm start without accounting for the difference. Retain the input set unchanged.

From the directory containing the example seed:

```bash
mkdir trial-01
cp seed/{POSCAR,POTCAR,INCAR,KPOINTS} trial-01/
# If using a scheduler:
cp seed/submit.sh trial-01/
```

Edit only the performance setting or resource allocation being tested. Choose MPI ranks/threads for the actual hardware and VASP build. The [parallelization discussion](README.md#7-match-parallelization-to-the-node-layout) explains the relevant work groups; input templates are starting examples rather than an allocation prescription.

## 3. Execute and record one pilot

Inside `trial-01`, record input identities and use your tested launcher:

```bash
cd trial-01
sha256sum POSCAR POTCAR INCAR KPOINTS > input.sha256
```

For a reviewed Slurm job script:

```bash
sbatch --parsable submit.sh > jobid.txt
```

For local execution, use the command appropriate to your build; for example, `mpirun -np RANKS vasp_std > vasp.out 2> vasp.err`, replacing `RANKS` with the chosen rank count. Do not run both launch alternatives for the same trial.

Verify normal termination and electronic convergence before treating timing as a valid result. Reaching an iteration limit quickly is not improved performance.

## 4. Compare a small set of resource choices

Return to the directory containing `seed` and `trial-01`. Copy the same seed into new trial directories, then change one factor at a time:

- Compare a few MPI rank counts supported by the machine/allocation.
- For CPU/MPI runs, test compatible `NCORE` values with `NPAR` absent; keep `KPAR` fixed initially, then test it where the irreducible k-point count permits.
- Test threading or `NSIM` separately when relevant to the chosen build/algorithm.
- Request memory and wall time from measured use plus reasonable headroom.

For a standard single-image CPU/MPI layout, require divisibility of total ranks by `KPAR` and of ranks per k-point group by `NCORE`. Keep `KPAR=1` for a single k point. Do not transfer those choices automatically to GPU, hybrid or response workloads. See the [VASP parallelization guide](https://vasp.at/wiki/Optimizing_the_parallelization).

Inspect actual `NBANDS`: VASP can adjust it for the parallel layout. If trials end up with different band counts, repeat using a common sufficiently large compatible count before attributing timing changes to parallelization alone.

## 5. Measure, repeat and choose

Inside each completed trial, inspect the output:

```bash
rg 'LOOP:|LOOP\+:|Elapsed time|Total CPU time|NBANDS' OUTCAR
```

If Slurm accounting is available, replace `JOB_ID` with the recorded ID:

```bash
sacct -j JOB_ID --format=JobID,State,Submit,Start,Elapsed,AllocCPUS,TotalCPU,MaxRSS,ReqMem
```

Record elapsed time, steady electronic-iteration times, SCF iteration count, actual `NBANDS`, memory and energy/force agreement. MaxRSS is generally a task-level peak, not total node memory; retain job-step information and use available node metrics for sizing. See [Slurm accounting fields](https://slurm.schedmd.com/sacct.html).

Repeat promising choices in fresh directories and compare their timing spread. Use the [log template](benchmark-log.md) locally or in an agreed shared location. Separate queue time from execution time. For charged core-hours, use allocated cores rather than just active ranks when the allocation reserves more resources.

Choose a layout that meets memory needs and gives useful throughput at acceptable cost. Verify it on the actual production task, then keep it fixed during numerical convergence comparisons. Preserve failed trials with their failure reason.
