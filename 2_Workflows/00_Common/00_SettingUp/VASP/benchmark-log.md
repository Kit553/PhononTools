# Shared VASP Benchmark Log

To get a better understanding about the resources and timing of VASP calculations on the cluster, please use the benchmarking logs to record the assigned resources, task and timing of your calculations. Make sure to include queue timings as well but explicitly differentiate between time on queue and time on the node.

This is meant to make our VASP workflows more resource efficient and enable higher throughput as well as provide a list of systems that are already at a stable state.

Make sure to include your name in the log in case any set up differences impact the timing systematically.

The shared log book is atm placed on my scratch directory at `/scratch/peckert/shared/sharedEdit/benchmarkLog.md`, I encourage you to contribute your timings.

Add one row per calculation. **WAVECAR/CHGCAR means present as input at the start**, not written at the end; note ISTART/ICHARG if restart files were present but unused. Please also include a note if you already wrote out the WAVECAR/CHGCAR and can provide it. Use `—` for unavailable measurements. Keep failed runs and briefly note the failure (OOM, nodeFail, ...).

| System | Atoms | Task | Name | VASP / node type | ENCUT (eV) | KPOINTS | EDIFF | ALGO | ISMEAR / SIGMA | LREAL | NBANDS | NCORE | KPAR | NSIM | WAVECAR (Yes/No) | CHGCAR (Yes/No) | Nodes | CPUs allocated | MPI tasks (total) | CPUs/task | Mem requested/node (GiB) | Queue time (s) | Run time (s) | CPU efficiency (%) | MaxRSS/task (GiB) | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
