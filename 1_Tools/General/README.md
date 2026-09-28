## Slurm Log Scanner

**Use:** Inspect job state, diagnostics, and resource usage from Slurm logs without requiring the job to remain in the scheduler queue.

**Prereqs:** Python 3 with its standard library; no CAT installation is needed.

**Setup:** Define the helper used by the Bash aliases, adjusting the path:

```bash
slurm_scan() {
    python3 /path/to/PhononTools/1_Tools/General/SlurmScan.py "$@"
}
```

[SlurmScan.py](SlurmScan.py) is the complete standalone implementation. You can copy it elsewhere and update the helper path, or invoke it directly with Python.

**Usage:** `slurm_scan [OPTIONS]`

```bash
slurm_scan --root /path/to/calculations
slurm_scan --all --recursive --details
slurm_scan --performance
slurm_scan --failed restart_dirs.txt
slurm_scan array 123456 --root /path/to/logs
```

The highest numeric job ID in each calculation directory determines its latest run. Clean completed jobs are hidden by default. `--failed` replaces a text file listing non-completed directories; it does not resubmit jobs.

See [the Tools reference](../README.md#slurm-log-commands) and the [Bash helper README](../../0_bashrc.d/README.md) for filters, diagnostics, array scanning, and the `sscan`, `sall`, `sperf`, and `sfailed` aliases.