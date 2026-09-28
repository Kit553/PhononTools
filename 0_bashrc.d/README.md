## Bash Helpers for Slurm Jobs

**Use:** Monitor slurm jobs, navigate to calculation directories, inspect job output and handle active jobs with out-of-memory errors directly from the terminal.

**Prerequisites:**

- Bash 4 or newer; the OOM helper uses associative arrays.
- Access to the cluster's slurm commands: `squeue`, `sacct`, `sstat`, `scontrol`, `scancel` and `sbatch`. Job history and historical working directories depend on the available accounting records. This is standard on Palma-II
- Standard Linux command-line tools, including `awk`, `sed`, `grep`, `find` and `realpath`.
- Python 3 for the optional [SlurmScan.py](../1_Tools/General/SlurmScan.py) aliases.

**Folder Setup:**

| File | Command | Purpose |
|:-----|:--------|:--------|
| [squeueClean.sh](squeueClean.sh) | `sq` | Display your queue with working directories and running/pending job counts. |
| [jobQueued.sh](jobQueued.sh) | `jq` | Check whether your current directory has running or pending jobs to prevent double starting. |
| [changeJobDir.sh](changeJobDir.sh) | `cdj` | Change to a job's working directory using its job ID. |
| [recentJobs.sh](recentJobs.sh) | `recentjobs` | Show accounting records for finished jobs in a selected time window. |
| [jobperf.sh](jobperf.sh) | `jp` | Inspect CPU and memory statistics and allocated resources for a running job. |
| [oomKill.sh](oomKill.sh) | `ohmu` | Find, cancel or resubmit active jobs whose output contains memory errors. |
| [slurm.sh](slurm.sh) | `sscan`, `sall`, `sperf`, `sfailed` | Shortcuts for the slurm output scanner. |


---

## Setup

Copy all shell scripts from this folder into `~/.bashrc.d`, replacing `/path/to/PhononTools` with your repository location:

```bash
mkdir -p ~/.bashrc.d
cp /path/to/PhononTools/0_bashrc.d/*.sh ~/.bashrc.d/
```

Make sure your `~/.bashrc` sources the scripts in that directory. If it does not already do so, add this block:

```bash
for helper in "$HOME"/.bashrc.d/*.sh; do
    [ -f "$helper" ] && source "$helper"
done
unset helper
```

Reload the configuration to make the commands available in your current terminal:

```bash
source ~/.bashrc
```

The scanner aliases require the additional `slurm_scan` setup described under `slurm.sh` below, make sure to source again after completion.

---

## Shell Commands

### `squeueClean.sh` — `sq`

**Purpose:** Display your slurm queue with working directories and a count of running and pending jobs in a human readable format.

The script removes `/scratch/tmp/$USER/$` from displayed working directories. Adapt that prefix in the script to your cluster's scratch directory if desired; this only changes the display. The separator width uses the interactive shell's `COLUMNS` variable.

**Usage:**

```bash
sq
```

`sq` takes no options. It lists your jobs with their IDs, names, states, elapsed times, time limits, submission times, nodes or pending reasons, expected start times, priorities and working directories. A final line reports the number of running and pending jobs.

---

### `jobQueued.sh` — `jq`

**Purpose:** Check whether you already have a running or pending job in the current calculation directory to prevent double starting your jobs if unsure.

Some versions of JSON processors take this command as their version of `sed`, make sure to either change the alias in case you use both.

**Usage:**

```bash
jq
```

`jq` takes no options. It compares your running and pending jobs' working directories with the resolved absolute path of the current directory. Subdirectories are not included.

- A match prints `YES`, followed by the matching jobs' IDs, names, states, elapsed times, time limits, submission times, nodes or pending reasons, start times, priorities and working directories. Exit status is `0`.
- No match prints `NO` and returns exit status `1`.
- This allows it to be used in automatic job starters while retaining sensible exit scenarios.

---

### `changeJobDir.sh` — `cdj`

**Purpose:** Change directly to a job's working directory using its slurm job ID.

**Usage:** `cdj JOB_ID`

```bash
cdj 123456
```

`JOB_ID` is required. The function first queries `scontrol`, then falls back to `sacct` if no working directory is found. Historical jobs therefore depend on available accounting records.

On success, it changes directory and prints the job ID and destination. A missing argument prints usage and returns `2`. An unavailable working directory or a directory that no longer exists returns `1`; a failed `cd` also returns a nonzero status. Does work on jobs no longer in the queue as long as you have their ID.

---

### `recentJobs.sh` — `recentjobs`

**Purpose:** Show accounting records for finished jobs within a selected time window.

The script removes `/scratch/tmp/$USER$/` from displayed working directories. Adapt this prefix to your own scratch directory if desired.

**Usage:** `recentjobs [DURATION]`

```bash
recentjobs          # Default time 1h
recentjobs 30min
recentjobs 24h
recentjobs 2d
```

The duration is a non-negative number immediately followed by one of these suffixes:

| Unit | Accepted suffixes |
|:-----|:------------------|
| Minutes | `m`, `min`, `mins`, `minute`, `minutes` |
| Hours | `h`, `hr`, `hrs`, `hour`, `hours` |
| Days | `d`, `day`, `days` |

The default is `1h`. Invalid durations print usage and return `1`.

The report contains job ID, state, runtime and working directory for records returned by `sacct` between `now - DURATION` and `now`. It includes `COMPLETED`, `FAILED`, `CANCELLED`, `TIMEOUT`, `OUT_OF_MEMORY`, `NODE_FAIL` and `PREEMPTED`, and excludes job steps. Cancellation states are displayed as `CANCELLED by user` without preserving the cancelling user's ID. Time zero works for whatever shenanigans you can think of.

---

### `jobperf.sh` — `jp`

**Purpose:** Inspect a running job's CPU and memory usage alongside its allocated resources. Useful to check for hanging jobs and running performance indicators.

**Usage:** `jp JOBID`

```bash
jp 123456
```

`JOBID` is required. A missing argument prints usage and returns `1`. The shared setup above loads the `jp` function from `jobperf.sh`.

The function queries `sstat` for `JobID`, `AveCPU`, `AveRSS`, `MaxRSS`, `AveVMSize` and `MaxVMSize`: the job step ID, average CPU time, average and maximum resident memory, and average and maximum virtual memory. Statistics depend on the cluster's job accounting configuration and the running steps selected by `sstat`.

It then prints the job's `AllocTRES` field from `scontrol show job` as `Resources: ...`, followed by the current date and time. Each invocation produces one snapshot; rerun it for an updated view. For finished jobs, use accounting records or the `sperf` scanner alias instead. Tested as part of automatic workflows.

---

### `oomKill.sh` — `ohmu`

**Purpose:** Find active jobs whose slurm output contains memory errors, with options to cancel or resubmit them. This is especially useful as VASP tends to throw oom error but then not cancel the job and return error.

**Usage:** `ohmu [OPTION] [DIR]`

| Option | Action |
|:-------|:-------|
| No option | Report matching active jobs. |
| `-c`, `--cancel` | Cancel matching active jobs. |
| `-r`, `--restart` | Cancel matching active jobs and submit their scripts again from their respective directory. |
| `-h`, `--help` | Print built-in help and return. |

`DIR` defaults to the current directory. Supply at most one option, before the directory argument.

```bash
ohmu                               # Report matches below the current directory
ohmu /path/to/calculations         # Report matches below a selected directory
ohmu -c /path/to/calculations      # Cancel matching active jobs
ohmu --restart /path/to/calculations
ohmu --help
```

The helper searches recursively for `slurm*.out` files containing case-insensitive memory-error patterns: `oomkill` (also with an underscore, space or hyphen), `out of memory`, `exceeded.*memory`, `memory limit` or `Killed process`. It extracts job IDs only from filenames shaped like `slurm-123456.out` or `slurm-123456_7.out` and processes each job ID once.

Only jobs still reported by `squeue` as `RUNNING` or `PENDING` are included. Finished jobs are excluded even if their logs contain memory errors. The report shows `JOBID`, `STATE`, `ACTION` and `DIRECTORY`.

**Cancellation and restart behavior:**

Both action options act immediately on all matches in the selected directory tree. Restart mode retrieves the original script name from `scontrol`, requires a file with that basename beside the Slurm output, and runs `sbatch ./SCRIPT` from that directory after cancellation succeeds. A missing script leaves the job running; a failed submission after cancellation leaves it cancelled.

Restart mode uses the script's current settings. It does not increase memory requests, restore options passed to the original `sbatch` command, or reconstruct a single array task submission. Adjust the submission script before restarting a job that needs different resources. A common error this solves is a too large `KPOINTS` setting that requires changing the input but not the job itself.

| Action shown | Meaning |
|:-------------|:--------|
| `-` | Report only; no action taken. |
| `CANCELLED` | Cancellation succeeded in cancel mode. |
| `FAILED` | Cancellation failed in cancel mode. |
| `NO SCRIPT` | The script name could not be retrieved, or the script was missing beside the output. |
| `CANCEL FAIL` | Cancellation failed in restart mode; no submission attempted. |
| `RESTARTED` | Cancellation and submission succeeded; submission output follows. |
| `SUBMIT FAIL` | Cancellation succeeded, but submission failed; error output follows. |

---

### `slurm.sh` — slurm output scanner aliases

**Purpose:** Provide shortcuts for [SlurmScan.py](../1_Tools/General/SlurmScan.py), which summarizes job states, diagnostics and resource usage from slurm output files.

**Setup:**

The aliases require a command named `slurm_scan`. If it is not already configured, add this function to `~/.bashrc.d/slurm.sh`, replacing the repository path with your own:

```bash
slurm_scan() {
    python3 /path/to/PhononTools/1_Tools/General/SlurmScan.py "$@"
}
```

Reload with `source ~/.bashrc`. The shared setup above loads both the function and aliases from `~/.bashrc.d/slurm.sh`.

**Usage:**

| Alias | Expansion | Behavior |
|:------|:----------|:---------|
| `sscan` | `slurm_scan` | Summarize latest runs, omitting clean completed jobs. |
| `sall` | `slurm_scan --all` | Include clean completed jobs. |
| `sperf` | `slurm_scan --performance` | Report resource/performance data for completed latest runs. |
| `sfailed` | `slurm_scan --failed` | Write directories whose latest run is not classified as completed. |

Run from the parent of your calculation folders, or select another directory with `--root`. By default, the scanner reads immediate subdirectories and treats the highest numeric job ID in each directory as the current run. Status and diagnostics come from log contents.

```bash
sscan
sall
sperf
sfailed                     # Write failed_subdirs.txt under the scan root
sfailed restart_dirs.txt    # Choose a different output filename
sscan --recursive --details
sscan --history
sscan --state FAILED,TIMEOUT
sscan --diag OOM --paths
sscan --root /path/to/calculations
```

`sfailed` can include unfinished jobs. It replaces the selected output file and does not resubmit calculations. Relative output filenames are resolved against the scan root.

**Scanner options:**

| Option | Purpose |
|:-------|:--------|
| `-a`, `--all` | Include clean completed latest runs. |
| `-p`, `--performance` | Show resource/performance information for completed latest runs. |
| `--history` | Also show older runs. |
| `--details` | Show representative diagnostic lines. |
| `-r`, `--recursive` | Search all descendant directories. |
| `--failed [FILE]`, `--failed-file [FILE]` | Write non-completed latest-run directories; default filename is `failed_subdirs.txt`. |
| `--root DIR` | Set the scan root; default is the current directory. |
| `--dir-regex REGEX` | Filter parent directory names with a regular expression; default is `.*`. |
| `--state STATE[,STATE...]` | Filter by Slurm state. |
| `--diag TYPE[,TYPE...]` | Filter by diagnostic type. |
| `--paths` | Print matching directory paths only; cannot be combined with performance mode. |
| `-h`, `--help` | Show command help. |

Supported diagnostic types are `OOM`, `TIME_LIMIT`, `HFI/PSM2`, `EAGER_MSG`, `MPI_ABORT`, `SEGFAULT`, `KILLED`, `CANCELLED`, `NODE_FAIL`, `VASP_ERROR` and `ERROR`.

For array output, use `slurm_scan array ARRAY_ID [ARRAY_ID ...]`:

```bash
slurm_scan array 123456
slurm_scan array 123456 123457 --performance
slurm_scan array 123456 --root /path/to/logs --details
slurm_scan array --help
```

Array mode scans `slurm-ARRAYID_TASKID.out` files and accepts `--performance`, `--details`, `--root`, `--state`, `--diag` and `--paths`. Its root is the directory containing the array output files. The other directory-scan options listed above do not apply to array mode.

Detected job failures do not themselves cause a non-zero scanner exit status; invocation, discovery or parsing errors do.
