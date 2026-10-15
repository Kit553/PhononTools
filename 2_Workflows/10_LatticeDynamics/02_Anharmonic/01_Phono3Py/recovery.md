# Audit, recovery, logs and boundary-MFP collection

**Use:** Recover incomplete sampling and reuse validated intrinsic scattering for multiple boundary lengths without mixing scientific inputs.

**Files:** native YAML, force constants, correct `BORN`, grid inventory, GP outputs, and original/generated job scripts. Start in the relevant calculation directory.

## 1. Establish which dataset is valid

Record the intended structure/primitive basis, FC2/FC3, `BORN`, q mesh, temperature list, isotope settings, integration method and native version. Preserve inputs and job logs with the results. Recalculation after a wrong `BORN` must use the corrected response consistently: previously computed phonon caches, PP tensors, gamma, JDOS and Grüneisen may depend on it. Reuse only data whose dependencies match. CAT's structural/output audit cannot reconstruct the `BORN` identity of an old file that never recorded it. Treat missing provenance explicitly, not as proof of reuse.

Do not generate recovery jobs into directories still being written by active jobs. Once they finish, retain the original files and audit a clearly identified dataset. Recursive discovery is useful only when its root excludes incompatible runs and previous complete collection copies.

## 2. Audit JDOS temperatures, then prepare only required reruns

```bash
cat-tools phono3py audit jdos 21 21 21 --root jdos \
  --ts "1 100 300 600 1000" --num-freq-points 601 --recursive --no-collect
cat-tools phono3py audit jdos 21 21 21 --root jdos \
  --ts "1 100 300 600 1000" --num-freq-points 601 --recursive \
  --rerun 3 --no-collect --nac --job-dir jobs-recovery --tag corrected-jdos
```

Run from the parent directory for these examples. Supply the same primitive matrix, module stack, resources and extra native arguments as the desired calculation. Reports include missing/invalid GPs and temperature-specific rerun lists. Reruns are split by missing/bad temperature, preserving already-good temperatures. Exact row counts, finite three-column spectra and monotonic frequencies are checked. Duplicate files are reported; resolve their provenance instead of arbitrarily choosing one. Reports are replaced on a subsequent audit, so copy any needed record.

The legacy native multi-temperature JDOS writer can exit successfully while writing only its first temperature. This is why Slurm COMPLETED alone was insufficient. The maintained CAT sampler invokes native JDOS once per temperature.

## 3. Recollect RTA/Wigner for several lengths

After corrected gamma GP sampling is complete, run in `scattering/`:

```bash
for length in 0.1 0.5 1 10; do
  cat-tools phono3py audit rta 21 21 21 --ts "100 300 600" \
    --isotope --nac --bmfp "$length" --tag corrected-rta \
    --job-dir "jobs-collect-${length}um"
done
```

Substitute `wigner` only for the corresponding Wigner calculation/runtime. Forward identical temperatures, isotope/mass-variance, primitive and integration choices. Lengths are in µm. Each generated RTA/Wigner job reads audited GP gamma in a fresh `collection/TAG_bmfp_LENGTHum/` directory, copies text inputs and links the selected GP data. It refuses an existing destination. Submit the generated scripts after reviewing them; no submission occurs above without `--submit`.

Original stored boundary metadata need not equal the new length: intrinsic gamma is reusable. A new `BORN`/force dataset is a different case. Do not collect from a stale aggregate `kappa-m*.hdf5`, and do not sum complete results from different boundary lengths. Point the analysis YAML at each complete output explicitly. Grüneisen uses a separate full-mesh recalculation; JDOS recovery is independent.

## 4. Locate the actual scheduler log and timing

An array index is an index into the submission wrapper's task list, **not a GP**. The printed `Array task 242: ...job_..._0-2.sh` means task 242 executed that script; it does not imply GP 242. Inspect the submitted array wrapper and its frozen task list, together with the job's `StdOut`, `StdErr`, work directory and command.

```bash
scontrol show job ARRAYJOBID_TASKID
sacct -j ARRAYJOBID --format=JobIDRaw,JobName,State,ExitCode,Elapsed,Timelimit,MaxRSS
```

Replace uppercase IDs. `scontrol` may no longer retain finished jobs; use accounting and archived submission records. Inspect accounting step rows for memory where the parent record is blank. For future arrays use unique `%A_%a` output/error paths and create their parent log directory before submission. Freeze the task list; do not rebuild it by globbing a directory while an array is running.

Native output is mirrored to Slurm and appended to CAT's calculation log. Elapsed walltime is not aggregate CPU time. A task that saved only one JDOS temperature cannot establish the time for CAT's complete per-temperature sequence; benchmark representative GPs and requested temperatures when sizing jobs.

The current CAT sampler generates individual batch scripts. If you wrap them in an external array, that wrapper owns array/log mapping and must be retained.
