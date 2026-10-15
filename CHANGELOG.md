# CAT changelog

## Unreleased

- Completed maintained anharmonic, quasi-harmonic, electronic-structure and mode-analysis procedures, with a task-by-task legacy coverage inventory. Added version-specific Phono3Py preparation/recovery, native Grüneisen/QHA, basis-derived LOBSTER setup, electronic bands/DOS and mode-viewing examples.
- Added a modern 1D/2D frozen-mode PES helper with explicit mass normalization, commensurability checks, structure-verified static energy collection and plots. Included workflow helpers and input templates in source distributions.
- Consolidated FreqBall calculation, family preparation, mass scans, ball viewing and static plotting under `cat-freqball`; retained existing entry points and positional YAML syntax.
- Added shared named element/site/direction projections for spectral tensors, descriptors and mode fatbands, with explicit Cartesian, direct-lattice and Miller-plane-normal conventions and strict selection validation.
- Mass collection now includes named-projection spectra and descriptors when stored definitions are available, and rejects inconsistent definitions.
- Added Matplotlib fatband/spectrum/mass plotting without requiring CATPlot. Ordinary projection columns can supply fatband size and color without families.
- Updated the standalone bundle from the same implementation, documented the unified workflow, and added an evidence-based release roadmap.
- Included the standalone scripts and synchronization utility in the source archive.

## 0.1.0

First release of `cat-tooling`: the preparation, calculation-management and scientific-analysis packages in this repository. CATPlot remains a separate distribution under `4_Plotting`.

### Phono3Py calculation management

- JDOS sampling runs one temperature per native invocation, avoiding the Phono3Py 2.7.0 early-return writer bug. Explicit lists and ranges are supported.
- JDOS auditing checks temperature coverage, finite three-column spectra, monotonic frequencies and requested row counts. Recursive scans detect duplicate filenames. Recovery generates only missing/bad temperatures.
- Grüneisen sampling produces one full-mesh job instead of competing GP jobs.
- Sampling and collection use an import-based runtime version probe, mirror native output to Slurm, append calculation logs, and propagate failures.
- Generated commands adapt force-constant reading flags to the Phono3Py 4 CLI. A version-specific workaround allows FC2-only JDOS on 4.4.0 despite its unrelated FC3 prerequisite check.
- Recovery preserves requested resources, modules and raw arguments.
- RTA/Wigner collection reuses intrinsic linewidths for different boundary lengths and writes into separate, fresh collection directories. It reads audited per-GP files rather than stale aggregate results.

### Installation and release checks

- Python 3.12 is selected by `.python-version`; the checked-in `uv.lock` is verified through a fresh installation with the documented extras.
- `cat-tools --version` and `cat-prepare --version` expose the package version.
- Release documentation covers installation, active-job upgrades, validation, and building/installing wheel and source distributions.

Upgrading CAT does not patch scripts already generated or queued. Generate new recovery jobs only after the corresponding existing jobs have finished.
