# CAT release roadmap

This is a code-and-test audit of the current working tree, not a promise of dates. The first release remains 0.1.0; the next feature milestone is provisionally 0.2. Packaging and installation instructions live in [RELEASE.md](RELEASE.md).

## Current priority — finish the existing workflows

The next-release proposal below is deferred until the maintained workflows cover at least the scientific tasks in `OLD_Workflows`. The active completion record is [2_Workflows/COVERAGE.md](2_Workflows/COVERAGE.md), including QHA, anharmonic, electronic-structure and mode-analysis procedures. Document actual CAT/upstream capabilities and their limits before adding new release scope.

## First-release scope

| Area | Available now | Evidence / procedure |
| --- | --- | --- |
| FreqBall | One command for spectra, ordinary projections, family fatbands, descriptors, balls and mass scans; compatible standalone bundle | [Workflow](2_Workflows/10_LatticeDynamics/01_Harmonic/02_DirectionalDOS/README.md), [reference](3_Analysis/docs/freqball.md), `test_freqball*.py` |
| Direction conventions | Cartesian vectors, direct-lattice axes and Miller plane normals share spectra/mode definitions; mass collection verifies those definitions | `test_freqball_unified.py`: oblique-cell geometry, native Phonopy PDOS, x+y+z sum rules, portable reanalysis, serial/parallel collection |
| Phono3Py jobs | Single-temperature native JDOS calls, strict audit, missing-temperature recovery and isolated boundary-MFP collection | [CAT utilities](1_Tools/CAT/README.md), `test_phono3py_jobs.py` |
| Installation | Checked-in uv lock, optional extras, wheel/source packaging and generated-source synchronization | [Release procedure](RELEASE.md), `.github/workflows/test.yml`, `test_freqball_standalone.py` |
| Broader toolbox | Harmonic/anharmonic analysis, structure, ordering, bonding, thermal and impedance tools with examples and tests | [Analysis index](3_Analysis/README.md), [preparation index](1_Tools/CAT/README.md) |

The tensor is the canonical raw directional-spectrum result. Mode tables still need force-data eigenvectors. Conduction fatbands describe harmonic mode character along chosen structural axes; they do not predict ionic conductivity or hopping rates.

## Validation before announcing a production release

These checks need deployment evidence beyond local synthetic fixtures. Record the commands, environment versions and result locations with the release record.

| Priority | Gap | Completion criterion |
| --- | --- | --- |
| P0 | Real MPI/Slurm execution has not been established by the simulated collective tests | Run one small serial and 2/4-rank calculation on the target cluster; compare tensors/tables within declared tolerances, test a rank-local input failure, and record memory/walltime. |
| P0 | Legacy 2.7 JDOS routing is tested with a reproducer; native numerical integration uses the locked 4.4 stack | Run a small two-temperature job in the actual Palma 2.7 environment; audit both spectra and logs, then exercise missing-temperature recovery. |
| P0 | Na3PSe4 production recovery is separate from toolbox unit tests | Once active jobs finish, audit corrected-`BORN` outputs and validate JDOS temperature coverage, Grüneisen output and boundary-MFP collections. Do not regenerate into directories used by running jobs. |
| P0 | CI configuration and local checks do not establish a successful remote release build | Obtain green Python 3.12/3.13 CI, build from the reviewed commit, inspect/install its wheel and source archive, then record checksums. |
| P1 | Headless 3D rendering depends on the deployment's VTK backend | Export a ball screenshot/GIF on the target workstation/cluster and document the working backend; ordinary table/2D workflows remain usable without it. |

## Proposed 0.2 release plan

**Theme: reproducible, restartable workflows from cluster jobs to comparable scientific results.** This is a mock scope proposal, not an implementation claim or a delivery-date commitment. Complete milestones in dependency order; use the first real calculation campaigns to refine the scope.

### Milestone 1 — Establish the release baseline

- Finish the production-validation checks above and record supported combinations of CAT, Python, native Phonopy/Phono3Py, MPI and the scheduler environment.
- Add a small redistributable polar-material fixture alongside synthetic tests. Include input hashes, expected numerical tolerances and a short reference run.
- Preserve the 0.1 CLI/configuration examples as compatibility tests. Define how configuration and artifact schema versions will evolve independently of package versions.

**Exit:** a clean checkout installs from the lock and reproduces the reference workflow; serial and real-MPI results agree within documented tolerances.

### Milestone 2 — Make jobs traceable and restartable

- Introduce a shared run record with workflow adapters for FreqBall mass scans and Phono3Py sampling/recovery. Record input/config hashes, resolved commands, working directories, executable versions, output expectations and log paths.
- Record the submission ID and immutable array-index-to-task mapping when CAT submits a job; support importing identifiers for externally submitted jobs. Resolve walltime, memory and exit status through scheduler accounting when available, with an explicit unknown state when records are unavailable.
- Add status and resume planning that distinguishes absent, incomplete, failed, stale and validated results. Show the proposed recovery set before an explicit run or submission action. Reuse the existing strict JDOS/GP audits.
- Define invalidation per scientific dependency: changed `BORN`/force data must invalidate affected numerical results; a new projection can reuse a compatible spectral tensor; a new boundary MFP can reuse compatible intrinsic linewidths.

**Exit:** an interrupted two-temperature/multi-mass fixture resumes only the work that needs it. Changed physical inputs cannot silently authorize old outputs. Every array task can be traced to its command, logs and accounting record.

### Milestone 3 — Give workflows a consistent front end

- Add early validation of versioned configurations: unknown keys, unavailable sites, contradictory options, paths and required dependencies fail before an allocation is used. Preserve legacy normalization through explicit migration.
- Add an environment report distinguishing the CAT analysis interpreter from the native executable/module stack used by generated jobs.
- Add FreqBall mode-only calculation using the existing projection model, so a fatband calculation does not require a frequency tensor first.
- Align error messages and common output/status conventions across the adapters; retain existing commands rather than requiring a package-wide CLI migration.

**Exit:** users can validate a workflow without running it, diagnose the selected environment, and produce identical mode tables through the full or mode-only path.

### Milestone 4 — Make comparisons reproducible

- Provide convergence summaries for mesh, Gaussian width, frequency window and directional sampling. Keep sampled-direction extrema distinct from exact named projections, which do not depend on sphere resolution.
- Report captured spectral weight and imaginary-frequency contributions with explicit definitions. Gaussian weight below zero alone must not be labeled an unstable mode; instability checks need the underlying mode frequencies.
- Add projection/mass filters and common plot scales for comparison figures. Record source-table hashes and plotting settings in a figure sidecar.
- Supply one complete walkthrough from harmonic inputs through projections, mass comparison and figures, plus one Phono3Py recovery example.

**Exit:** a colleague can reproduce a comparison from its inputs and manifests, identify the selected numerical settings and regenerate the figures.

### Milestone 5 — Qualify and package 0.2

- Run the compatibility suite and cluster smoke tests; exercise interruption, changed-input and partial-result cases explicitly.
- Build/install the wheel and source distribution outside the checkout. Check the standalone bundle against shared sources and run its copied-bundle test.
- Validate documentation links/examples and ship a migration guide, tested-version matrix, release notes, checksums and benchmark record.

**Exit:** a release candidate passes the documented checks and a real campaign uses status/recovery successfully. Tagging and publishing are separate release actions; this roadmap does not perform them.

### Scope boundary

Milestones 1–3 and release qualification are the core 0.2 scope. For milestone 4, the minimum is reproducible comparison plots and the worked examples; a general convergence-sweep orchestrator can move to 0.3 if it threatens the release.

Keep these larger scientific additions as **0.3 candidates**:

- Mode tracking across mass changes, including crossings and degenerate subspaces.
- Directional group-velocity/conductivity views with explicit vector/tensor transformations and validated identity matching between analysis modules.
- Ordering extensions such as vacancies, shared motifs or external enumeration backends; select one independently testable extension rather than bundling all.

Avoid making a graphical application, a general scheduler framework or additional solver backends prerequisites for 0.2. Revisit them after the common workflow and artifact contracts have been exercised in practice.

## Detailed backlog and acceptance tests

| Priority | Work | Current limitation and acceptance test |
| --- | --- | --- |
| P1 | Mode-only FreqBall action | `run` computes a tensor first; `analyze` expects one. Add a mode-only action using the existing projection definitions, with equality to current band/mesh tables and no HDF5 prerequisite. |
| P1 | Validated mass-scan resume and status | `run` executes every job and collection rejects recorded failures. Resume should select only missing/failed points after checking input/config hashes; interruption must never mark stale outputs complete. Test interrupted and changed-input runs. |
| P1 | Early configuration validation | FreqBall checks projection selections and families before computation, but does not have a complete versioned YAML schema. Reject unknown settings consistently and provide a cheap validation command before allocations. Keep v6 normalization regression coverage. |
| P1 | Explicit convergence reports | Users currently inspect mesh, broadening, frequency-window and direction-count effects manually. Report clipped spectral weight, instability weight and descriptor convergence on a reproducible reference series. Keep physical tolerance choices explicit. |
| P1 | Small shareable scientific reference dataset | Current tests mainly use synthetic fixtures. Add a redistributable polar material fixture with known input hashes, tensor/fatband references, NAC limits and a mass scan; record expected tolerances and runtime. |
| P2 | Mass-dependent mode tracking | Each point has its own fatband table. A common band index is not a tracked physical mode across crossings or degeneracies. Add overlap/subspace tracking and explicit ambiguous assignments before exporting cross-mass mode trends. |
| P2 | Velocity and transport direction views | FreqBall projects displacement eigenvectors; harmonic group-velocity and anharmonic conductivity analysis remain in their respective modules. Reuse their table keys/units and define vector/tensor transformations before offering a unified directional-transport view. |
| P2 | Plot provenance and comparison controls | Tables have manifests; quick Matplotlib figures currently do not have source/settings sidecars. Record source hashes, plot options and common scales, and add selection controls for dense named-projection mass plots. |
| P2 | Automated documentation checks | Add local-link and CLI-example checks, a docs landing page and an explicit tested-version matrix. The standalone source-sync check already prevents code drift. |

## Broader toolbox follow-up

- **Ordering:** vacancies, positional/orientational disorder, shared-polyhedron constraints, Supercell import and a DFT+U compatibility policy remain explicit limits in the [ordering reference](1_Tools/CAT/ordering/README.md). Each needs its own physical model and fixtures; do not broaden accepted inputs without them.
- **Compatibility:** distinguish the locked CAT analysis environment from native executables selected by generated cluster jobs. Expand a tested compatibility matrix before widening supported version claims.
- **Distribution:** keep the standalone FreqBall exporter, compatibility launchers and packaged CLI synchronized. Publish generated bundles from a reviewed commit so users can identify the exact source used for a calculation.

Real MPI, legacy native runs and production Na3PSe4 validation are outstanding deployment checks. They are not replaced by passing unit tests or a successful build.
