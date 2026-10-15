## CAT — Computational Analysis Toolbox

CAT brings together calculation preparation, scientific analysis, and visualization for materials research. This PhononTools repository contains the toolbox and its workflows.

The first release is **CAT 0.1.0**. See [CHANGELOG.md](CHANGELOG.md) for changes and [RELEASE.md](RELEASE.md) for the locked installation, verification and build procedure. The [release roadmap](ROADMAP.md) distinguishes pending validation from next-version work. For directional phonons, start with the unified [FreqBall workflow](2_Workflows/10_LatticeDynamics/01_Harmonic/02_DirectionalDOS/README.md): `cat-freqball --help` covers spectra, projections, fatbands, balls and mass scans.

The folders describe each program's purpose, independently of whether it runs locally or on an HPC system:

| Folder | Contents |
| --- | --- |
| [0_bashrc.d](0_bashrc.d/) | Interactive shell helpers |
| [1_Tools](1_Tools/README.md) | Input preparation, calculation management, structure conversion, and shared Python utilities |
| [2_Workflows](2_Workflows/README.md) | Procedures connecting preparation, calculation, and analysis |
| [3_Analysis](3_Analysis/README.md) | Scientific post-processing, fitting, visualization, examples, and plotting recipes |
| [4_Plotting](4_Plotting/README.md) | CATPlot: table-driven scientific figures and plotting recipes |

The Python code copied from `consol` is one installable distribution, `cat-tooling`. Its calculation utilities and preparation packages live under `1_Tools/CAT`; the analysis package lives under `3_Analysis/src`. Install from this repository's root:

```bash
uv sync --locked --extra structure --extra bonding --extra ordering --extra visualization
source .venv/bin/activate
cat-tools --version
cat-tools --help
cat-prepare --help
cat-phonopy --help
```

Install CATPlot in a separate tool environment while keeping CAT's locked dependencies:

```bash
uv tool install --python 3.12 --editable './4_Plotting[io]'
uv tool update-shell
```

Open a new shell if the tool directory was added to PATH, activate CAT's `.venv` again, and run `catplot --help`. CAT's lock selects pandas 3; CATPlot currently requires pandas 2. The separate tool environment lets both commands use their supported dependencies. CATPlot reads the exported CSV/Parquet tables. The [bonding README](3_Analysis/src/cat_tools/bonding/README.md) walks through LOBSTER extraction and plotting together.

For a single environment without the uv lock, install both distributions together with `python -m pip install -e '.[structure,bonding,ordering,visualization]' -e './4_Plotting[io]'`. The resolver selects versions satisfying both packages; this environment differs from `uv.lock`.

The example YAML files are templates: update their input paths and scientific settings for your calculation. Paths inside a YAML file are relative to that file. MPI is optional; install the `mpi` extra only in an environment with an MPI runtime.

On a cluster, run `uv python install 3.12` and then `uv sync --locked` from this repository. The checked-in `.python-version` selects Python 3.12; `uv.lock` fixes the dependency versions. `uv run --no-sync cat-tools ...` uses that installation without changing it. Add extras using the same `--locked` command. Do not sync an environment while running jobs are importing code from it; use a separate checkout/environment for an upgrade.

CAT's environment generates and audits jobs. The generated job loads its own Phono3Py module (or the runtime you select with `--no-module`). The locked analysis environment uses Phono3Py 4.4; generated jobs also support the legacy 2.7 command line. JDOS runs one temperature per invocation automatically, so the 2.7 multi-temperature writer bug does not lose higher-temperature output.

`cat-tools` groups VASP utilities, Phono3Py job preparation/auditing, and the preparation commands. `cat-prepare` remains a direct entry point for structure conversion, LOBSTER input generation, SoftBV/BVPA launching, and FreqBall mass-scan setup. Existing `cat-*` analysis commands remain available. The original `consol` directory is not needed at runtime. [SlurmScan](1_Tools/General/README.md) remains a standalone tool.

LOBSTER preparation uses an explicit YAML configuration for basis functions and interaction requests. `cat-prepare topology` derives periodic contacts from selected framework atoms; `cat-prepare lobsterin --config config.yaml -o lobsterin` renders the resulting input. See [the configuration guide](1_Tools/README.md#lobster-configuration-and-framework-detection).

The former standalone VASP, Phono3Py, and bonding implementations are retained in [1_Tools/legacy](1_Tools/legacy/README.md) for comparison during the migration to CAT. These archived copies belong in the repository for now, but are excluded from installation and receive no ongoing updates. Once CAT has been verified against the workflows they supported, the archive can be removed in a later commit; Git history will preserve it.

Run the tests from the repository root after installation:

```bash
python -m unittest discover -s tests -v
```

CATPlot recipe validation is skipped when the separate `catplot` package is absent.

Occupational ordering is available through `cat-ordering` or `cat-tools ordering`. Install the `ordering` extra for icet enumeration and SQS. The [ordering workflow](2_Workflows/00_Common/05_Ordering/README.md) covers independent substitutions, intact polyhedral motifs, supercell patterns, symmetry reduction, candidate selection, and mixing hulls built from calculated energies.

The [ordering command reference](1_Tools/CAT/ordering/README.md) documents model options and limitations; the workflow gives the sequential calculation procedure.
