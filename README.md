## CAT — Computational Analysis Toolbox

CAT brings together calculation preparation, scientific analysis, and visualization
for materials research. This PhononTools repository contains the toolbox and its workflows.

The folders describe each program's purpose, independently of whether it runs locally or on an HPC system:

| Folder | Contents |
| --- | --- |
| [0_bashrc.d](0_bashrc.d/) | Interactive shell helpers |
| [1_Tools](1_Tools/README.md) | Input preparation, calculation management, structure conversion, and shared Python utilities |
| [2_Workflows](2_Workflows/README.md) | Procedures connecting preparation, calculation, and analysis |
| [3_Analysis](3_Analysis/README.md) | Scientific post-processing, fitting, visualization, examples, and plotting recipes |
| [4_Plotting](4_Plotting/README.md) | CATPlot: table-driven scientific figures and plotting recipes |

The Python code copied from `consol` is one installable distribution, `cat-tooling`.
Its calculation utilities and preparation packages live under `1_Tools/CAT`;
the analysis package lives under `3_Analysis/src`. Install from this repository's root:

```bash
uv sync --extra structure --extra bonding --extra visualization
source .venv/bin/activate
cat-tools --help
cat-prepare --help
cat-phonopy 3_Analysis/examples/phonopy.yaml
```

Alternatively, install with `python -m pip install -e '.[structure,bonding,visualization]'`.
To include CATPlot in the same environment, install both local packages together:

```bash
python -m pip install -e '.[bonding]' -e ./4_Plotting
catplot --help
```

For a uv-managed environment, use `uv pip install -e ./4_Plotting` after syncing CAT.
CATPlot is a separate local distribution; `uv sync` may remove it unless reinstalled.
The [bonding README](3_Analysis/src/cat_tools/bonding/README.md) walks through
LOBSTER extraction and plotting together.

The example YAML files are templates: update their input paths and scientific settings
for your calculation. Paths inside a YAML file are relative to that file.
MPI is optional; install the `mpi` extra only in an environment with an MPI runtime.

`cat-tools` groups VASP utilities, Phono3Py job preparation/auditing, and the preparation commands. `cat-prepare` remains a direct entry point
for structure conversion, LOBSTER input generation, SoftBV/BVPA launching, and
FreqBall mass-scan setup. Existing `cat-*` analysis commands remain
available. [MIGRATION.md](MIGRATION.md) records the copied code's ownership and
compatibility entry points. The original `consol` directory is not needed at runtime.
[SlurmScan](1_Tools/General/README.md) remains a standalone tool.

LOBSTER preparation uses an explicit YAML configuration for basis functions and
interaction requests. `cat-prepare topology` derives periodic contacts from selected
framework atoms; `cat-prepare lobsterin --config config.yaml -o lobsterin` renders
the resulting input. See [the configuration guide](1_Tools/README.md#lobster-configuration-and-framework-detection).

The former standalone VASP, Phono3Py, and bonding implementations are retained in
[1_Tools/legacy](1_Tools/legacy/README.md) for comparison during the migration to CAT.
These archived copies belong in the repository for now, but are excluded from
installation and receive no ongoing updates. Once CAT has been verified against
the workflows they supported, the archive can be removed in a later commit;
Git history will preserve it.

Run the tests from the repository root after installation:

```bash
python -m unittest discover -s tests -v
```

CATPlot recipe validation is skipped when the separate `catplot` package is absent.
