## CAT Implementation and Python Interfaces

**Use:** Maintain and reuse the code behind `cat-tools`, `cat-prepare`, and the analysis commands. For normal command usage, start with the [1_Tools command reference](../README.md).

**Folder Setup:**

| Package | Purpose | Entry point |
|:--------|:--------|:------------|
| [cat_utilities](cat_utilities/) | VASP preparation/restarts, Phono3Py job generation/auditing. | `cat-tools` |
| [cat_prepare](cat_prepare/) | Structure conversion, config-driven LOBSTER inputs and periodic framework topology, BVPA, and FreqBall mass-scan setup. | `cat-prepare` or `cat-tools prepare` |
| [cat_shared](cat_shared/) | Structure/frame I/O, table writing, mode joins, topology, and shared parsers. | Python API |
| [cat_tools](../../3_Analysis/src/cat_tools/) | Scientific analysis and visualization. | Domain commands such as `cat-phonopy` |

---

## Installation and Ownership

`CAT` replaces the previous `1_Tools/src` source directory. Install from the repository root using its `pyproject.toml` and `uv.lock`; there is no separate installation step for this folder. Package discovery includes `1_Tools/CAT` and `3_Analysis/src` in the same `cat-tooling` distribution.

VASP and Phono3Py shell implementations live beside their related Python modules under `cat_utilities`. They are included in both wheels and source distributions and are invoked with Bash. Former standalone implementations are preserved in [legacy](../legacy/README.md), outside package discovery. The q-mesh evaluator uses its bundled sampler by default, so no separate `qmesh_sampling` alias is needed. [SlurmScan](../General/README.md) remains independently usable of CAT.

The `cat-tools` dispatcher loads a Python command only when selected and passes arguments to its `main(argv)` function. It forwards shell arguments as an argument list and preserves the selected tool's exit status and working directory.

---

## Reusing Functions

The original numerical functions remain importable from their owning modules:

```python
from pathlib import Path
from cat_utilities.vasp.compare import read_poscar, match_by_species
from cat_utilities.phono3py.evaluate import read_expected_gps
from cat_prepare.lobsterin import generate_lobsterin, LobsterinSettings
from cat_prepare.mass_scan import setup_mass_scan
from cat_shared.structure.io import read_structure, write_structure

reference = read_poscar(Path("POSCAR"))
relaxed = read_poscar(Path("CONTCAR"))
comparison = match_by_species(reference, relaxed, reference.lattice)
```

To invoke the same parser used by the CLI without changing `sys.argv`:

```python
from cat_utilities.vasp.rescale import main

status = main(["POSCAR", "--iso", "5", "--dry-run"])
```

Argument errors and `--help` follow normal argparse behavior and may raise `SystemExit`. For numerical work, prefer the module's individual functions.

The `POSCAR` comparison module uses bottleneck matching to pair atomic sites. `cat_tools.structure.compare` uses Pymatgen's structure matcher. Both use the structure representations in `cat_shared`.

---

## Development Checks

After installing the required extras, run from the repository root:

```bash
python -m unittest discover -s tests -v
uv lock --check
cat-tools --help
cat-tools vasp --help
cat-tools phono3py --help
```

Tests use temporary calculation fixtures and do not submit jobs. Shell resources must be retained in `pyproject.toml` package data and `MANIFEST.in` when adding a new shell command. Add dispatch entries in `cat_utilities/cli.py` and document the command's inputs, writes, overwrite behavior, and submission behavior in the parent README.

`tests/test_phono3py_jobs.py` executes generated scripts against a fake native program that reproduces the 2.7 JDOS writer bug. It also checks error propagation, temperature recovery, runtime version handling, and isolation of conductivity collections. [The release procedure](../../RELEASE.md) covers a fresh locked environment and installation of the built wheel.

This documentation was automatically generated based on the code using gen AI (GPT-6 Astra).

LOBSTER configuration and the topology-to-YAML workflow are documented in the [Tools reference](../README.md#lobster-configuration-and-framework-detection). Framework roles and basis functions are supplied explicitly for each system.

The topology helper also supports nearest or all ligand contacts between connected components, distinguishing bonded polymer repeats from neighboring chains. See the [inter-component guide](../README.md#inter-polyanion-and-inter-chain-contacts) and [polymer seed](examples/lobsterin.polymer.seed.yaml).

Ordering generation, symmetry reduction, SQS, energy collection and hull schemas are documented in the [ordering reference](ordering/README.md). Follow the [numbered ordering workflow](../../2_Workflows/00_Common/05_Ordering/README.md) for a procedure using your own structures, directories and execution environment. The [JCH material inventory](ordering/JCH_INVENTORY.md) is kept beside that reference.

Optional [VASP site notes](vasp-site-notes.md) and [benchmark case studies](examples/vasp-benchmarks/README.md) preserve concrete installation and material examples separately from the general workflows.
