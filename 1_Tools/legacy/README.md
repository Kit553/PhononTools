## Archived Tools

These full standalone implementations preserve the former script layout for reference. They are snapshots at retirement, including the CLI fixes made during CAT integration; they are not compatibility launchers or byte-for-byte pre-migration originals. They are excluded from CAT installation and receive no ongoing updates. Use [CAT](../README.md) for maintained commands.

Keep this directory in the repository during the migration to support comparison with CAT. Make ongoing fixes in the maintained implementations. Once CAT has been verified against the workflows these scripts supported, this directory will be removed.

| Folder | Contents | Maintained replacement |
| --- | --- | --- |
| [VASP_Workflow](VASP_Workflow/) | `POTCAR`, `KPOINTS`, convergence, restarts, comparison, and QHA rescaling | `cat-tools vasp` |
| [Phono3Py](Phono3Py/) | Grid-point sampler and evaluator | `cat-tools phono3py` |
| [consol_bonding](consol_bonding/README.md) | Historical LOBSTER generators | `cat-prepare lobsterin` |

The archived Phono3Py evaluator uses its adjacent sampler. These scripts retain original calculation and submission behavior: convergence generation submits by default. 

SlurmScan remains maintained independently in [General](../General/README.md).
