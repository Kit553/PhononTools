# LOBSTER preparation and analysis

Prepare `lobsterin` from the POSCAR used by the LOBSTER calculation:

```bash
cat-prepare lobsterin POSCAR -o lobsterin
```

Set basis functions, species, and interaction ranges for your chemistry using the
generator's options or a YAML/JSON configuration. Keep the generated
`lobsterin.bonding.json` beside the input: the analysis uses its topology and provenance.

Run your VASP/LOBSTER calculation with the appropriate external programs and job
settings. Then adapt [bonding.yaml](../../3_Analysis/examples/bonding.yaml) to the
result directory and run:

```bash
cat-bonding bonding.yaml
```

The analysis exports interaction tables and validates the requested interactions
against the available outputs. See [BondingTools](../../3_Analysis/docs/bonding.md).
