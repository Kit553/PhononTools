## System Convergence

**Versioning:** This workflow is built around the Cluster installations of `VASP`(v6.2.0) for electronic and structural convergence, with `PhonoPy` (v2.20.0) used for phonon calculations.

**Use:** Establish numerically stable calculation parameters while retaining a reasonable computational cost.

Convergence should be established as the step 0 of any computational project. This is especially true for phonon calculations, which depend on the derivatives of the PES and are therefore particularly sensitive to insufficiently converged forces.

The main convergence axes considered here are:
- Plane-wave basis completeness via `ENCUT`
- Brillouin-zone sampling via electronic *k*-point density
- Electronic Convergence via `EDIFF`
- Ionic convergemnce via `EDIFFG`
- Electronic smearing via `SIGMA`; this is most relevant for metallic systems
- Number of calculated bands via `NBANDS`; this is relevant for projection quality and bonding analysis
- Finite-size effects and supercell convergence for property-specific calculations

---

**Folder Setup:**
The working directory should contain:
- `POSCAR`   &rarr; your initial strucutre from crystallographic data etc.
- `POTCAR`   &rarr; the corresponding potentials
- `INCAR.d`  &rarr; generic base settings file
- `KPOINTS`  &rarr; generic *k*-mesh settings file

**Helper Scripts:**
These should be callable for your workflow e.g. in your `/bin/` folder:

`KPOINTSGen.sh`

Generates a reasonable first-pass $Gamma$ centered *k*-mesh from an arbitrary `POSCAR` for the initial convergence of `ENCUT`

The default reciprocal-density criterion used is approximaetly equivalent to

[
N_iL_i \approx 50~\mathrm{\AA}
]

for orthogonal cells. For non-orthogonal cell the reciprocal lattice vectors are evaluated explicitly. 
Using this script, different density cirterions and `POSCAR`/`CONTCAR` files can be specified via

```bash
./KPOINTSGen.sh INPUTFILE DENSITY
```

be mindfull of heavily skewed triclinic systems as these can be problematic with regard to a defined density.
This is only meant as a first reasonable *k*-mesh to enable convergence of earlier parameters, the mesh itself needs to be converged at a later point.

---

`ConvGen.sh`

Generates and submits the convergence calculations. 
The general syntax for range convergence is:

`ConvGen.sh --PARAMETER --range LOW HIGH STEP`

or for manual value convergence:

`ConvGen.sh --PARAMETER --man VALUES`

For example the convergence of `ENCUT` can be done via:

```bash
ConvGen.sh --encut --range 1.0 2.0 0.25
```

For `ENCUT`, range values are multiples of the largest `ENMAX` detected automatically in the `POTCAR`.

Alternatively manual values can be supplied as:

```bash
ConvGen.sh --encut --man 300 350 400 450 500 550
```

This script creates one directory per calculation, copies the necessary `VASP` inputs, modifies the scanned parameter, adjusts the SLURM job name and submits the calculations.

**Additional options:**
`--calc static`    &rarr; static calculation; default
`--calc ion`       &rarr; relax internal coordinates
`--calc full`      &rarr; relax internal coordinates and cell parameters
`--calc keep`      &rarr; explicitly preserve settings from the `INCAR`

`--set TAG=VALUE`  &rarr; override an arbitrary `INCAR` tag; may be supplied repeatedly
`--no-submit`      &rarr; generate calculation directories without submitting them
`-v, --verbose`    &rarr; print detailed setup/submission information

**Supported Convergence Axes:**

| Flag        | Quantity                    | `--range` / `--man` values       |
|:------------|:----------------------------|:---------------------------------|
| `--encut`   | Plane-wave cutoff           | ENMAX multipliers / absolute eV  |
| `--kpoints` | Reciprocal *k*-point density| Density criterion                |
| `--ediff`   | Electronic SCF convergence  | `EDIFF` values                   |
| `--ediffg`  | Ionic-force convergence     | `EDIFFG` values                  |
| `--sigma`   | Smearing width              | eV                               |
| `--nbands`  | Number of bands             | Integer band counts              |

`EDIFFG` should normally be scanned with `--calc ions` or `--calc full`; it has no effect on static calculations and will be ignored by `VASP`.

---

`ConvCheck.sh`

Reads completed convergence directories and prints a compact comparison table.

For example:

```bash
./ConvCheck.sh --prefix ENCUT
```

The static table contains the quantities:

`Directory | ENCUT | ETOT | ETOT/Atom | MaxForceAtom | MaxCellStress | WallTime`

The directory prefix makes the tool reusable for all convergence checks:

```bash
./ConvCheck.sh --prefix ENCUT
```
```bash
./ConvCheck.sh --prefix KPOINTS
```
```bash
./ConvCheck.sh --prefix EDIFF
```
```bash
./ConvCheck.sh --prefix SIGMA
```

For relaxtion scans:

```bash
./ConvCheck.sh --prefix EDIFFG --relax
```

The script instead compares the starting `POCAR` with the final `CONTCAR` and reports the structural shift.

---

## Workflow

### 1. Plane-Wave Basis Convergence

Generate a sensible starting grid:

```bash
./KPOINTS_Generate.sh POSCAR 50
```

This grid will only be used as a provisional mesh while converging the plane-wave basis and will be discarded afterwards.

Next run a **static** `ENCUT` scan: 

```bash
./ConvGen.sh --encut --range 1.0 2.5 0.25
```

and inspect the output via:

```bash
./ConvCheck.sh --prefix ENCUT
```

Monitor:
- total energy per atom
- atomic forces
- cell stress
- Pullay stress
- wall time

and pick the lowest `ENCUT` value for which these quantities plateau. You'll notice dramatic differences between the wall time and memory usage of higher `ENCUT` jobs.
`VASP` officially recommends a value around $1.3\times\max(\mathrm{ENMAX})$ while group standard is a more conservative $2\times\max(\mathrm{ENMAX})$.

If your goal are phonon calculations the **forces** and **stress** are more important than total energy convergence alone because an incomplete basis may introduce Pullay stress on the cell.

---

### 2. *k*-Point Density:

Fix the converged `ENCUT` value and vary only the reciprocal-space density:

```bash
./ConvGen.sh --kpoints --range 30 70 10
```

Inspect with

```bash
./ConvCheck.sh --prefix KPOINTS
```

Choose the lowest density grid for which energy, forces and stress are essentially stable.

---

### 3. Remaining Electronic Parameters:

With `ENCUT` and the *k*-mesh fixed, converge parameters as required by the job you have in mind.

Typcial examples with ranges are:

**EDIFF:**
This is particularly usefull if you want to min-max the timing of your calculations

```bash
./ConvGen.sh --ediff --man 1e-4 1e-5 1e-6 1e-7 1e-8
```

**Sigma:**
This is essential if you are working on a metallic/semiconducting system, non-converged `SIGMA` can make electronic structure calculations inaccurate.

```bash
./ConvGen.sh --sigma --man 0.20 0.10 0.05 0.02
```

**NBANDS**
This usually requires little attention in ordinary ground-state DFT but should be explicitly checked for workflows requiring a substantial unoccupied-state manifold, such as LOBSTER or spectroscopy calculations

```bash
./ConvGen.sh --nbands --range 100 400 50
```

---

### 4. Preliminary Relxation

After the static problem is sufficiently optimized, relax the system with those settings.

```bash
./ConvGen.sh --ediffg --man -0.02 -0.01 -0.005 -0.002 \
             --calc full
```

and check after completion

```bash
./ConvCheck.sh --prefix EDIFFG --relax
```

Compare how strongly the resulting structure differs as the force criterion is thightened.

---

### 5. Parameter Validation

After the first-pass relaxation, verify that the chosen parameters remain well converged. Usually it is sufficient to compare the selected production settings against one stricter calculation rather than repeating the complete intial scans. 

The final structure can then be relaxed with the selected production settings. Most of the time the settings are transferable to chemically similar systems and there is no need to re-converge absolutely everything for each system.

---

### 6. Property Specific Convergence Checks

After this additional property-specific convergence tests are required based on your workflow. Common examples include:#
- phonon supercell size
- phonon displacement amplitude
- *q*-mesh density
- third-order interaction cutoff
- NEB image count
- MD cell size; trajectory lenght or step size
- spectroscopy-specific empty-states requirements

These should only be considered once the underlying electronic problem is stable. All of these are highly specific to the tools you use and the quantity you have in mind. As such no general criteria or workflows are provided here. 





