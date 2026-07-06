#### 1.1 Building your `INCAR` file:
The INCAR gives you full control over what VASP does with your supplied structure. This section presents a small overview over possible tags, a full list can be found at the [VASP Wiki](https://vasp.at/wiki/Category:INCAR_tag).

**System and Output parsing:**
- `NWRITE` &rarr; how much information is written to the `OUTCAR` file, leave at `NWRITE=2` for most purposes [Wiki](https://vasp.at/wiki/NWRITE)
- `SYSTEM` &rarr; sets the title string of the run; use your system identifier [Wiki](https://vasp.at/wiki/SYSTEM)
- `LCHARG` &rarr; <Bool>; determines wheter the `CHGCAR` is written. Use when restarting similar systems from a reused `CHGCAR`
- `LWAVE`  &rarr;
- `GGA`    &rarr;

**Parallelization:**
- `NCORE`
- `KPAR`

**Electronic Structure:**
- `Prec`
- `EDIFF`
- `ENCUT`
- `ISMEAR`
- `ALOG`
- `NELM`
- `NELMIN`
- `NELMDL`
- `LASHP`
- `ADDGRID`

**Ionic Relaxation:**
- `EDIFFG`
- `IBRION`
- `NSW`
- `ISIF`
- `IWAVPR`
