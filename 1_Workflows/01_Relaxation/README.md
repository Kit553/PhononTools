## Relaxation:
**Use**: Obtain the local athermal minimum of the structure as a requirement for most other workflows.

## Workflow
### 0. Obtain the initial structure guess:
To start the relaxation a structure must be supplied as the starting point, preferrably use already computed structures known to work or low-temperature structural data. The less thermal effects need to be reversed the better.

##### 0.1 From Experimental data:
 - Get a cif from Pearson Crystal Database (PCD), Materials Project, Inorganic Crystal Structure Database or experimental data
 - Open the cif in VESTA, then under File > Export > VASP structure file

#### 0.2 From previous calculations:
 - Quantum Materials Open Database directly provides the POSCAR necessary from previous users
 - Using the /Structure/Conversion/POSCAR-ToCif.py tool can convert it back to cif format to inspect it or adjust parameters

**Folder Setup**
- `POSCAR`  &rarr; POSition file containing all atom species and their fractional coordinates
- `POTCAR`  &rarr; POTentials used in the simulation, can be obtained using the `apgu` script in the [tools]([https://github.com/user/repo/blob/branch/other_file.md) section
- `INCAR`   &rarr; INput control file; see below for a list of the most important tags
- `KPOINTS` &rarr; Brillouin zone integration mesh parameters; recommended to use roughly 50 points per $A^{-1}$
- `job.sh`  &rarr; job submission script for PALMA-II

### 1. INCAR settings:
Below a list of the most important settings for ionic relaxation in the `INCAR` is given:
- `EDIFFG` &rarr; force convergence criterion, usually set to `EDIFFG = -1E-03` so all forces are below $-1E-03$ eV per Angström
- `NSW`    &rarr; maximum number of ionic steps for one relaxation run
- `IBRION` &rarr; determines structure change during the relaxation, standard is conjugant gradient with `IBRION = 2`
- `ISIF`   &rarr; determines how the stress tensor is calculated and which degrees of freedom are considered, standard is `ISIF = 3` for all degrees of freedom
- `IWAVPR` &rarr; determines how the wavefunctions are extrapolated from one step to the next, standard is `IWAVPR = 1` for simple extrapolation based on charge densities

As the relaxed structure is the foundation of all further characterization it is generally recommended to go for a high degree of accuracy of this calculations, set sufficiently tight `EDIFF` and `EDIFFG` criteria.
Additionally, it can be helpful to consider aspherical contributions (`LASPH = .TRUE.`), use the support grid for augementation charges (`ADDGRID = .TRUE.`) and playing around with the smearing method (`ISMEAR`) to get the optimal starting point.

### 2. Determine the *k*-mesh: 
A good relaxation is strongly dependent on thorough sampling of the Brillouin zone, controled via the `KPOINTS` file. For insulators and semiconductors a good guideline is to construct the *k*-mesh based on reciprocal space.\

$$k-point lenght / reciprocal density \appox 40-50$$

While this *k*-spacing is sufficient in most cases, it is still necessary to check for convergence against denser grids to make absolutely sure that the calculation run is valid. See Workflow 00_Convergence for more details.

Make sure that any anisotropy of you unitcell is also reflected in your sampling grid.

#### 3. Running and Checking the calculation(s)
Once the `INCAR` and `KPOINTS` file have been adapted to suit your system, the calculations can be run via the job script `job.sh`. After the run is concluded the new structure will be written to the `CONTCAR` file while any other information on stresses etc. will be written to the `OUTCAR`.

To ensure the structure reaches an appropriate minimum, check said `OUTCAR` for the stress blocks shown below. Ensure Pullay stress is negliable as this will ruin any attempt at accurate phonon calculations.

Stress versus deformation of the cell:
```bash
  FORCE on cell =-STRESS in cart. coord.  units (eV):
  Direction    XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------------
  Alpha Z   150.37359   150.37359   150.37359
  Ewald   -1512.88780 -1512.88780 -1537.42255    -0.00000    -0.00000     0.00000
  Hartree   558.85252   558.85252   548.65389    -0.00000    -0.00000    -0.00000
  E(xc)    -424.72695  -424.72695  -424.72834     0.00000    -0.00000     0.00000
  Local    -459.99241  -459.99241  -426.45844     0.00000    -0.00000     0.00000
  n-local   551.88505   552.19628   551.56581    -0.22710    -0.18552    -0.33117
  augment    87.94374    87.94374    87.89045     0.00000    -0.00000     0.00000
  Kinetic  1049.61328  1047.18459  1050.12782    -0.44034    -0.51076    -0.44682
  Fock        0.00000     0.00000     0.00000     0.00000     0.00000     0.00000
  -------------------------------------------------------------------------------------
  Total       0.00229     0.00229     0.00224     0.00000     0.00000     0.00000
  in kB       0.01097     0.01097     0.01075     0.00000     0.00000     0.00000
  external pressure =        0.01 kB  Pullay stress =        0.00 kB
```

Stress versus displacement of the atoms:
```bash
FORCES acting on ions
    electron-ion (+dipol)            ewald-force                    non-local-force                 convergence-correction
 -----------------------------------------------------------------------------------------------
   -.442E+02 -.353E+02 -.376E+02   0.452E+02 0.368E+02 0.389E+02   -.103E+01 -.157E+01 -.128E+01   -.168E-03 0.137E-03 0.179E-03
   0.442E+02 0.353E+02 -.376E+02   -.452E+02 -.368E+02 0.389E+02   0.103E+01 0.157E+01 -.128E+01   0.168E-03 -.137E-03 0.179E-03
   -.353E+02 0.442E+02 0.376E+02   0.368E+02 -.452E+02 -.389E+02   -.157E+01 0.103E+01 0.128E+01   0.137E-03 0.168E-03 -.179E-03
   0.353E+02 -.442E+02 0.376E+02   -.368E+02 0.452E+02 -.389E+02   0.157E+01 -.103E+01 0.128E+01   -.137E-03 -.168E-03 -.179E-03
   0.442E+02 -.353E+02 0.376E+02   -.452E+02 0.368E+02 -.389E+02   0.103E+01 -.157E+01 0.128E+01   0.168E-03 0.137E-03 -.179E-03
   -.442E+02 0.353E+02 0.376E+02   0.452E+02 -.368E+02 -.389E+02   -.103E+01 0.157E+01 0.128E+01   -.168E-03 -.137E-03 -.179E-03
   0.353E+02 0.442E+02 -.376E+02   -.368E+02 -.452E+02 0.389E+02   0.157E+01 0.103E+01 -.128E+01   -.137E-03 0.168E-03 0.179E-03
   -.353E+02 -.442E+02 -.376E+02   0.368E+02 0.452E+02 0.389E+02   -.157E+01 -.103E+01 -.128E+01   0.137E-03 -.168E-03 0.179E-03
   -.116E-12 0.557E-12 0.427E+01   -.334E-13 0.108E-13 -.348E+01   -.154E-19 0.182E-18 -.788E+00   0.183E-12 0.355E-12 0.354E-04
   0.172E-12 -.418E-12 -.427E+01   0.439E-13 -.450E-13 0.348E+01   0.305E-19 -.875E-19 0.788E+00   0.275E-12 -.418E-12 -.354E-04
   0.718E-12 0.129E-13 -.427E+01   0.255E-13 -.928E-14 0.348E+01   0.332E-19 -.153E-18 0.788E+00   0.712E-13 0.531E-12 -.354E-04
   -.802E-12 0.845E-12 0.427E+01   -.581E-13 0.380E-13 -.348E+01   -.463E-19 0.618E-19 -.788E+00   0.993E-14 -.370E-12 0.354E-04
   0.191E-13 0.523E-12 -.336E-12   -.258E-13 -.291E-13 0.711E-11   -.593E-20 0.474E-19 -.489E-19   -.222E-12 0.634E-12 0.101E-11
   0.622E-12 -.241E-11 0.296E-11   0.395E-13 0.422E-13 -.107E-11   0.251E-19 0.457E-19 -.111E-19   0.585E-13 0.560E-12 -.915E-12
   0.159E-12 0.354E-12 -.127E-12   -.650E-14 0.876E-15 0.989E-12   0.219E-17 -.353E-17 0.462E-18   -.605E-13 -.256E-12 -.584E-12
   -.317E-13 -.250E-11 -.568E-12   0.223E-13 0.451E-13 0.103E-11   -.171E-17 -.772E-18 -.316E-18   -.625E-13 -.315E-12 0.548E-12
 -----------------------------------------------------------------------------------------------
   0.151E-06 0.151E-06 0.296E-07   0.217E-15 0.182E-13 0.113E-13   -.444E-15 0.440E-15 0.222E-15   0.156E-12 -.181E-13 0.191E-12
```

For thorough relaxation it is also good to backup and resubmit the first relaxed structure and check if further displacement takes place. To backup use the `varc` script [tools]([https://github.com/user/repo/blob/branch/other_file.md). This moves the structure of the prior relaxation run `CONTCAR` to `POSCAR` and saves a copy of the `OUTCAR` and `vasprun.xml`.
To compare the two structures use `vimdiff` or a similar feature to see if more movement took place, if so, resubmit the job and repeat.

### 4. Get the static energy:
After the structure is completely relaxed, change the `INCAR` settings to perform one static run to get the final energy of the system. In most cases it is sufficient to set `NS" = 0`, but to be fully sure, also change the settings of `IBRION`, and remove `ISIF`, `EDIFFG`, and `IWAVPR`.

Before running the job, make sure you have saved the last relaxation run properly by backing up with `varc`. Then resubmit the job and obtain the total energy at the bottom of the `OUTCAR` under: 

```bash
  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
 ---------------------------------------------------
 free  energy   TOTEN  =       -69.37340465 eV

 energy  without entropy=      -69.37340465  energy(sigma->0) =      -69.37340465
```

This same value can also be used for convex hull analysis by comparing against mixtures of the competing phases.
