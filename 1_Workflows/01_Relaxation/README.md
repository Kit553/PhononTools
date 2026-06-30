## Relaxation:
**Use**: Obtain the local athermal minimum of the structure

## Workflow
#### 0. Obtain the initial structure guess:

###### 0.1 From Experimental data:
 - Get a cif from Pearson Crystal Database (PCD), Materials Project, Inorganic Crystal Structure Database or experimental data
 - Open the cif in VESTA, then under File > Export > VASP structure file

#### 0.2 From previous calculations:
 - Quantum Materials Open Database directly provides the POSCAR necessary from previous users
 - Using the /Structure/Conversion/POSCAR-ToCif.py tool can convert it back to cif format to inspect it or adjust

#### 1. Setting up the working directory:
 **Your working directory should now contain:**
    - `POSCAR` &rarr; POSition file containing all atom species and their fractional coordinates
    - `POTCAR` &rarr; POTentials used in the simulation, can be obtained using the apgu script in 




















## Parameters:

INCAR:
  SYSTEM = RELAX
  NWRITE = 2
  PREC = Accurate
  KPAR = 1
  NCORE = 4

  ENCUT = 520
  GGA = PS

  LCHARG = .FALSE.
  LWAVE = .FALSE.

  LASPH = .TRUE.
  ADDGRID = .TRUE.
  EDIFF = 1E-08
  ALGO = Normal
  NELM = 100
  NELMIN = 5
  NELMDL = -5

  EDIFFG = -1E-03
  NSW = 100
  IBRION = 2
  ISIF = 3
  IWAVPR = -5

  ISMEAR = -5
  #SIGMA = 0.04
