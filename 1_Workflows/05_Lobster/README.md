## Lobster
**Use:** Obtain information on the electronic structure of the material and get the relvant bonding descriptors.

## Workflow:
**Folder Setup:**
**Static Calculation Run**
- `POSCAR`  &rarr;
- `POTCAR`  &rarr;
- `INCAR`   &rarr;
- `KPOINTS` &rarr; 
- `job.sh`  &rarr;

**Lobster Specific**
- `lobsterin` &rarr; Input Parameters for Lobster; see below on details to write this file

## Workflow:



---

#### 1. Creation of the Inpout file:
Creation of the `lobsterin` file depends strongly on the usecase of your calculation. The file provided here is an all-rounder
providing basically all metric needed for thorough bonding analysis, an explanation of the keywords is given below, <opt> marks optional tags.

**Mandatory Keywords**
- `basisSet` &rarr; Defines the electronic basis set used for Lobster; recommended to use `pbeVaspFit2015` for PBE functionals, be sure to adjust this when using other funtionals though
- `basisFuntions` &rarr; defines what funtions are used and what orbitals are considered syntax like `basisFunctions Element orbitals` e.g. `basisFunctions S  3s 3p`

**Recommended Keywords**
- `COHPstartEnergy` and `COHPendEnergy` &rarr; Define the considered energy window in eV; this applies to all metrics
- `COHPsteps` &rarr; Sets the energy grid resolution; this applies to all metrics
- `saveProjectionToFile` &rarr; Saves the projection to a local file for reruns
- `writeBasisSetFunctions` &rarr; Writes out the used basis funtions; important for debugging and reproducibility
- `printTotalSpilling` &rarr; Explicitly prints out charge spilling for diagnostics

· · ·

**Descriptor Specific Keywords**

**Projection Quality:** More explicit reconstruction of the projected basis set; use for reproducibility and validation
- `kpointwiseSpilling`     &rarr; reports the charge spilling per *k*-point
- `bondwiseSpilling`       &rarr; reports the charge spilling per bond type
- `loadProjectionFromFile` &rarr; used to restart desriptor calculation from already existing calculation

. . .

**COHP:** Energy-resolved bonding analysis:
- `cohpGenerator` &rarr; defines the spatial extension in which COHP is generated between atoms syntax like `cohpGenerator from start to end type Element1  type Element2 <opt>:orbitalWise` e.g. `cohpGenerator from 1.8 to 2.4 type P  type S orbitalWise`
- `cohpBetween`   &rarr; requests COHP calculation between an explicit pair of atoms syntax like `cohpBetween atom 1 atom 2 cell n1 n2 n3 <opt>:orbitalWise` where n1 n2 n3 denote the unitcell, use negative for wrapping e.g. cohpBetween atom 11 atom 12 cell -1 1 0 `orbitalWise`
- `skipCOHP`      &rarr; requests to skip COHP calculaiton for a pair of atoms

**COOP:** Overlap population analysis, created alongside the COHPs
- `skipCOHP` &rarr; explicilty skips COOP calculation for a pair of atoms, same syntax as `cohpGenerator`

**COBI:** Bond order analysis, created alongside the COHPs but can be explicitly computed
- `cobiBetween` &rarr; explicitly requests COBI calculations, same syntax as `cohpGenerator`
- `skipCOBI`    &rarr; explicilty skips COBI calculation for a pair of atoms, same syntax as `cohpGenerator`

**DOS/pDOS:** Electronic Density of states, can be element- and/or orbital-resolved, created alongside the COHPs
- `skipDOS` &rarr; explicitly requests to skip DOS/pDOS creation
- `LSODOS`  &rarr; writes the DOS 


. . .

****

---

## Theory and Further Reading:
\[1\]
\[1\]
\[1\]
\[1\]
\[1\]
\[1\]
\[1\]
\[1\]
