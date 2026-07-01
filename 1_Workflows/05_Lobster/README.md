## Lobster
**Use:** Obtain information on the electronic structure of the material and get the relvant bonding descriptors from a finished static VASP run.

## Workflow:
**Folder Setup:**
**Static Calculation Run**
- `POSCAR`  &rarr; the final relaxed structure of your system
- `POTCAR`  &rarr; the corresponding potential file
- `INCAR`   &rarr; check that the Lobster relevant tags are set to the specifications below
- `KPOINTS` &rarr; 
- `job.sh`  &rarr;

**Lobster Specific**
- `lobsterin` &rarr; Input Parameters for Lobster; see below on details to write this file

## Workflow:

#### 1. Static Calculation
The first step of this workflow is obtaining a high-quality PAW output from VASP, to achieve this make sure to adjust these tags of the `INCAR`:

**Mandatory Tags**
These tags are required for LOBSTER to work on your run
- Write the `WAVECAR` with `LWAVE = .TRUE.`
- Disable Symmetry with `ISYM = 0`
- Set the appropriate amount of calculated bands with `NBANDS = N`, `N` should be
- `NBANDS` needs to be set to include enough states, a general rule of thumb to setting this tag is $NBANDS \gtrapprox \sum_{Ions} n_{\mathrm{Ion}}\cdot n_{\mathrm{Orbitals}}$

**Recommended Tags**:
These are common strategies to achieve the highest possible level of electronic precision.
- Set `LASPH = .TRUE.` to include non-spherical contributions
- Use the high-precision algorithm `PREC = Accurate`
- Set very strict electronic convergence with `EDIFF` e.g. `EDIFF = 1E-8`
- Use the additional *k*-grid with `ADDGRID = .TRUE.`

---

#### 2. Creation of the Inpout file:
Creation of the `lobsterin` file depends strongly on the usecase of your calculation. The file provided here is an all-rounder
providing basically all metric needed for thorough bonding analysis, an explanation of the keywords is given below, <opt> marks optional tags.

**Mandatory Keywords:** Basis set definition
- `basisSet` &rarr; Defines the electronic basis set used for Lobster; recommended to use `pbeVaspFit2015` for PBE functionals, be sure to adjust this when using other funtionals though; non-specified will use a standard LOBSTER set.
- `basisFuntions` &rarr; defines what funtions are used and what orbitals are considered syntax like `basisFunctions Element orbitals` e.g. `basisFunctions S  3s 3p`.
The chosen basis is extremely important for any descriptor you calculate, it should match well with the functional set used for the DFT runs**!**

---

**Strongly Recommended Keywords**  
- `COHPstartEnergy` and `COHPendEnergy` &rarr; define the considered energy window in eV; this applies to all metrics  
- `COHPsteps` &rarr; sets the energy grid resolution important for smooth curves without gaussian smearing; this applies to all metrics  
- `saveProjectionToFile` &rarr; saves the projection to a local file for reruns  
- `writeBasisSetFunctions` &rarr; writes out the used basis funtions; important for debugging and reproducibility  
- `printTotalSpilling` &rarr; explicitly prints out charge spilling for diagnostics  

---  

**Descriptor Specific Keywords**  

**Projection Quality:** More explicit reconstruction of the projected basis set; use for reproducibility and validation.  
- `kpointwiseSpilling`      &rarr; reports the charge spilling per *k*-point  
- `bondwiseSpilling`        &rarr; reports the charge spilling per bond type  
- `loadProjectionFromFile`  &rarr; used to restart desriptor calculation from already existing calculation

Advanced
- `basisRotation`           &rarr; rotates the basis set read from the rerun file
- `autoRotate`              &rarr; allows LOBSTER to rotate the basis set where needed in reruns
- `doNotOrthogonalizeBasis`, `skipReOrthonormalization` and `doNotUseAbsoluteSpilling` &rarr; skip validation of the projected basis 

---  

**COHP:** Energy-resolved bonding analysis, the main descriptor for bonding in LOBSTER.
- `cohpGenerator` &rarr; defines the spatial extension in which COHP is generated between atoms syntax like `cohpGenerator from rmin to rmax type Element1  type Element2 <opt>:orbitalWise` 
- `cohpBetween`   &rarr; requests COHP calculation between an explicit pair of atoms syntax like `cohpBetween atom 1 atom 2 cell n1 n2 n3 <opt>:orbitalWise` where n1 n2 n3 denote the zero-based indices of the unitcell, if a primitive cell is used, wrapping between periodic images can be done by providing negative indices. The atom numeration follows the one-based `POSCAR` syntax.  
- `skipCOHP`      &rarr; requests to skip COHP calculaiton for a pair of atoms  
- `kSpaceCOHP`    &rarr; requests *k*-resolution of the COHP  
The provided script ... can read in your POSCAR and automatically generate the `cohpBetween` lines to match your usecase. 

**COOP:** Overlap population analysis, created alongside the COHPs.  
- `skipCOHP` &rarr; explicilty skips COOP calculation for a pair of atoms, same syntax as `cohpGenerator`  

**COBI:** Bond order analysis, created alongside the COHPs but can be explicitly computed.
- `cobiBetween` &rarr; explicitly requests COBI calculations, same syntax as `cohpGenerator`  
- `skipCOBI`    &rarr; explicilty skips COBI calculation for a pair of atoms, same syntax as `cohpGenerator`  

**DOS/pDOS:** Electronic Density of states, can be element- and/or orbital-resolved, created alongside the COHPs. 
- `skipDOS` &rarr; explicitly requests to skip DOS/pDOS creation  
- `LSODOS`  &rarr; writes the DOS   
  
**Population Analysis:** Wavefunction-based atomic charges and electron population at the atoms, created alongside COHPs.  
- `skipPopulationAnalysis` &rarr; skips *Mulliken* and *Löwdin* atom charge calculation  
- `skipGrossPopulation`    &rarr; skips Population analysis entirely  
Population analysis is needed for estimation of Madelung energies *vide infra*  
  
**Madelung Energies:** Lattice energy and site potentials estimated from population analysis.  
- `EwaldSum`             &rarr; changes scaling parameters for long-range electrostatic interaction syntax like `EwaldSum Param1 Param2`
- `skipMadelungEnergies` &rarr; skips lattice energy calculation  

--- 
 
**Fatbands:** Electronic branch diagram resolved via *k*- element- and orbital-contribution.
- `createFatband` &rarr; requests calculation of the Fatband plot, syntax like `createFatband Element Orbital1 OrbitalN`
This feature is **highly** dependent on the projection that is being used, make sure it is appropriate! Due to their tendency to get extremely crowded these are best used sparingly 

---  

**Distribution Funtions:** Weighted spatial RDFs used to reconstruct bonding structure.
- `BWDF`     &rarr; calculate ond weighted distribution funtion
- `BWDFCOHP` &rarr; colculate COHP weighted distribution funtion
These are especially useful for large/disordered/defective systems where the bonding information is spatially variant.
For small/ordered phases this information is usually redunant with COHP/COBI/COOP

---

**Density of Energies:** COHP-like decomposition of the band-structure energy.
- `densityOfEnergy` &rarr; calculate the density of energies (DOE)
This is very informative in systems where bonding isnt necesarrily pair-wise but multi-centered.

---

**Real-Space Analysis:** matrix-level analysis, visualization tools and advanced debugging.
- `printLCAORealSpaceWavefunction` &rarr; prints the real-space LCAO wavefuntion for a selected *k*-point, synatax like `printLCAORealSpaceWavefunction kpoint 1 coordinates 0 0 0 coordinates 1 1 1 box pointsPerAngstrom 25 bandList 10 11`
- `printPAWRealSpaceWavefunction`  &rarr; prints the real-space PAW wavefuntion for a selected *k*-point, syntax like `printPAWRealSpaceWavefunction kpoint 1 coordinates 0 0 0 coordinates 1 1 1 box pointsPerAngstrom 25 bandList 10 11`
- `writeAtomicOrbitals`            &rarr; writes out what orbitals are used in the local basis
- `gridDensityForPrinting`         &rarr; controlls real-space grid spacing, syntax like `gridBufferForPrinting dist`
- `gridBufferForPrinting`          &rarr; adds a buffer region in angström around printed orbitals to prevent boundary cutoff
- `noFFTforVisualization`          &rarr; enforces direct, (potentally) more accurate visualization instead of FFT-based approach
- `realspaceHamiltonian`           &rarr; wrties out the real-space Hamiltonian matrix, can be used for tight-binding style analysis, syntax like `realspaceHamiltonian layers 2`
- `realspaceOverlap`               &rarr; writes the real-space overlapp matrix, synatx like `realspaceOverlap layers N`
- `writeMatricesToFile`            &rarr; extract all matrices to file

---

**Molecular Fragment Analysis:** deeper analysis of molecular fragments in the material
This requires the molecules to be defined first via:
`molecule atom i atom j atom k ... <opt> cell n1 n2 n3>`
- `printLmosOnAtoms`                     &rarr; prints out localized molecular orbitals on the fragment, syntax like `printLmosOnAtoms atom1 atom2 ...`
- `printLmosOnAtomswriteAtomicDensities` &rarr; *vide supra*, additionally writes out the atom-centered densities
- `printMofeAtomWise`                    &rarr; only prints out the atom-centred densities
- `printMofeMoleculeWise`                &rarr; only prints out the molecule-centred densities
- `skipMOFE`                             &rarr; skips printing out the molecular formation energies
- `skipMolecularOrbitals`                &rarr; skips writing molecular-orbital cube files

#### 3. Running the Lobster Calculation

## Theory and Further Reading:
\[1\]
\[1\]
\[1\]
\[1\]
\[1\]
\[1\]
\[1\]
\[1\]
