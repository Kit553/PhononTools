# Inventory of the JCH reference collection

Inspected 2026-10-05 at `~/projects/from_JCH`. This inventory is based on archive member listings, reading the scripts and input files, extracting PDF text, reading spreadsheet labels/formulas, and inspecting scanned title/rotation pages. The books are catalogued by title/topic; this is not a claim to have read every book. Original files and archives were left intact. Supplied executables and job scripts were not run, and no contents were uploaded to an external service.

## Directly useful for ordering

### `supercell (2).zip`

Contains the **Supercell** program by Kirill Okhotnikov, Thibault Charpentier and Sylvian Cadars, associated with the 2016 paper “Supercell program: a combinatorial structure-generation approach for the local-level modeling of atomic substitutions and partial occupancies in crystals” (DOI 10.1186/s13321-016-0129-3).

- `supercell/supercell`: supplied executable; its runtime/platform compatibility was not tested. It was not installed or executed.
- `supercell_manual.pdf`: command-line documentation.
- `Supercell program tutorial.pdf`: worked examples of partial occupancies, substitution enumeration, symmetry merging, charge restrictions and selection. The opening example uses Al/Si disorder in Ca2Al2SiO7.
- `art%3A10.1186%2Fs13321-016-0129-3.pdf`: the research article.
- `README`: installation/launch guidance, including OpenBabel's `BABEL_DATADIR`.
- `atomtyp.txt`, `bondtyp.txt`, `element.txt`, `space-groups.txt`, `types.txt`: accompanying chemistry/crystallographic data.
- `ExtractTotalEnergies.py`, attributed to J. M. Skelton: walks numerically named folders, reads `OSZICAR`, and writes `ExtractTotalEnergies.csv`.
- `Xrun.txt`: legacy launcher that copies `XPOSCAR-NNN` and shared VASP input files into numbered folders and invokes `vrun-external`.

Supercell is relevant for exact occupational enumeration and symmetry reduction. Its electrostatic ranking is a heuristic based on assigned charges, not a DFT formation energy or an SQS criterion. CAT's first implementation uses icet for library-backed enumeration/SQS and pymatgen/spglib for structure handling rather than depending on this archived executable.

Two concrete script issues:

1. `ExtractTotalEnergies.py` matches **`F=`**, although its CSV header says `E0 [eV]`. It retains the last match, warns on multiple matches, uses a narrow number regex, and does not validate SCF/ionic convergence, composition, method or normalization. It is useful as a historical sketch, not a reliable hull importer.
2. `Xrun.txt` executes `rm -rf` on each numbered destination before recreating it, and reuses a shared `POTCAR`. Existing jobs could be destroyed, and changed species/orderings need a matching potential sequence. CAT does not use this launcher; its new output directories refuse overwrites.

## Phonons and thermodynamics

### `Computational chemistry - Phonon.pdf`

Five-page exported OneNote practical guide for a VASP/Phonopy workflow. It refers to local shell helpers such as `scell`, `phrun` and `mkForces`, supercell/displacement generation, force collection and thermal/phonon analysis. These aliases are not self-contained software distributed by the PDF. Its example mesh/cell suggestions are starting points, not convergence criteria for arbitrary materials.

### `The physical significance of imaginary phonon.pdf`

Pallikara, Kayastha, Skelton and Whalley, “The physical significance of imaginary phonon modes in crystals,” *Electronic Structure* **4**, 033002 (2022), DOI 10.1088/2516-1075/ac78b3. Review of the interpretation of imaginary modes, numerical artifacts versus structural instabilities, anharmonicity and associated structure prediction. Relevant to assessing the stability of generated orderings, and to the earlier cubic/tetragonal Na3PSe4 question.

### `Helmholtz.xlsx`

Three sheets, with thermal data labels for T, F, S, Cv and E; no worksheet formulas were found. Primarily a stored thermal-data table rather than an automated solver.

### `Helmholtz_v1.xlsx`

Three sheets; Sheet1 has 608 stored formula cells, two parallel thermal datasets, static total energies, 12-atom normalization and a kJ/mol-to-eV divisor of 96.485. Specific formulas show what the comparison actually computes:

- `H8 = G8/$G$1` converts the column labeled **E (kJ/mol)**.
- `I8 = ($D$4+H8)/$F$4` adds that phonon energy to a static cell energy and divides by atom count. The second dataset uses analogous formulas in Q/R.
- `L8 = R8-I8` compares the datasets; B's reference column subtracts I from itself.
- The Helmholtz **F** columns exist separately in C and M and are not used in those displayed E(T) formulas.

Thus the displayed comparison is static-plus-vibrational **internal energy**, not Helmholtz free energy at finite T. A free-energy hull would need a consistent F(T) reference and normalization, plus any required configurational contributions. The “Total E/Unit Cell” label on a quantity divided by atom count is also misleading.

### `thesis_Dengg.pdf`

Thomas Dengg, *Computation of Temperature Dependent Elastic Constants within the Framework of Density-Functional Theory*, PhD thesis (MCL/KFU Graz). Theory and applications cover elasticity, stress/strain, DFT, chemical disorder, phonons, Helmholtz/electronic/vibrational free energies, quasi-harmonic and Debye–Grüneisen models, thermal expansion and temperature-dependent elasticity. Sections 2.3 and 9.3.2 address chemical disorder; the text discusses SQS, virtual crystal and coherent-potential approaches. Useful theoretical background for future finite-temperature alloy/elasticity workflows, not a ready-made ordering API.

## Bonding and electronic structure

### `COHP.zip`

Four files:

- `DFT Calculations.docx`: introductory computing/DFT notes covering SSH and the University of Houston cluster environment, MobaXterm/VESTA and calculation setup. Contains site-specific instructions rather than a portable workflow package.
- `lobsterin`: explicit Sc/Ru/B/Ir basis functions and ten pair-family COHP requests over 1–4 Å, using an energy window from −20 to +10 eV relative to EF.
- `lobstersubmit`: Slurm launcher with an absolute path to the author's LOBSTER binary; needs local adaptation before use.
- `NBANDS.xlsx`: orbital/atom-count worksheet for estimating LOBSTER basis/band counts; 31 formula cells with example systems.

Useful alongside CAT's existing structure-aware LOBSTER preparation, particularly for understanding historical basis and interaction choices. The chemistry and distance windows should not be copied unchanged to Na3PnX4.

### `Lobster_Users_Guide_4.1.0.pdf` and `lobster 4.1.zip`

LOBSTER 4.1.0 manual plus a software/example archive. The archive contains Linux 4.1.0 and macOS 4.1.0 executables, a Windows executable labeled **4.0.1**, a FAQ, changelog, and examples for **VASP, ABINIT and Quantum ESPRESSO**. Example systems include carbon nanotubes, diamond, FeNi3, fullerene, GaAs and Ti; non-VASP examples also contain potential data. The archived user guide duplicates the standalone manual by filename/version. This is a useful historical distribution, not evidence that 4.1.0 is the currently preferred LOBSTER version.

### `The Intermetallic Reactivity Datab.pdf`

Van Buskirk, Kraus and Fredrickson, *The Intermetallic Reactivity Database: Compiling Chemical Pressure and Electronic Metrics toward Materials Design and Discovery*, DOI 10.1021/acs.chemmater.3c00214. Describes a database focused on bonding/DFT chemical-pressure indicators and structural reactivity, alongside electronic metrics. Relevant to future structure-selection descriptors and bonding interpretation, but not an ordering/SQS generator or supplied database dump.

## Elasticity, shear and ideal strength

### `Stress_Strain_Template.xlsx`

Sheets: MnB4, PLQY-Debye, Optimize, Results, Cubic, Hexagonal, Tetragonal, Orthorhombic and Monoclinic. Includes elastic-matrix input layouts, kbar/GPa conversion, Voigt/Reuss/Hill elastic moduli and quasi-harmonic Debye-related fields. Useful as a reference for a future tested elasticity module. The workbook has not been comprehensively validated numerically.

### `stress_strain/`

- `rotations.pdf`: four scanned handwritten pages with plane/direction labels, cube sketches and Euler-angle rotation instructions/checklists. The first page was visually inspected; it is not machine-readable rotation software.
- `image (5).png`: visually inspected slide summarizing crystallographic notation: points, directions `[hkl]`, direction families `<hkl>`, planes `(hkl)` and plane families `{hkl}`. It is a teaching reference, not a dataset or executable.
- `shear  starting files.zip`: legacy Python 2 `strength.py`, VASP `INCAR`/`KPOINTS`, Slurm `job`, `input.dat`, a blank rotations spreadsheet, and many orientation directories at strain steps **0.01, 0.02 and 0.04**.

`strength.py` rotates lattice vectors, applies incremental tensile or shear deformation, executes the command from `input.dat`, extracts an `OUTCAR` stress component, converts kbar to GPa, writes `strength.dat`, archives each step, and copies `CONTCAR` to `POSCAR` for the next increment. The sample command names a nonstandard **`vasp_shear`** executable; its implementation/constraint modifications are not supplied. The accompanying `INCAR` has `ISIF=3`, so simply replacing that executable with standard unconstrained VASP can relax away the imposed strain. The script is Python-2-specific, runs shell commands, overwrites working inputs, and performs no convergence validation. Treat it as a historical method to study, not something to launch directly in a calculation directory.

## `Books.zip`

About 894 MB compressed, with 30 nondirectory members: 29 PDFs plus a nested ZIP. This is a reference library rather than a collection of ordering scripts. Its subjects include:

| Area | Titles/files |
| --- | --- |
| Crystallography | *Crystal Structure Refinement: A Crystallographer's Guide*; *Crystal Structure Analysis—Principles and Practice*; *Introduction to X-Ray Powder Diffraction Data Analysis*. |
| DFT and solid-state bonding | Koch/Holthausen's *A Chemist's Guide to Density Functional Theory*; Dronskowski's *Computational Chemistry of Solid State Materials* (two differently named copies/editions); *Chemical Bonding*; Hoffmann's *Solids and Surfaces*. |
| General solid-state chemistry/synthesis | Frank Hoffmann's *Solid State Chemistry*; Tyagi/Ningthoujam (eds.), *Handbook on Synthesis Strategies for Advanced Materials*, Volume III: Materials Specific Synthesis Strategies (the truncated `2021_Book_HandbookOnSynthesisStrategiesF.pdf`, identified from its title page). |
| Mechanical behavior | *Electronic Basis of the Strength of Materials*; *Mechanical Behavior of Materials*; *The Hardness of Solids*; a fracture-toughness chapter. |
| Poorly named PDFs, identified from title pages | `Chemistry and Physics of2009-05-26 15-39-2554.PDF` is John J. Gilman's *Chemistry and Physics of Mechanical Hardness*. `Adobe Acrobat 7.02010-07-27 11-35-2853.PDF` is Richard W. Hertzberg's *Deformation and Fracture Mechanics of Engineering Materials*, fourth edition. |
| Ceramics | *Handbook of Advanced Ceramics and Composites*; *Handbook of Ceramic Hard Materials*. |
| Rare-earth materials | *Alloy Systems and Compounds Containing Rare Earth Metals and Carbon*; two *Handbook on the Physics and Chemistry of Rare Earths* volumes, one including actinides. |
| Nuclear/radiation materials | *Comprehensive Nuclear Materials*, volumes 1–5; an additional nested `Comprehensive Nuclear Materials 1st Edition.zip`; *Structural Materials for Generation IV Nuclear Reactors*. |
| Machine learning | *Ensemble Machine Learning: Methods and Applications*; *XGBoost with Python*. |

The nested nuclear-materials ZIP was catalogued as an archive; its full contents were not extracted or read.

`Thumbs.db` is a Windows thumbnail cache, with no scientific workflow role.

## Reuse priorities

1. For current ordering work: Supercell paper/tutorials and Dengg's discussion of disorder. Use established maintained libraries for generation and SQS.
2. For finite-temperature extensions: imaginary-phonon review, Dengg thesis, and the thermal spreadsheets after correcting quantity/normalization assumptions.
3. For bonding/descriptor extensions: LOBSTER examples and the IRD paper, alongside CAT's existing bonding tools.
4. For a later elasticity module: workbook formulas and rotation/strength notes, with a new validated strain-constraint implementation and numerical tests.
