## External Tools used in the Workflows
It is recommended to set up a `bin/` folder containing the necessary files on the cluster, which specific funtionalities are used is listed in the repsective workflows.

### ALAMODE

ALAMODE is used for harmonic and anharmonic lattice-dynamics workflows,
including force-constant extraction and anharmonic phonon-property analysis.

- Project: ALAMODE
- Upstream repository: <https://github.com/ttadano/alamode>
- Documentation: <https://alamode.readthedocs.io>
- Local status: linked only; no scripts copied here

---

### Lobster and LobsterPy

These are the implementation of COHP and electronic DOS calculations.

- Project: LobsterPy
- Upstream repository: <https://github.com/JaGeo/LobsterPy>
- Documentation: <https://jageo.github.io/LobsterPy/>

---

### ModeMap

ModeMap is used for mapping potential-energy surfaces along selected phonon
modes, especially nice for soft-mode or anharmonic-mode analysis.

- Project: ModeMap
- Upstream repository: <https://github.com/JMSkelton/ModeMap>
- Main scripts: `ModeMap.py`, `ExtractTotalEnergies.py` and `ModeMap_PostProcess.py`
- Local status: linked only; no scripts copied here

---

### PhonoPy

PhonoPy is used as Frozen Phonon method implementation, it is the base for any
harmonic phonon calculaiton.

- Project: PhonoPy
- Upstream repository: <https://github.com/phonopy/phonopy/tree/master>
- Documentation: <https://phonopy.github.io/phonopy/>
- Most of the workflows and analysis methods here are wrappers for this package

---

### Phono3Py

Phono3Py is the Frozen Phonon implementation up to third order force constants,
wrappers for workflow and analysis are provided here.

- Project: Phono3Py
- Upstream repository: <https://github.com/phonopy/phono3py>
- Documentation: <https://phonopy.github.io/phono3py/>

---

