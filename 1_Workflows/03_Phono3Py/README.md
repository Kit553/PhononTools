## Phono3Py

**Versioning:** This tutorial is build around the cluster version of `Phono3Py` (v2.7.0) and `VASP` (6.2.0) to generate the third-order force constants, the subsequent analysis is then build around a local installation of `Phono3Py` (v4.3.3).

**Use:** Obatin the thrid-order forcce constants of your system to perform phonon lifetime analysis or calculate thermal properties. Properties inlcuded in this tutorial are: (Spectral) Thermal Conductivity (RTA,BTE,Wigner), Grüneisen Parameter, Phonon Lifetime/Linewidth, Joint Density of States (JDOS), Phonon Interaction Strenght, and the Spectral Funtion.

**Disclaimer:** It is extremely important for this workflow that the structure you supplied is dynamically stable under DFT conditions i.e. no negative frequencies exist. If that is not the case the obtained thermal conductivity will not be reliable**!**

**Folder Setup:**
Your working directory should contain the following files:
- `POSCAR` &rarr; final relaxed structure of your system
- `POTCAR` &rarr; corresponding potentials file
- 
-
