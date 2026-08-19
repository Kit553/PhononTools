P## Phono3Py

**Versioning:**

**Use:** Obatin the thrid-order forcce constants of your system to perform phonon lifetime analysis or calculate thermal properties. Properties inlcuded in this tutorial are: (Spectral) Thermal Conductivity (RTA,BTE,Wigner), Grüneisen Parameter, Phonon Lifetime/Linewidth, Joint Density of States (JDOS), Phonon Interaction Strenght, and the Spectral Funtion.

**Disclaimer:** It is extremely improtant for this workflow that the structure you supplied is dynamically stable under DFT conditions i.e. no negative frequencies exist. If that is not the case the obtained thermal conductivity will not be reliable**!**

**Folder Setup:**
Your working directory should contain the following files:
- `POSCAR` &rarr; final relaxed structure of your system
- `POTCAR` &rarr; corresponding potentials file
- 
