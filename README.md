# TBFIT
Tight-binding FITting package (temporary name)

Now you can fit your tight-binding parameters via Slatet-Koster method with a very little effort. 

TBFIT is a scientific Fortran program for numerical tight-binding parameter fitting mainly based on Slater-Koster scheme and tight-binding calculations for the electronic band structures of given atomic and electronic configurations with an simple input interfaces. Basically you can fit Slater-Koster parameters including scaling factors to your target first-principles band structure by TBFIT. You can taste how TBFIT works and its performance in Example section. Once you get proper tight-binding parameters, you can also calculate various Berry phase related quantities, such as Berry curvature, Zak phase, Wilson's loop, first Chern number, and so on. In addition, density of states, eigenstate charge density or wavefunction plot, STM simulations (integerated eigenstate density within certain energy window), band structures for edge/surface geometries, E-field effects, etc.

You do not need to specify all the bond connections, instead, just provide hopping classes (whether it is 1st-, 2nd-, or 3rd-nearest-neighbor hopping) between atomic species in your input file. Then TBFIT automatically set up tight-binding hamiltonian based on your input geometry and tight-binding parameter setup as defined in input interfaces. 

In the future release, I will add some routines for the Green function approach to get the surface state (or edge spectrum as well) and the routines for the spin/mirror Chern number evaluation within the given tight-binding parameter and model hamiltonian.

For the details and examples you can find documents that describes the input tags in /MANUAL/ folder and several input/output files in EXAMPLE/ folder, respectively.

