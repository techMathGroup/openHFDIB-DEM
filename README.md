# openHFDIB-DEM

Implementation of the Hybrid Fictitious Domain-Immersed Boundary (HFDIB)
coupled with the Discrete Element Method (DEM) by

Martin Isoz         (https://github.com/MartinIsoz)
Martin Šourek       (https://github.com/MSourek)
Ondřej Studeník     (https://github.com/OStudenik)
and
Petr Kočí           (http://vscht.cz//monolith/)

Initial HFDIB implementation spans from

https://github.com/fmuni/openHFDIB

by Federico Municchi

however, the code was completely rewritten.

The DEM implementation for arbitrarily shaped solids into OpenFOAM is
original.

Code capabilities
-----------------
* simulations with either spherical or STL-defined particles (see HFDIBDEM/geomModels)
* simulations of two-phase (solid-fluid) flow via a proof-of-concept
  pimpleHFDIBDEMFoam (see pimpleHFDIBFOAM)
* adaptive mesh refinement based on the particles position
* spring-dashpot contact model based on particle elastic modulus and
  a damping constant (see HFDIBDEM/contactModels)
* extension to contact model for inclusion of adhesive forces
* multiple solid phase initialization options such as random spatial
  distribution of uniformly sized bodies (see HFDIBDEM/addModels)
* instructive tutorials for a single particle falling through a fluid and
  for interaction between a particle and a complex-shaped impeller (see Tutorials)
  
Compatibility
-------------
The code is prepared for compilation with OpenFOAMv8 (https://openfoam.org/version/8/)

Compilation
-----------
Note: the scripts have to be ran from terminals with sourced OpenFOAMv8

* compileAll.sh     -> compiles openHFDIB-DEM and pimpleHFDIBDEM solver
* compileLib.sh     -> compiles openHFDIB-DEM library only
* compileSol.sh     -> compiles pimpleHFDIBDEM solver only 
    (!! requires the library to be compiled !!)


Cite this work as
-----------------

Isoz, M.; Šourek, M.; Studeník, O.; Kočí, P.: Hybrid fictitious domain-immersed boundary solver coupled with discrete element method for simulations of flows laden with arbitrarily-shaped particles. Submitted to Computers & Chemical Engineering, 2021.
