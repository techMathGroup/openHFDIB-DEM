# !!! this is not a production version of the code, assume that everything is broken and use the official repository: https://github.com/techMathGroup/openHFDIB-DEM !!!

# openHFDIB-DEM
Preparation for porting openHFDIB-DEM to openfoam.com

About openHFDIB-DEM
-------------------
openHFDIB-DEM is a free, open-source CFD-DEM library based on OpenFOAM (https://openfoam.org) and capable of simulating flows laden with arbitrarily-shaped solids. The code is developed mostly by members of the techMathGroup of the Institute of Thermomechanics of the Czech Academy of Sciences (https://www.it.cas.cz/) and members of the Monolith group of the University of Chemistry and Technology, Prague (https://monolith.vscht.cz).

The main contributors are:
* Martin Isoz         (https://github.com/MartinIsoz)
* Martin Kotouč Šourek       (https://github.com/MartinKotoucSourek)
* Ondřej Studeník     (https://github.com/OStudenik)
* Petr Kočí           (https://monolith.vscht.cz)

The implementation itself is based on the Hybrid Fictitious Domain-Immersed Boundary (HFDIB)
coupled with the Discrete Element Method (DEM). The initial HFDIB implementation spans from the work of Federico Municchi (https://github.com/fmuni/openHFDIB). However, the code was heavily modified. The DEM implementation for arbitrarily shaped solids into OpenFOAM is original.

![alt text](https://github.com/techMathGroup/openHFDIB-DEM/blob/main/Images/openHFDIB-DEM_IntroImage.png?raw=true)

Code capabilities
-----------------
* simulation settings are defined in HFDIBDEMdict (see DOCUMENTATION)
* simulations with either spherical or STL-defined particles (see HFDIBDEM/geomModels)
* simulations of two-phase (solid-fluid) flow (see pimpleHFDIBFoam)
* simulations of solid phase in standard DEM mode (see HFDIBDEMFoam)
* adaptive mesh refinement based on the particles position
* spring-dashpot contact model based on particle elastic modulus and
  a damping based on coefficient of restitution(see HFDIBDEM/contactModels)
* extension to contact model for inclusion of adhesive forces
* multiple solid phase initialization options such as random spatial
  distribution of uniformly sized bodies (see HFDIBDEM/addModels)
* instructive tutorials for a single particle falling through a fluid and
  for interaction between a particle and a complex-shaped impeller (see Tutorials)
* possible applications added to examples; see for settings of large-scale simulations
* fully parallel version

Compatibility
-------------
The code is prepared for compilation with OpenFOAMv8 (https://openfoam.org/version/8/)

Compilation
-----------
Note: the scripts have to be ran from terminals with sourced OpenFOAMv8

* compileAll.sh     -> compiles openHFDIB-DEM and pimpleHFDIBDEM solver
* compileLib.sh     -> compiles openHFDIB-DEM library only
* compileSol.sh     -> compiles pimpleHFDIBDEM solver only (!! requires the library to be compiled !!)

For users
---------
Please note that the code is developed by a small research team. The openHFDIB-DEM library is distributed, similarly to OpenFOAM itself, in the hope that it will be useful but without any warranty. Furthermore, we are glad for any remarks and comments on any potential bugs. Please, report them to us via the "Issues" (https://github.com/techMathGroup/openHFDIB-DEM/issues).

For any questions regarding specific cases or code capabilities, we kindly ask you to use the "Discussions" page (https://github.com/techMathGroup/openHFDIB-DEM/discussions). We are trying to be present there and to communicate with you as users. However, due to our team size and workload, we cannot provide individual support to each openHFDIB-DEM user.

Finally, doxygen documentation for the library is available at http://docs.isoz.eu.

For the installation guide and detailed tutorials, see our "Wiki" page (https://github.com/techMathGroup/openHFDIB-DEM/wiki/)

If you need to contact the authors in matters regarding openHFDIB-DEM, please do so via email: openhfdib-dem@it.cas.cz


Cite this work as
-----------------
* Isoz, M.; Kotouč Šourek, M.; Studeník, O.; Kočí, P.: Hybrid fictitious domain-immersed boundary solver coupled with discrete element method for simulations of flows laden with arbitrarily-shaped particles. Computers & Fluids, 244 (2022) 105538, DOI: 10.1016/j.compfluid.2022.105538
* Kotouč Šourek, M.; Studeník, O.; Isoz, M.; Kočí, P.; York, A.P.E.: Viscosity prediction for dense suspensions of non-spherical particles based on CFD-DEM simulations. Powder Technology, 444 (2024) 120067, DOI:10.1016/j.powtec.2024.120067
* Studeník, O.; Isoz, M.; Kotouč Šourek, M.; Kočí, P.: OpenHFDIB-DEM: An extension to OpenFOAM for CFD-DEM simulations with arbitrary particle shapes. SoftwareX, 27 (2024), DOI: 10.1016/j.softx.2024.101871

License
-------
openHFDIB-DEM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software  Foundation, either version 3 of the License, or (at your option) any later version.  See http://www.gnu.org/licenses/, for a description of the GNU General Public License terms under which you can copy the files.
