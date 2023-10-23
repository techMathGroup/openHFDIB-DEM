#!/bin/sh

# compile the library
cd src/HFDIBDEM
wclean
wmake libso
cd ../..

# compile the solver
cd application/solvers/pimpleHFDIBFoam
wclean
wmake
cd ../..

cd application/solvers/HFDIBDEMFoam
wclean
wmake
cd ../..

exit 0
