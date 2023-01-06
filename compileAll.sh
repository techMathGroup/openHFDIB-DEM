#!/bin/sh

# compile the library
cd src/HFDIBDEM
wclean
wmake libso
cd ../..

# compile the solver
cd application/pimpleHFDIBFoam
wclean
wmake
cd ../..

cd application/HFDIBDEMFoam
wclean
wmake
cd ../..

exit 0
