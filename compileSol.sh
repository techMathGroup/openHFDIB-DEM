#!/bin/sh

# compile the solver
cd application/pimpleHFDIBFoam
wclean
wmake
cd ../..

# compile the solver
cd application/pimpleHFDIBFoam
wclean
wmake
cd ../..

exit 0
