#!/bin/sh

# compile the solver
cd application/solvers/pimpleHFDIBFoam
wclean
wmake
cd ../..

# compile the solver
cd application/solvers/pimpleHFDIBFoam
wclean
wmake
cd ../..

exit 0
