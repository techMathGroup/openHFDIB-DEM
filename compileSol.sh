#!/bin/sh

# compile the solver
cd pimpleHFDIBFoam
wclean
wmake
cd ..

exit 0
