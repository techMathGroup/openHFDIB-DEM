#!/bin/sh

# compile the library
cd HFDIBDEM
wclean
wmake libso
cd ..

# compile the solver
cd pimpleHFDIBFoam
wclean
wmake
cd ..

exit 0
