#!/bin/sh

# compile the library
cd HFDIBDEM
wclean
wmake libso
cd ..

exit 0
