#!/bin/sh

# compile the library
cd src/HFDIBDEM
wclean
wmake libso
cd ..

exit 0
