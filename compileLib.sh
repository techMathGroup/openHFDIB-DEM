#!/bin/sh

# compile the library
cd scr/HFDIBDEM
wclean
wmake libso
cd ..

exit 0
