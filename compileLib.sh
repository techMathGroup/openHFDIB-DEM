#!/bin/sh

rootDir=$(pwd)

echo "Working from directory: $rootDir"

# compile the library
cd src/HFDIBDEM
wclean
wmake libso
cd $rootDir

exit 0