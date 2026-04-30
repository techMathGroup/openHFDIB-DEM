#!/bin/sh

if [ "$1" = "clean" ]; then
    do_clean=true
else
    do_clean=false
fi

rootDir=$(pwd)

echo "Working from directory: $rootDir"

# compile the library
echo "compiling the library"
cd src/HFDIBDEM
if [ "$do_clean" = true ]; then
    wclean
fi
wmake libso
cd $rootDir

exit 0
