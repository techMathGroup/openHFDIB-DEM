#!/bin/sh

rootDir=$(pwd)

echo "Working from directory: $rootDir"

# compile the solvers
cd applications/solvers
for pd in ./*/ ; do
    [ -L "${pd%/}" ] && continue
    cd $pd
    for sd in */ ; do
        [ -L "${sd%/}" ] && continue
        cd $sd
        wclean
        wmake
        cd ..
    done
    cd ..
done
cd $rootDir

exit 0