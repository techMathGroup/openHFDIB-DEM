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

# directories to skip during compilation
skip_dirs="VoFHFDIB multiphaseInterHFDIBFoam"

# compile the solvers
cd applications/solvers
for pd in ./*/ ; do
    [ -L "${pd%/}" ] && continue
    echo "compiling solvers in $pd"
    cd $pd
    for sd in */ ; do
        [ -L "${sd%/}" ] && continue
        dir_name=$(basename "$sd")
        for skip in $skip_dirs; do
            if [ "$dir_name" = "$skip" ]; then
                continue 2
            fi
        done
        echo "... solver: $dir_name"
        cd $sd
        if [ "$do_clean" = true ]; then
            wclean
        fi
        wmake
        cd ..
    done
    cd ..
done
cd $rootDir

exit 0
