#!/bin/bash

source /opt/openfoam8/etc/bashrc

. $WM_PROJECT_DIR/bin/tools/RunFunctions

rm -rf 0

cp -r 0.org 0

paraFoam -touch
paraFoam -builtin -touch


runApplication blockMesh
runApplication setFields
runApplication decomposePar -force
application=`getApplication`

# runApplication $application
runParallel $application

# runApplication reconstructPar 


# ----------------------------------------------------------------- end-of-file
