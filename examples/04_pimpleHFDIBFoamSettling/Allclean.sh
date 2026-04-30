#!/bin/sh
source /opt/openfoam8/etc/bashrc

cd ${0%/*} || exit 1    # run from this directory

# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

cleanCase

rm -rf bodiesInfo

# ----------------------------------------------------------------- end-of-file
