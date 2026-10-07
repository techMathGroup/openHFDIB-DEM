#!/bin/bash
# clean the wall-contact benchmark case (keeps the mesh-agnostic
# inputs; removes blockMesh output, logs, plot data and time dirs)
cd "$(dirname "$0")/case"

rm -rf constant/polyMesh 0 0.* log.* plots

# ----------------------------------------------------------------- end-of-file
