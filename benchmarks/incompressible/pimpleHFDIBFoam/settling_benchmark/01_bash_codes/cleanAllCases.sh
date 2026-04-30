#!/bin/bash

module load openfoam/2312

# Source tutorial clean functions
source $WM_PROJECT_DIR/bin/tools/CleanFunctions

for dir in */; do
  if [[ -d "$dir" ]]; then
    echo "Entering directory: $dir"
    
    if [[ -f "$dir/Allclean.sh" ]]; then
      echo "Cleaning simulation data in $dir"
      (
        cd "$dir" || { echo "Failed to enter $dir"; exit 1; }
        cleanCase || echo "Error cleaning the case in $dir"
        rm -rf bodiesInfo || echo "Error removing bodiesInfo in $dir"
      )
    else
      echo "Warning: some other problem encountered"
    fi

  fi
done

module unload openfoam/2312
