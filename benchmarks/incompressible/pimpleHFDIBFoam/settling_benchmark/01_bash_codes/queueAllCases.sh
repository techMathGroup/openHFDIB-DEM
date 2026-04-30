#!/bin/bash
for dir in */; do
  if [[ -d "$dir" ]]; then
    echo "Entering directory: $dir"
    
    if [[ -f "$dir/Allclean.sh" ]]; then
      echo "Queuing the case $dir"
      (
        cd "$dir" || { echo "Failed to enter $dir"; exit 1; }
        sbatch Allrun-slurm || echo "Error queuing the case $dir"
      )
    else
      echo "Warning: some other problem encountered"
    fi

  fi
done
