#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

python3 -m unittest discover -s tests/virtualMesh -p "test_*.py" -v

# solver-level regression tests (need OpenFOAM env: blockMesh, decomposePar,
# HFDIBDEMFoam, mpirun); skipped automatically if the tools are missing
python3 -m unittest discover -s tests/HFDIBDEMFoam -p "test_*.py" -v
