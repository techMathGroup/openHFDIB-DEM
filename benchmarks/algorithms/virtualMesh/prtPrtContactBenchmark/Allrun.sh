#!/bin/bash
# prt-prt contact virtual-mesh accuracy benchmark: build, run, gate.
# run from this directory with the OpenFOAM environment sourced.
set -e

cd "$(dirname "$0")"

# build the benchmark utility (no-op when up to date)
wmake >/dev/null

cd case
blockMesh > log.blockMesh 2>&1

vmPrtPrtContactBenchmark > log.vmPrtPrtContactBenchmark 2>&1

cd ..
python3 evalPrtPrtContactBenchmark.py case/log.vmPrtPrtContactBenchmark
