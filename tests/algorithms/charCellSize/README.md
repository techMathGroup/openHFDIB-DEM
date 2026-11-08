# Characteristic cell size lifecycle regression

This suite runs three short pure-DEM cases against the real
`HFDIBDEMFoam` executable. It checks the lifecycle introduced by
`dc6ecd1` without duplicating the line-interpolation regression:

- `exit`: a represented sphere accelerates out of the mesh. Its logged
  size must be 0.01 before it becomes unavailable, and the body must
  switch off once without a stale contact estimate.
- `fallback`: a prescribed, motionless sphere fully contains the only
  mesh cell. With no surface cells but one internal cell, the size must
  be the cube root of the cell volume. Three body recreations must emit
  only one fallback warning.
- `stlWall`: an ordinary sphere collapses with an STL-wall body over
  two steps. The ordinary body supplies the size; the intentionally
  unsized wall must not reach the guarded getter.

The fallback variant forces legacy body creation. The connectivity
path requires a surface seed and cannot represent this degenerate
fixture. A one-cell mesh also isolates cell-size selection from halo
construction, which is outside this regression's scope.

## Run

Source the OpenFOAM v2412 environment and ensure the rebuilt
`HFDIBDEMFoam` is on `PATH`, then run from anywhere:

```sh
python3 -I tests/algorithms/charCellSize/run_regression.py
```

To retain generated cases and logs for diagnosis:

```sh
python3 -I tests/algorithms/charCellSize/run_regression.py \
    --keep /tmp/charCellSize-regression
```

The keep directory must not already exist. Otherwise, temporary cases
are removed in a `finally` block. Every command has a timeout, and logs
larger than 20 MiB fail the suite.

The shared fixture has no saved result times. `stepDEM 1` means one
complete CFD step, rather than an absolute timestep of one second.
