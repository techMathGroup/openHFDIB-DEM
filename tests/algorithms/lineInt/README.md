# lineInt processor-handoff regression

Run from the repository root in an OpenFOAM environment containing the
current HFDIB user library:

```bash
python3 -I tests/algorithms/lineInt/run_regression.py --strict
```

The runner uses the real `blockMesh`, `decomposePar`, `reconstructPar`, and
`pimpleHFDIBFoam` applications. It builds deterministic 128-cell cases for
three comparisons:

- an np1/np2 comparison with a processor split at x=2;
- an np1/np4 comparison with a 2x2x1 split and a sphere near the processor
  intersection, exercising targets received on two consecutive ranks;
- the same multi-hop comparison on a bounded, seeded nonorthogonal mesh.

`deform_mesh.py` adapts the shared-interior-point deformation idea from
`nonOrthogonalizeMeshV2.py`. The runner first generates a fresh blockMesh,
then moves its x/y coordinates with seed 17 and an amplitude limit of 0.25
cell spacings per coordinate. Outer boundary points and mesh topology stay
unchanged, coordinates retain full output precision, and the source fixture
is never modified. `checkMesh -allGeometry -allTopology` must report a valid
mesh, and the test measures a nonzero nonorthogonality angle to ensure that
the deformation is actually active. No third-party Python packages are
required.

Each comparison requires:

1. clean solver completion without unresolved-point warnings;
2. nonempty, equal serial/parallel final-owner completion counts;
3. every valid owner-completed point contained in its registered cell;
4. at least one actual handoff, and at least one two-hop handoff in the np4
   case, so a loss of coverage fails instead of silently passing;
5. identical physical surface-cell sets and completed-order histograms;
6. identical final validity and coordinates on matched physical
   cell/initialization/order keys, with tolerance 1e-10 on aligned meshes
   and 1e-8 on the nonorthogonal mesh;
7. identical reconstructed U, p, Ui, f, and lambda production fields,
   ignoring only the generating time-directory location in the header.

Physical cell centers are paired geometrically within 1e-8, not frozen as
expected coordinates or split into unstable rounded dictionary keys.
Local cell IDs and origin labels are not compared across decompositions.
The serial case runs directly because OpenFOAM rejects `mpirun -np 1` with
`-parallel`.

`HFDIB_MARCH_DUMP=1` enables one default-off `marchProbe.<rank>.dat`
stream per rank. The trace contains only the fields consumed by the runner:

- `MARCH_OWNER time order origin_rank origin_label owner_rank owner_cell
  contained`;
- `MARCH_RECV time order origin_rank origin_label`;
- `MARCH_FINAL time order center_x center_y center_z final_rank point_x
  point_y point_z`.

A negative `final_rank` marks an invalidated point. Compare these files,
not interleaved MPI stdout: concurrent writes can corrupt individual
output lines and produce false differences.

Without `--strict`, unresolved points are invalidated and excluded from
later orders; the master emits one warning per run. `--strict` instead
requests a collective abort after the global unresolved count is known.
This normal-path regression treats any such warning as a failure.

Generated cases are deleted automatically. To retain them, supply a new
empty directory:

```bash
python3 -I tests/algorithms/lineInt/run_regression.py \
  --strict --keep /tmp/lineInt-handoff-regression
```
