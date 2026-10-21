# surface-loading law: table check

The regime-aware normal force of `contactModels/surfaceLoadingLaw.H` is a
set of free functions over plain `double`s, so it can be checked without
OpenFOAM, a mesh or a dictionary.

```
g++ -O2 -I../../src/HFDIBDEM/contactModels \
    -o surfaceLoadingLawTable surfaceLoadingLawTable.C
./surfaceLoadingLawTable > table.txt
python3 -I reference.py table.txt
```

`surfaceLoadingLawTable.C` evaluates the law for a kapton flyer on an
aluminium sample at velocities that bracket every regime (elastic, blend
window, full shock, both sides of the kapton melt kink) and prints
`v F pMean pHug pStr pEdge f Z kEff`, plus one legacy row.

`reference.py` recomputes the same quantities from the equations with a
different algorithm — the impedance match is solved in pressure with a
dense scan and regula falsi, where the C++ scans in particle velocity and
bisects — and fails on any relative difference above 1e-9, or above 1e-12
on the legacy row (which must be bit-exact).

`table.txt` is generated output, not a golden file.

The solver-side checks (inertness regression and a physics run) live with
the cases, see `1000_loading_models/ai_reports/surface-loading-model-design.md` §13.
