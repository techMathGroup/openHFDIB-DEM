#!/usr/bin/env python3
"""Structural contract checks for the WC2 commit (plane source for
the wall-path exact evaluation).

Checks (local-only, tests/ is gitignored):
1. geomModel declares virtual getBoxSurfacePlane(const boundBox&,
   point&, vector&); base returns false SILENTLY (no InfoH print).
2. stlBased declares the override (final).
3. stlBased::getLeafSubVolumePlane now delegates to the shared
   getShapesSurfacePlane helper (no duplicated facet logic).
4. stlBased::getBoxSurfacePlane routes through the private
   boxOverlappingFacets broad phase (facet binning when valid
   and current, tree.findBox otherwise) - no shapesIn_
   precondition.
5. the helper keeps the validated dominant-facet + probe logic
   (bbox-overlap metric, degenerate-facet guard, both-side probes,
   posInside == negInside fallback).
6. the probe offset uses the box average dimension (leaf.avgDim()),
   matching the prt-path scale (sv.mag() was the subVolume edge).
"""
import re
import sys
import pathlib

root = pathlib.Path(__file__).resolve().parents[3]

failures = []

def check(name, cond):
    print(f"  {name}: {'OK' if cond else 'FAIL'}")
    if not cond:
        failures.append(name)

gm = (root / "src/HFDIBDEM/geomModels/geomModel.H").read_text()
m = re.search(
    r"virtual bool getBoxSurfacePlane\s*\((.*?)\)\s*\{(.*?)\n\s*\}",
    gm, re.S)
check("geomModel declares virtual getBoxSurfacePlane", m is not None)
if m:
    body = m.group(0)
    check("base takes (const boundBox&, point&, vector&)",
          "const boundBox& leaf" in m.group(0))
    check("base returns false silently (no InfoH)",
          "return false" in body and "InfoH" not in body)

sb_h = (root / "src/HFDIBDEM/geomModels/stlBased/stlBased.H").read_text()
check("stlBased overrides getBoxSurfacePlane (final)",
      re.search(r"bool getBoxSurfacePlane\s*\([^)]*boundBox&[^)]*\)\s*final override",
                sb_h) is not None)
check("stlBased declares private getShapesSurfacePlane",
      re.search(r"bool getShapesSurfacePlane\s*\(", sb_h) is not None)

sb_c = (root / "src/HFDIBDEM/geomModels/stlBased/stlBased.C").read_text()

leaf_fn = re.search(
    r"bool stlBased::getLeafSubVolumePlane.*?^}", sb_c, re.S | re.M)
check("getLeafSubVolumePlane body present", leaf_fn is not None)
if leaf_fn:
    b = leaf_fn.group(0)
    check("getLeafSubVolumePlane delegates to helper",
          "getShapesSurfacePlane" in b)
    check("no duplicated facet loop in getLeafSubVolumePlane",
          "forAll(shapesIn" not in b)

box_fn = re.search(
    r"bool stlBased::getBoxSurfacePlane.*?^}", sb_c, re.S | re.M)
check("getBoxSurfacePlane body present", box_fn is not None)
if box_fn:
    b = box_fn.group(0)
    check("getBoxSurfacePlane routes through boxOverlappingFacets",
          "boxOverlappingFacets(leaf)" in b)
    check("getBoxSurfacePlane delegates to helper",
          "getShapesSurfacePlane" in b)
    check("no direct findBox left in getBoxSurfacePlane",
          "findBox" not in b)

# the box-overlap broad phase: binning when its build is current
# and valid, tree.findBox otherwise (the degenerate fallback)
broad_fn = re.search(
    r"labelList stlBased::boxOverlappingFacets.*?^}", sb_c, re.S | re.M)
check("boxOverlappingFacets body present", broad_fn is not None)
if broad_fn:
    b = broad_fn.group(0)
    check("boxOverlappingFacets lazily builds/rebuilds the binning",
          "binningVersion_ != surfaceVersion_" in b
          and "new facetBinning" in b)
    check("boxOverlappingFacets filters via the octree shapes predicate",
          "tree().shapes()" in b and "facetsNear" in b)
    check("boxOverlappingFacets falls back to tree.findBox",
          "findBox(treeBoundBox(queryBox))" in b)

helper = re.search(
    r"bool stlBased::getShapesSurfacePlane.*?^}", sb_c, re.S | re.M)
check("helper body present", helper is not None)
if helper:
    b = helper.group(0)
    check("dominant-facet bbox-overlap metric kept",
          "overlap.volume() > bestArea" in b)
    check("degenerate-facet guard kept", "aNmag < VSMALL" in b)
    check("both-side pointInside probes kept",
          "posInside" in b and "negInside" in b)
    check("non-separating fallback kept",
          "posInside == negInside" in b)
    check("probe offset from box average dimension",
          "avgDim()/8.0" in b)

# 6: prt path unchanged semantically - getLeafSubVolumePlane still
# guards on shapesIn_ validity before delegating
if leaf_fn:
    check("shapesIn_ validity guard kept in getLeafSubVolumePlane",
          "shapesIn_.valid()" in leaf_fn.group(0))

# 7: structural binning invalidation - every mutator that resets
# triSurfSearch_ must bump the surface version (setBodyPosition,
# initVars, synchronPos, bodyMovePoints, bodyRotatePoints,
# bodyScalePoints). the constructor is the one reset without a
# touch: it initializes the counters instead (surfaceVersion_ 0,
# binningVersion_ -1, forcing a first-query build)
n_resets = len(re.findall(r"triSurfSearch_\.reset", sb_c))
n_touches = len(re.findall(r"touchSurface\(\)", sb_c))
check("every non-ctor triSurfSearch_ reset is followed by touchSurface()",
      n_resets == n_touches + 1 and n_resets > 1)
check("ctor initializes the version counters stale",
      re.search(r"surfaceVersion_\(0\)", sb_c) is not None
      and re.search(r"binningVersion_\(-1\)", sb_c) is not None)

sb_h2 = sb_h
initvars_fn = re.search(r"void initVars\s*\([^)]*\)\s*\{.*?\}", sb_h2, re.S)
check("initVars (header-inline) bumps the surface version",
      initvars_fn is not None and "touchSurface()" in initvars_fn.group(0))

print()
if failures:
    print(f"test_wc2_contracts: FAILED ({len(failures)} checks)")
    sys.exit(1)
print("test_wc2_contracts: all OK")
