#!/usr/bin/env python3
"""Structural contract checks for the WC4 commit (exact wetted area
on the plane virtual meshes).

Checks (local-only, tests/ is gitignored):
1. planePolyClip::polygonArea promoted to the library (shoelace,
   centroid-anchored fan, magnitude of the area vector, < 3
   points -> 0).
2. virtualMeshWall::evaluateContact persists the inside-leaf
   indices (insideLeaves_) for the area pass; truncation clears
   them.
3. evaluateContactAreaExact: guard on getExactSubVolume and a
   non-empty flood; band leaves clipped with the wall plane as
   hs1 and the body plane as hs2; interior leaves and no-plane
   leaves keep the full face area; no kept material behind the
   wall -> no footprint; clipped-to-zero falls back to the
   count-based area with a warning.
4. wallContact.C: the plane-VM loop passes the BEHIND-wall half-
   space (stored normal flipped); the exact area multiplies by
   emptyScale; the count-based (V/svVol)^(2/3) fallback is
   retained on -1.
5. band rule and leaf box shared with WC3 via isOnContactBand /
   leafBox (no duplicated neighbor walk).
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

# --- planePolyClip.polygonArea ---------------------------------------
pp_h = (root / "src/HFDIBDEM/virtualMesh/planePolyClip.H").read_text()
check("polygonArea declared in planePolyClip",
      "static scalar polygonArea" in pp_h)

pp_c = (root / "src/HFDIBDEM/virtualMesh/planePolyClip.C").read_text()
pa = re.search(r"scalar planePolyClip::polygonArea.*?^}", pp_c, re.S | re.M)
check("polygonArea implemented", pa is not None)
if pa:
    b = pa.group(0)
    check("polygonArea guards < 3 vertices",
          "poly.size() < 3" in b)
    check("polygonArea fans from the polygon centroid",
          "c /= scalar(poly.size())" in b)
    check("polygonArea returns the area-vector magnitude",
          "mag(areaVec)" in b)

# --- virtualMeshWall ---------------------------------------------------
vmw_h = (root / "src/HFDIBDEM/virtualMesh/virtualMeshWall.H").read_text()
check("insideLeaves_ member declared", "insideLeaves_" in vmw_h)
check("evaluateContactAreaExact declared",
      "evaluateContactAreaExact" in vmw_h)
check("isOnContactBand declared (shared band rule)",
      "isOnContactBand" in vmw_h)
check("leafBox declared (shared leaf geometry)",
      "leafBox" in vmw_h)

vmw_c = (root / "src/HFDIBDEM/virtualMesh/virtualMeshWall.C").read_text()

ev = re.search(r"scalar virtualMeshWall::evaluateContact\(\).*?^}",
               vmw_c, re.S | re.M)
check("evaluateContact persists inside leaves",
      ev is not None and "insideLeaves_.append" in ev.group(0))
check("evaluateContact exact pass reads insideLeaves_",
      ev is not None and "insideLeaves_,lI" in ev.group(0))
check("truncation clears insideLeaves_",
      ev is not None
      and re.search(r"insideLeaves_\.clear\(\);\s*break;", ev.group(0))
      is not None)

# the shared helpers have no duplicated neighbor walk
band_fn = re.search(r"bool virtualMeshWall::isOnContactBand.*?^}",
                    vmw_c, re.S | re.M)
check("isOnContactBand implemented", band_fn is not None)
if band_fn:
    b = band_fn.group(0)
    check("band rule: facet-overlap via boxIntersectsSurface",
          "boxIntersectsSurface" in b)
    check("band rule: no neighbour walk left in the band test",
          "nI <" not in b and "isCBody" not in b)

leaf_fn = re.search(r"boundBox virtualMeshWall::leafBox.*?^}",
                    vmw_c, re.S | re.M)
check("leafBox implemented", leaf_fn is not None)

area = re.search(r"scalar virtualMeshWall::evaluateContactAreaExact.*?^}",
                 vmw_c, re.S | re.M)
check("evaluateContactAreaExact implemented", area is not None)
if area:
    b = area.group(0)
    check("guard on getExactSubVolume",
          "virtualMeshLevel::getExactSubVolume()" in b)
    check("guard on non-empty flood (insideLeaves_.size() == 0)",
          "insideLeaves_.size() == 0" in b)
    check("reuses the persisted flood (no re-flood)",
          "while" not in b)
    check("wall plane as hs1 with volumeCentroidAndFace",
          re.search(r"volumeCentroidAndFace\s*\(\s*leaf,\s*wallPlane", b)
          is not None)
    check("body plane as hs2",
          "getBoxSurfacePlane" in b)
    check("interior leaves keep the full face area",
          "leafFaceArea" in b)
    check("clipped-to-zero falls back to the count-based area",
          "clipped-to-zero" in b and "WarningInFunction" in b)

# --- wallContact.C -----------------------------------------------------
wc = (root / "src/HFDIBDEM/contactModels/wallContact.C").read_text()
# the plane-VM loop: from its header to the aggregation comment
# that follows the loop (inner braces make brace-matching
# unreliable)
plane_loop = re.search(
    r"for\(label i = 0; i< vMPlaneInfoSize; i\+\+\).*?"
    r"volume-weighted mean over every contribution",
    wc, re.S
)
check("plane-VM loop found", plane_loop is not None)
if plane_loop:
    b = plane_loop.group(0)
    check("plane VM paired with contactPatches[i]",
          "i < contactPatches.size()" in b)
    check("behind-wall half-space (stored normal, into the wall)",
          # the stored wall normal points out of the fluid, into
          # the wall - which is the planePolyClip kept-side
          # convention as-is; flipping it (the original form)
          # kept the fluid side and was wrong
          "planeInfo[1],\n                            planeInfo[0]" in b
          and "-planeInfo[0]" not in b)
    check("exact area multiplies by emptyScale",
          "exactArea*vmWInfo->getEmptyScale()" in b)
    check("count-based fallback retained on -1",
          "(vmVolume/vmWInfo->getSVVolume())" in b)
    check("detectFirstFaceContactPoint branch retained",
          "detectFirstFaceContactPoint" in b)

print()
if failures:
    print(f"test_wc4_contracts: FAILED ({len(failures)} checks)")
    sys.exit(1)
print("test_wc4_contracts: all OK")
