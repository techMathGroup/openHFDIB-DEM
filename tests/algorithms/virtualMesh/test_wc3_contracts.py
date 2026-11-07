#!/usr/bin/env python3
"""Structural contract checks for the WC3 commit (exact contact
volume in virtualMeshWall::evaluateContact).

Checks (local-only, tests/ is gitignored):
1. evaluateContact collects inside-leaf indices during the flood and
   keeps the legacy count return when the exact guard is off.
2. the exact second pass runs only when the guard is on AND wall
   planes were set (plane virtual meshes stay legacy).
3. band rule: an inside leaf with any non-inside face neighbour, or
   fewer than 6 face neighbours (out-of-matrix seam counts as
   outside), is on the contact band.
4. band leaves go through evaluateContactLeafExact; interior leaves
   keep their full sub-volume.
5. evaluateContactLeafExact: caches via clippedVolume_ sentinel
   (-1 = unevaluated), leaf box from centre +/- half the sub-volume
   edge, getBoxSurfacePlane (WC2), list-variant planePolyClip with
   wall planes + body plane, full-sub-volume fallback when no plane.
6. wallContact.C passes the sub-contact wall half-spaces into the
   contact VMs via setWallPlanes; the half-space keeps the stored
   patch normal (points into the fluid, the planePolyClip kept
   side).
7. the aggregation in getWallContactVars_ArbShape is volume-weighted
   (contactCenterSum accumulated with volumes, divided by
   intersectVolume).
8. planePolyClip list variant: declaration, implementation, and
   default-constructible halfSpace (needed by List<halfSpace>).
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

# --- planePolyClip list variant --------------------------------------
pp_h = (root / "src/HFDIBDEM/virtualMesh/planePolyClip.H").read_text()
check("halfSpace has a default constructor",
      re.search(r"halfSpace\(\)", pp_h) is not None)
check("list variant declared (List<halfSpace> planes)",
      re.search(r"const List<halfSpace>& planes", pp_h) is not None)
check("list core declared (clipAndIntegrateList)",
      "clipAndIntegrateList" in pp_h)

pp_c = (root / "src/HFDIBDEM/virtualMesh/planePolyClip.C").read_text()
m = re.search(r"void planePolyClip::clipAndIntegrateList.*?^}",
              pp_c, re.S | re.M)
check("list core implemented", m is not None)
if m:
    b = m.group(0)
    check("list core loops over all planes",
          "forAll(planes, plI)" in b)
    check("list core keeps the empty-clip early return",
          "everything cut away" in b)
    check("list core integrates from the box midpoint anchor",
          "box.midpoint()" in b)

# --- virtualMeshWall ---------------------------------------------------
vmw_h = (root / "src/HFDIBDEM/virtualMesh/virtualMeshWall.H").read_text()
check("wallPlanes_ member declared",
      "List<planePolyClip::halfSpace> wallPlanes_" in vmw_h)
check("setWallPlanes setter declared", "setWallPlanes" in vmw_h)
check("evaluateContactLeafExact declared",
      "evaluateContactLeafExact" in vmw_h)
check("virtualMeshWall.H includes planePolyClip.H",
      "#include \"planePolyClip.H\"" in vmw_h)

vmw_c = (root / "src/HFDIBDEM/virtualMesh/virtualMeshWall.C").read_text()

ev = re.search(r"scalar virtualMeshWall::evaluateContact\(\).*?^}",
               vmw_c, re.S | re.M)
check("evaluateContact body present", ev is not None)
if ev:
    b = ev.group(0)
    check("inside leaves collected during the flood",
          "insideLeaves_.append" in b)
    check("legacy count return kept",
          "return volumeCount*bbMatrix_.getSubVolumeV()" in b)
    check("exact pass gated on getExactSubVolume()",
          "virtualMeshLevel::getExactSubVolume()" in b)
    check("exact pass gated on wallPlanes_.size() > 0",
          "wallPlanes_.size() > 0" in b)
    check("exact pass gated on volumeCount > 0",
          re.search(r"volumeCount > 0\s*&&", b) is not None)

    band = re.search(r"if \(volumeCount > 0.*?return volumeCount\*bbMatrix_",
                     b, re.S)
    check("exact pass body present", band is not None)
    if band:
        bb = band.group(0)
        check("band test shared via isOnContactBand",
              "isOnContactBand" in bb)
        check("band leaves use evaluateContactLeafExact",
              "evaluateContactLeafExact" in bb)
        check("interior leaves clipped by the wall half-spaces",
              re.search(r"else\s*\{[^}]*planePolyClip::volumeAndCentroid"
                        r"[^}]*wallPlanes_", bb, re.S) is not None)
        check("contact centre volume-weighted in the exact pass",
              "exactCenter/exactVolume" in bb)
        check("detected-but-clipped-to-zero falls back to the count",
              "clipped-to-zero" in bb and "WarningInFunction" in bb)

band_fn = re.search(r"bool virtualMeshWall::isOnContactBand.*?^}",
                    vmw_c, re.S | re.M)
check("isOnContactBand body present", band_fn is not None)
if band_fn:
    bb = band_fn.group(0)
    check("band = facet-overlap leaf (boxIntersectsSurface)",
          "boxIntersectsSurface" in bb)
    check("band test reads the leaf box",
          "leafBox" in bb)

leaf_fn = re.search(r"boundBox virtualMeshWall::leafBox.*?^}",
                    vmw_c, re.S | re.M)
check("leafBox body present", leaf_fn is not None)
if leaf_fn:
    check("leaf box from centre +/- half sub-volume edge",
          "0.5*svEdge" in leaf_fn.group(0))

leaf = re.search(r"void virtualMeshWall::evaluateContactLeafExact.*?^}",
                 vmw_c, re.S | re.M)
check("evaluateContactLeafExact body present", leaf is not None)
if leaf:
    b = leaf.group(0)
    check("cache sentinel guard (clippedVolume_ > -0.5)",
          "clippedVolume_ > -0.5" in b)
    check("leaf box built by the shared leafBox helper",
          "leafBox(subVolume)" in b)
    check("body plane from getBoxSurfacePlane",
          "getBoxSurfacePlane" in b)
    check("list-variant clip with wall + body planes",
          "clipPlanes_" in b and "volumeAndCentroid" in b)
    check("full-sub-volume fallback when no plane",
          "getSubVolumeV()" in b)

# --- wallContact.C -----------------------------------------------------
wc = (root / "src/HFDIBDEM/contactModels/wallContact.C").read_text()
check("contact model builds wall half-spaces from contactPatches",
      re.search(r"forAll\(contactPatches,cP\).*?wallPlanes\.append",
                wc, re.S) is not None)
check("setWallPlanes called on contact VMs",
      "virtMeshWall.setWallPlanes(wallPlanes)" in wc)
check("half-space keeps the stored patch normal (fluid side)",
      re.search(r"halfSpace\(planeInfo\[1\], planeInfo\[0\]\)", wc)
      is not None)
check("volume-weighted aggregation (contactCenterSum)",
      "contactCenterSum" in wc)
check("unweighted mean removed (no contactCenters list)",
      "contactCenters" not in wc)
check("plane VMs get no wall planes (legacy)",
      re.search(r"virtMeshPlane.*setWallPlanes", wc, re.S) is None)

print()
if failures:
    print(f"test_wc3_contracts: FAILED ({len(failures)} checks)")
    sys.exit(1)
print("test_wc3_contracts: all OK")
