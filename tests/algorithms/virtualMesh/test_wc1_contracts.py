#!/usr/bin/env python3
"""Structural contract checks for the WC1 groundwork commit
(wall-path exact volume/area extension).

Checks (local-only, tests/ is gitignored):
1. planePolyClip declares volumeCentroidAndFace (public) and
   clipAndIntegrate (private core) with a facePolygon / faceOut arg.
2. volumeAndCentroid is a thin wrapper: delegates to
   clipAndIntegrate with nullptr, no duplicated clip logic.
3. the hs1 closing face is captured before later half-spaces clip
   it, and a second-half-space re-clip of the captured ring exists.
4. subVolumeProperties carries clippedVolume_ / clippedCentroid_
   with the -1 sentinel, resetClip() called from the ctor.
5. resetSubVolume also calls resetClip().
6. no existing call site of volumeAndCentroid is modified (legacy
   prt path unchanged).
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

# 1-3: planePolyClip declarations and wrapper structure
pp_h = (root / "src/HFDIBDEM/virtualMesh/planePolyClip.H").read_text()
check("volumeCentroidAndFace declared public",
      re.search(r"static void volumeCentroidAndFace", pp_h) is not None)
check("clipAndIntegrate declared private",
      re.search(r"private:.*static void clipAndIntegrate", pp_h, re.S) is not None)
check("facePolygon arg present",
      "DynamicList<point>& facePolygon" in pp_h)
check("faceOut pointer arg present",
      "DynamicList<point>* faceOut" in pp_h)

pp_c = (root / "src/HFDIBDEM/virtualMesh/planePolyClip.C").read_text()
# 2: wrapper contains no clip logic of its own
wrapper = re.search(
    r"void planePolyClip::volumeAndCentroid.*?^}", pp_c, re.S | re.M)
check("volumeAndCentroid delegates to clipAndIntegrate",
      wrapper is not None and "clipAndIntegrate" in wrapper.group(0)
      and "closing" not in wrapper.group(0))
wrapper_face = re.search(
    r"void planePolyClip::volumeCentroidAndFace.*?^}", pp_c, re.S | re.M)
check("volumeCentroidAndFace delegates to clipAndIntegrate",
      wrapper_face is not None and "clipAndIntegrate" in wrapper_face.group(0)
      and "closing" not in wrapper_face.group(0))

# 3: capture-then-reclip structure in the core
core = re.search(
    r"void planePolyClip::clipAndIntegrate.*?\n}\n", pp_c, re.S)
check("clipAndIntegrate body present", core is not None)
core_b = core.group(0) if core else ""
check("face captured at pass 0",
      re.search(r"pass == 0 && faceOut", core_b) is not None)
check("second-half-space re-clip of the face exists",
      re.search(r"faceOut->size\(\) > 0 && pass == 0 && hasSecond", core_b)
      is not None)
check("face cleared when no face requested or produced",
      "faceOut->clear()" in core_b)

# 4: subVolumeProperties cache fields
svp = (root / "src/HFDIBDEM/virtualMesh/subVolumeProperties.H").read_text()
check("clippedVolume_ field with -1 sentinel",
      re.search(r"scalar clippedVolume_ = -1;", svp) is not None)
check("clippedCentroid_ field",
      re.search(r"vector clippedCentroid_ = vector::zero;", svp) is not None)
check("resetClip declared and defined inline",
      re.search(r"void resetClip\(\)", svp) is not None)
ctor = re.search(r"subVolumeProperties\(point center\):(.*?)\}", svp, re.S)
check("ctor calls resetClip",
      ctor is not None and "resetClip()" in ctor.group(0))

# 5: resetSubVolume calls resetClip
vmw = (root / "src/HFDIBDEM/virtualMesh/virtualMeshWall.C").read_text()
rs = re.search(
    r"void virtualMeshWall::resetSubVolume\(subVolumeProperties&[^)]*\)\s*\{(.*?)\}",
    vmw, re.S)
check("resetSubVolume calls resetClip",
      rs is not None and "resetClip()" in rs.group(0))

# 6: legacy call sites unchanged: volumeAndCentroid callers in src
vm_c = (root / "src/HFDIBDEM/virtualMesh/virtualMesh.C").read_text()
n_calls = len(re.findall(r"planePolyClip::volumeAndCentroid", vm_c))
check("virtualMesh.C keeps its 1 volumeAndCentroid call (prt path)",
      n_calls == 1)

print()
if failures:
    print(f"test_wc1_contracts: FAILED ({len(failures)} checks)")
    sys.exit(1)
print("test_wc1_contracts: all OK")
