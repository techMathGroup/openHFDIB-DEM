#!/usr/bin/env python3
"""T2 contract tests for the exactSubVolume (approach B) wiring.

these are structural checks in the grep style of
test_subVolume_contracts.py: the opt-in keyword must be parsed, the
plane accessor must be declared and overridden, the leaf cache must be
deep-copied, and the virtualMesh leaf paths must guard on the switch.
"""
import re
import unittest
from pathlib import Path


class TestExactSubVolumeContracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = Path(__file__).resolve().parents[3]
        cls.src = cls.root / "src" / "HFDIBDEM"

    @classmethod
    def read(cls, *parts):
        return cls.src.joinpath(*parts).read_text()

    def test_keyword_parsed_in_openhfdibdem(self):
        text = self.read("openHFDIBDEM.C")
        self.assertIn(
            '"exactSubVolume"',
            text,
            "virtualMesh dict must parse the exactSubVolume keyword",
        )
        self.assertIn(
            "setExactSubVolume",
            text,
            "parsed keyword must reach virtualMeshLevel",
        )

    def test_switch_stored_in_virtualmeshlevel(self):
        text = self.read("parameters", "virtualMeshLevel.H")
        self.assertIn(
            "static bool exactSubVolume_",
            text,
            "virtualMeshLevel must hold the exactSubVolume switch",
        )
        self.assertIn(
            "getExactSubVolume",
            text,
            "virtualMeshLevel must expose getExactSubVolume()",
        )

    def test_plane_accessor_declared_in_geommodel(self):
        text = self.read("geomModels", "geomModel.H")
        self.assertRegex(
            text,
            r"virtual bool getLeafSubVolumePlane",
            "geomModel must declare the getLeafSubVolumePlane virtual",
        )

    def test_plane_accessor_overridden_in_stlbased(self):
        text = self.read("geomModels", "stlBased", "stlBased.H")
        self.assertRegex(
            text,
            r"bool getLeafSubVolumePlane[\s\S]*?final override",
            "stlBased must override getLeafSubVolumePlane",
        )

    def test_plane_clip_utility_registered(self):
        text = (self.root / "src" / "HFDIBDEM" / "Make" / "files").read_text()
        self.assertIn(
            "virtualMesh/planePolyClip.C",
            text,
            "planePolyClip.C must be in the library build",
        )

    def test_leaf_cache_deep_copied(self):
        text = self.read("virtualMesh", "subVolume.H")
        self.assertRegex(
            text,
            r"ibSubVolumeInfo\(const ibSubVolumeInfo& other\)[\s\S]*?"
            r"clippedVolume_\(other\.clippedVolume_\)",
            "ibSubVolumeInfo copy ctor must copy the clipped-volume cache",
        )

    def test_subcontact_uses_clipped_volume(self):
        text = self.read("virtualMesh", "subContact.C")
        self.assertIn(
            "clippedVolume_",
            text,
            "subContact::addSubVolume must accumulate the cached exact"
            " leaf volume when available (W6)",
        )

    def test_leaf_paths_guard_on_switch(self):
        text = self.read("virtualMesh", "virtualMesh.C")
        self.assertEqual(
            text.count("getExactSubVolume()"),
            2,
            "both inspectSubVolume and detectFirstVolumeInContact"
            " leaf paths must guard on the exactSubVolume switch",
        )
        self.assertIn(
            "evaluateLeafExact",
            text,
            "the shared exact leaf evaluation helper must exist",
        )


if __name__ == "__main__":
    unittest.main()
