#!/usr/bin/env python3
import re
import unittest
from pathlib import Path


class TestSubVolumeContracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.repo_root = Path(__file__).resolve().parents[3]
        cls.header = (
            cls.repo_root
            / "src"
            / "HFDIBDEM"
            / "virtualMesh"
            / "subVolume.H"
        ).read_text()
        cls.source = (
            cls.repo_root
            / "src"
            / "HFDIBDEM"
            / "virtualMesh"
            / "subVolume.C"
        ).read_text()

    def test_default_ctor_initializes_isedge_false(self):
        pattern = re.compile(
            r"subVolume::subVolume\(\)\s*:\s*[\s\S]*?isEdge_\(false\)",
            re.MULTILINE,
        )
        self.assertRegex(
            self.source,
            pattern,
            "subVolume default constructor must initialize isEdge_ to false",
        )

    def test_ibsubvolumeinfo_ctor_has_safe_fallback(self):
        has_init = bool(
            re.search(
                r"ibSubVolumeInfo\(volumeType\s+vType\)\s*:\s*"
                r"[\s\S]*?volumeType_\(volumeType::UNKNOWN\)",
                self.header,
                re.MULTILINE,
            )
        )
        self.assertTrue(
            has_init,
            "ibSubVolumeInfo(volumeType) should initialize volumeType_ in init list",
        )

        ctor_start = self.header.find("ibSubVolumeInfo(volumeType vType)")
        self.assertNotEqual(ctor_start, -1)
        ctor_block = self.header[ctor_start: ctor_start + 1200]
        self.assertIn(
            "default:",
            ctor_block,
            "ibSubVolumeInfo(volumeType) switch should include default case",
        )

    def test_parent_pointer_is_weak_ptr(self):
        self.assertIn(
            "std::weak_ptr<subVolume> parentSV_;",
            self.header,
            "parentSV_ should use weak_ptr to avoid ownership cycles",
        )


if __name__ == "__main__":
    unittest.main()
