#!/usr/bin/env python3
"""analytic validation of vmWallReplica: sphere penetrating a wall.

the spherical-cap contact has closed-form volume and wetted footprint
area; the replica (both leaf rules) is compared against them.
"""
import math
import sys
import tempfile
from pathlib import Path

import numpy as np
import trimesh

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vmReplica import StlBody  # noqa: E402
import vmWallReplica as vwr  # noqa: E402


def cap_volume(R, d):
    """spherical cap volume (penetration depth d, radius R)"""
    return math.pi * d * d * (R - d / 3.0)


def cap_area(R, d):
    """wetted footprint area: circle of radius sqrt(R^2-(R-d)^2)"""
    return math.pi * (2.0 * R * d - d * d)


def main():
    R = 0.005
    d = 0.0008

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        m = trimesh.creation.icosphere(subdivisions=4, radius=R)
        # bottom wall at y = 0, fluid above, normal into the wall
        m.apply_translation([0.0, R - d, 0.0])
        stl = td / "sphere.stl"
        m.export(stl)

        body = StlBody(stl, "body", None)
        walls = {
            "bottom": (
                np.array([0.0, 0.0, 0.0]),
                np.array([0.0, -1.0, 0.0]),
            ),
        }

        v_exact = cap_volume(R, d)
        a_exact = cap_area(R, d)
        print("analytic: volume {:.6e} area {:.6e}".format(
            v_exact, a_exact
        ))

        for level in (1, 2, 3):
            for exact in (False, True):
                settings = {
                    "level": level,
                    "charCellSize": 0.0012,
                    "maxSubVolumes": 5e6,
                    "exactSubVolume": exact,
                }
                out = td / "out_L{}_{}".format(level, exact)
                vwr.run(settings, {
                    "stl": stl, "transform": None,
                }, walls, out)

                # re-run the single sub-contact to read the results
                wall_ci = vwr.WallContactInfo(
                    body, walls, settings["charCellSize"],
                    exact_sub_volume=exact,
                )
                wall_ci.detect_wall_contact()
                sub_contacts = wall_ci.find_contact_areas()
                res, _ = vwr.get_wall_contact_vars_arb_shape(
                    body, settings, sub_contacts[0], walls
                )
                tag = "exact" if exact else "legacy"
                print(
                    "level {} {:>6}: V {:.4f}x  A {:.4f}x  "
                    "n = {}".format(
                        level, tag,
                        res["volume"] / v_exact,
                        res["area"] / a_exact,
                        np.round(res["normal"], 6),
                    )
                )


if __name__ == "__main__":
    main()
