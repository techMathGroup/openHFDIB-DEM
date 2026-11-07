#!/usr/bin/env python3
"""analytic validation of vmReplica: two overlapping spheres.

the sphere-sphere overlap (lens) has closed-form volume and contact
patch area; the replica (both leaf rules) is compared against them.
"""
import math
import sys
import tempfile
from pathlib import Path

import numpy as np
import trimesh

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vmReplica import StlBody, VirtualMesh  # noqa: E402


def lens_volume(R, r, d):
    """volume of the intersection of two spheres (radii R, r,
    center distance d)"""
    return (
        math.pi * (R + r - d) ** 2
        * (d * d + 2 * d * r - 3 * r * r + 2 * d * R + 6 * r * R
           - 3 * R * R)
        / (12 * d)
    )


def main():
    R = 0.005
    d = 0.008

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        spheres = []
        for i, off in enumerate((0.0, d)):
            m = trimesh.creation.icosphere(
                subdivisions=3, radius=R
            )
            m.apply_translation([off, 0.0, 0.0])
            p = td / f"sphere{i}.stl"
            m.export(p)
            spheres.append(StlBody(p, f"s{i}"))

        c_body, t_body = spheres

        char_cell = 0.002
        for level in (2, 3):
            sub_len = char_cell / int(2 ** level)
            sub_v = sub_len ** 3

            # pair bbox: c bbox clipped to t bbox
            cmin, cmax = c_body.bounds()
            tmin, tmax = t_body.bounds()
            bmin = np.maximum(cmin, tmin)
            bmax = np.minimum(cmax, tmax)

            for exact in (True, False):
                vm = VirtualMesh(
                    c_body, t_body, bmin, bmax, sub_v, exact
                )
                det = vm.detect_first_contact_point()
                vol = vm.evaluate_contact() if det else 0.0
                if len(vm.edge_points) <= 4:
                    vol = 0.0
                area, normal = (0.0, np.zeros(3))
                if vol > 0:
                    area, normal = vm.get_3d_contact_normal_and_surface()

                ana_v = lens_volume(R, R, d)
                # contact patch: circle of radius sqrt(R^2 - (d/2)^2)
                ana_a = math.pi * (R * R - (d / 2) ** 2)

                print(
                    f"level {level} ({'exact' if exact else 'legacy'}):"
                    f" V={vol:.4e} (analytic {ana_v:.4e},"
                    f" {vol / ana_v:6.1%})"
                    f" A={area:.4e} (analytic {ana_a:.4e},"
                    f" {area / ana_a:6.1%})"
                    f" n={normal}"
                    f" visited={vm.iter_count}"
                )


if __name__ == "__main__":
    main()
