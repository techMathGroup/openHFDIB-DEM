#!/usr/bin/env python3
"""python replica of the openHFDIB-DEM particle-particle virtualMesh.

mirrors the branch 3.0_com implementation:
  src/HFDIBDEM/virtualMesh/virtualMesh.{H,C}    octree + contact eval
  src/HFDIBDEM/virtualMesh/subVolume.C          octree node
  src/HFDIBDEM/virtualMesh/planePolyClip.C      exact leaf clipping
  src/HFDIBDEM/geomModels/stlBased/stlBased.C   getVolumeType,
      limitFinalSubVolume, getShapesSurfacePlane,
      getClosestPointAndNormal
  src/HFDIBDEM/contactModels/prtContactInfo.C   pair bbox + subVolumeV
  src/HFDIBDEM/contactModels/prtContact.C       the <=4 edge gate

leaf rules:
  legacy (first order) : limitFinalSubVolume bbox clamp
  exact  (second order): exactSubVolume tangent-plane + planePolyClip
                         partial volumes (exactSubVolume switch)

writes legacy ascii vtk for: internal sub-volumes, mixed
sub-volumes, contact plane (swept hull), contact center, contact
normal.

deps: numpy, trimesh (rtree for contains). run:
  python3 vmReplica.py
"""
import math
import sys
from pathlib import Path

import numpy as np
import trimesh

# OpenFOAM double precision constants (doubleScalar.H):
# SMALL = eps, VSMALL = smallest normalized double
SMALL = float(np.finfo(np.float64).eps)
VSMALL = float(np.finfo(np.float64).tiny)
ROOTVSMALL = math.sqrt(VSMALL)

# volumeType (OpenFOAM volumeType.H order)
UNKNOWN = 0
INSIDE = 1
OUTSIDE = 2
MIXED = 3


def cross(a, b):
    return np.cross(a, b)


def dot(a, b):
    return float(np.dot(a, b))


def mag(v):
    return float(np.linalg.norm(v))


def normed(v):
    n = mag(v)
    if n < VSMALL:
        return np.zeros(3)
    return v / n


# ------------------------------------------------------------------ #
# triangleFuncs (OpenFOAM triangleFuncs.C)                            #
# ------------------------------------------------------------------ #
def intersect_axes_bundle(V0, V10, V20, i0, origins, max_length):
    """intersect axis-i0 parallel lines (through origins) with the
    triangle (V0, V0+V10, V0+V20); returns the hit point or None"""
    i1 = (i0 + 1) % 3
    i2 = (i1 + 1) % 3

    u1 = V10[i1]
    v1 = V10[i2]
    u2 = V20[i1]
    v2 = V20[i2]

    local_scale = abs(u1) + abs(v1) + abs(u2) + abs(v2)
    det = v2 * u1 - u2 * v1

    if local_scale < VSMALL or abs(det) / local_scale < SMALL:
        # triangle parallel to the axis direction
        return None

    with np.errstate(divide="ignore", invalid="ignore"):
        for P in origins:
            u0 = P[i1] - V0[i1]
            v0 = P[i2] - V0[i2]

            alpha = 0.0
            beta = 0.0
            inter = False

            if abs(u1) < ROOTVSMALL:
                beta = u0 / u2
                if 0.0 <= beta <= 1.0:
                    alpha = (v0 - beta * v2) / v1
                    inter = alpha >= 0.0 and alpha + beta <= 1.0
            else:
                beta = (v0 * u1 - u0 * v1) / det
                if 0.0 <= beta <= 1.0:
                    alpha = (u0 - beta * u2) / u1
                    inter = alpha >= 0.0 and alpha + beta <= 1.0

            if inter:
                p = V0 + alpha * V10 + beta * V20
                s = p[i0] - P[i0]
                if 0.0 <= s <= max_length:
                    return p

    return None


def seg_intersects_box(p0, p1, bmin, bmax):
    """treeBoundBox::intersects(p0, p1) - closed box segment test"""
    d = p1 - p0
    t0 = 0.0
    t1 = 1.0
    for i in range(3):
        if abs(d[i]) < VSMALL:
            if p0[i] < bmin[i] or p0[i] > bmax[i]:
                return False
        else:
            ta = (bmin[i] - p0[i]) / d[i]
            tb = (bmax[i] - p0[i]) / d[i]
            if ta > tb:
                ta, tb = tb, ta
            t0 = max(t0, ta)
            t1 = min(t1, tb)
            if t0 > t1:
                return False
    return True


def tri_intersects_bb(tri, bmin, bmax):
    """triangleFuncs::intersectBb(p0, p1, p2, cubeBb): the 12 cube
    edges versus the triangle, then the 3 triangle edges versus the
    box"""
    p0, p1, p2 = tri
    p10 = p1 - p0
    p20 = p2 - p0

    mn = np.asarray(bmin, dtype=float)
    mx = np.asarray(bmax, dtype=float)

    cube0 = mn
    cube1 = np.array([mn[0], mn[1], mx[2]])
    cube2 = np.array([mx[0], mn[1], mx[2]])
    cube3 = np.array([mx[0], mn[1], mn[2]])
    cube4 = np.array([mn[0], mx[1], mn[2]])
    cube5 = np.array([mn[0], mx[1], mx[2]])
    cube7 = np.array([mx[0], mx[1], mn[2]])

    if intersect_axes_bundle(p0, p10, p20, 0,
                             [cube0, cube1, cube5, cube4],
                             mx[0] - mn[0]) is not None:
        return True

    if intersect_axes_bundle(p0, p10, p20, 1,
                             [cube0, cube1, cube2, cube3],
                             mx[1] - mn[1]) is not None:
        return True

    if intersect_axes_bundle(p0, p10, p20, 2,
                             [cube0, cube3, cube7, cube4],
                             mx[2] - mn[2]) is not None:
        return True

    for a, b in ((p0, p1), (p1, p2), (p2, p0)):
        if seg_intersects_box(a, b, mn, mx):
            return True

    return False


# ------------------------------------------------------------------ #
# stlBased replica                                                    #
# ------------------------------------------------------------------ #
class StlBody:
    """stlBased geomModel replica: facet search, inside test, leaf
    plane, legacy bbox clamp"""

    def __init__(self, stl_path, name, transform=None):
        self.name = name
        mesh = trimesh.load(str(stl_path), process=True)
        if transform:
            if "scale" in transform:
                mesh.apply_scale(transform["scale"])
            if "rotate" in transform:
                axis, angle_deg = transform["rotate"]
                mesh.apply_transform(
                    trimesh.transformations.rotation_matrix(
                        math.radians(angle_deg), axis
                    )
                )
            if "translate" in transform:
                mesh.apply_translation(transform["translate"])

        self.mesh = mesh
        self.tri = np.asarray(mesh.triangles, dtype=float)
        self.n_facets = len(self.tri)
        # facet area normals (outward per the stl winding)
        self.area_n = np.asarray(mesh.face_normals, dtype=float) * \
            np.asarray(mesh.area_faces, dtype=float)[:, None]
        self.centres = self.tri.mean(axis=1)
        self.tri_bb_min = self.tri.min(axis=1)
        self.tri_bb_max = self.tri.max(axis=1)

        if mesh.is_watertight:
            self.com = np.asarray(mesh.center_mass, dtype=float)
        else:
            self.com = np.asarray(mesh.vertices, dtype=float).mean(axis=0)
        # dC_: max vertex distance from the CoM (the C++ takes the max
        # over surf-cell centres; used only as a search span here)
        self.dC = float(
            np.linalg.norm(
                np.asarray(mesh.vertices, dtype=float) - self.com,
                axis=1,
            ).max()
        )

    def bounds(self):
        v = np.asarray(self.mesh.vertices, dtype=float)
        return v.min(axis=0), v.max(axis=0)

    def point_inside(self, p):
        return bool(
            self.mesh.contains(np.array([p], dtype=float))[0]
        )

    def facet_overlaps(self, fi, bmin, bmax):
        """treeDataTriSurface::overlaps(index, cubeBb)"""
        tmin = self.tri_bb_min[fi]
        tmax = self.tri_bb_max[fi]
        if not (
            tmax[0] >= bmin[0] and tmin[0] <= bmax[0]
            and tmax[1] >= bmin[1] and tmin[1] <= bmax[1]
            and tmax[2] >= bmin[2] and tmin[2] <= bmax[2]
        ):
            return False
        tri = self.tri[fi]
        if ((tri >= bmin) & (tri <= bmax)).all(axis=1).any():
            return True
        return tri_intersects_bb(tri, bmin, bmax)

    def find_box(self, bmin, bmax):
        """indexedOctree findBox: all facets overlapping the box"""
        ov = (self.tri_bb_max >= bmin) & (self.tri_bb_min <= bmax)
        cand = np.where(ov.all(axis=1))[0]
        if len(cand) == 0:
            return []
        tri = self.tri[cand]
        inb = ((tri >= bmin) & (tri <= bmax)).all(axis=2)
        hit = list(cand[inb.any(axis=1)])
        rest = cand[~inb.any(axis=1)]
        for fi in rest:
            if tri_intersects_bb(self.tri[fi], bmin, bmax):
                hit.append(int(fi))
        return [int(f) for f in hit]

    def get_volume_type(self, sv, c_ib):
        """stlBased::getVolumeType with shapesIn caching through the
        parent chain"""
        info = sv.c_info if c_ib else sv.t_info
        if info.shapes_in is None:
            parent = sv.parent
            if parent is not None:
                p_info = parent.c_info if c_ib else parent.t_info
                if p_info.shapes_in is not None:
                    if len(p_info.shapes_in) == 0:
                        info.shapes_in = []
                    else:
                        info.shapes_in = [
                            f
                            for f in p_info.shapes_in
                            if self.facet_overlaps(
                                f, sv.bmin, sv.bmax
                            )
                        ]
                else:
                    info.shapes_in = self.find_box(sv.bmin, sv.bmax)
            else:
                info.shapes_in = self.find_box(sv.bmin, sv.bmax)

        if len(info.shapes_in) > 0:
            return MIXED
        return INSIDE if self.point_inside(sv.midpoint()) else OUTSIDE

    def get_intersection_points(self, fi, bmin, bmax):
        """stlBased::getIntersectionPoints: facet vertices inside the
        box, the three axes-bundle hits and the edge-box hits"""
        pts = []
        tri = self.tri[fi]
        for v in tri:
            if (v >= bmin).all() and (v <= bmax).all():
                pts.append(v.copy())

        p0, p1, p2 = tri
        p10 = p1 - p0
        p20 = p2 - p0
        mn = np.asarray(bmin, dtype=float)
        mx = np.asarray(bmax, dtype=float)
        cube0 = mn
        cube1 = np.array([mn[0], mn[1], mx[2]])
        cube2 = np.array([mx[0], mn[1], mx[2]])
        cube3 = np.array([mx[0], mn[1], mn[2]])
        cube4 = np.array([mn[0], mx[1], mn[2]])
        cube5 = np.array([mn[0], mx[1], mx[2]])
        cube7 = np.array([mx[0], mx[1], mn[2]])

        p = intersect_axes_bundle(p0, p10, p20, 0,
                                  [cube0, cube1, cube5, cube4],
                                  mx[0] - mn[0])
        if p is not None:
            pts.append(p)
        p = intersect_axes_bundle(p0, p10, p20, 1,
                                  [cube0, cube1, cube2, cube3],
                                  mx[1] - mn[1])
        if p is not None:
            pts.append(p)
        p = intersect_axes_bundle(p0, p10, p20, 2,
                                  [cube0, cube3, cube7, cube4],
                                  mx[2] - mn[2])
        if p is not None:
            pts.append(p)

        for a, b in ((p0, p1), (p1, p2), (p2, p0)):
            if seg_intersects_box(a, b, mn, mx):
                # the C++ appends the clipped entry point; the mean of
                # the collected points is what matters downstream
                if seg_box_entry(a, b, mn, mx) is not None:
                    pts.append(seg_box_entry(a, b, mn, mx))
        return pts

    def limit_final_sub_volume(self, sv, c_ib):
        """stlBased::limitFinalSubVolume: legacy (first order) leaf
        rule; returns the clamped box or None"""
        info = sv.c_info if c_ib else sv.t_info
        bmin = sv.bmin.copy()
        bmax = sv.bmax.copy()
        shapes_in = info.shapes_in
        if shapes_in is None or len(shapes_in) == 0:
            return None

        normal = np.zeros(3)
        pts = []
        for fi in shapes_in:
            pts.extend(self.get_intersection_points(fi, sv.bmin, sv.bmax))
            normal = normal + self.area_n[fi]

        if len(pts) == 0:
            return None

        closest = np.zeros(3)
        for p in pts:
            closest = closest + p
        closest /= len(pts)

        # round the accumulated normal to the closest axis and clamp
        mag_x = abs(normal[0])
        mag_y = abs(normal[1])
        mag_z = abs(normal[2])

        if mag_x > mag_y and mag_x > mag_z:
            if normal[0] < 0:
                bmin[0] = closest[0]
            else:
                bmax[0] = closest[0]
        elif mag_y > mag_x and mag_y > mag_z:
            if normal[1] < 0:
                bmin[1] = closest[1]
            else:
                bmax[1] = closest[1]
        elif mag_z > mag_x and mag_z > mag_y:
            if normal[2] < 0:
                bmin[2] = closest[2]
            else:
                bmax[2] = closest[2]
        else:
            return None

        return bmin, bmax

    def get_shapes_surface_plane(self, shapes_in, bmin, bmax):
        """stlBased::getShapesSurfacePlane: dominant-facet tangent
        plane; returns (p, n) with n the unit inward normal or None"""
        best_facet = -1
        best_area = 0.0
        leaf_min = np.asarray(bmin, dtype=float)
        leaf_max = np.asarray(bmax, dtype=float)
        for fi in shapes_in:
            ov_min = np.maximum(leaf_min, self.tri_bb_min[fi])
            ov_max = np.minimum(leaf_max, self.tri_bb_max[fi])
            span = ov_max - ov_min
            if (span > 0).all():
                vol = float(np.prod(span))
                if vol > best_area:
                    best_area = vol
                    best_facet = fi

        if best_facet == -1:
            return None

        area_n = self.area_n[best_facet]
        a_mag = mag(area_n)
        if a_mag < VSMALL:
            return None

        n_out = area_n / a_mag
        p = self.centres[best_facet].copy()

        # avgDim = mean of the box spans
        avg_dim = float(np.mean(leaf_max - leaf_min))
        offset = n_out * (avg_dim / 8.0 + VSMALL)

        pos_inside = self.point_inside(p + offset)
        neg_inside = self.point_inside(p - offset)
        if pos_inside == neg_inside:
            return None

        # n points into the body: planePolyClip keeps the
        # dot(x - p, n) >= 0 side, which is the body interior
        n = n_out if pos_inside else -n_out
        return p, n

    def get_leaf_sub_volume_plane(self, sv, c_ib):
        info = sv.c_info if c_ib else sv.t_info
        if info.shapes_in is None or len(info.shapes_in) == 0:
            return None
        return self.get_shapes_surface_plane(
            info.shapes_in, sv.bmin, sv.bmax
        )

    def get_closest_point_and_normal(self, start_point):
        """stlBased::getClosestPointAndNormal (span from dC)"""
        point, _dist, tri_idx = trimesh.proximity.closest_point(
            self.mesh, np.array([start_point], dtype=float)
        )
        return point[0], normed(self.area_n[tri_idx[0]])


def seg_box_entry(p0, p1, bmin, bmax):
    """entry point of segment (p0, p1) into the closed box"""
    d = p1 - p0
    t = 1.0
    hit = None
    for i in range(3):
        if abs(d[i]) < VSMALL:
            continue
        for face in (bmin[i], bmax[i]):
            ti = (face - p0[i]) / d[i]
            if 0.0 <= ti <= 1.0:
                q = p0 + ti * d
                ok = True
                for j in range(3):
                    if j == i:
                        continue
                    if q[j] < bmin[j] or q[j] > bmax[j]:
                        ok = False
                        break
                if ok and ti <= t:
                    t = ti
                    hit = q
    return hit


# ------------------------------------------------------------------ #
# planePolyClip replica                                               #
# ------------------------------------------------------------------ #
# boundBox::points() corner order (boundBox.C)
# 0=(min,min,min) 1=(max,min,min) 2=(max,max,min) 3=(min,max,min)
# 4=(min,min,max) 5=(max,min,max) 6=(max,max,max) 7=(min,max,max)
FACE_CORNERS = (
    (0, 3, 2, 1),   # zmin, outward -z
    (4, 5, 6, 7),   # zmax, outward +z
    (0, 1, 5, 4),   # ymin, outward -y
    (2, 3, 7, 6),   # ymax, outward +y
    (1, 2, 6, 5),   # xmax, outward +x
    (0, 4, 7, 3),   # xmin, outward -x
)


def clip_and_integrate(bmin, bmax, hs1, hs2, faces_out=None):
    """planePolyClip::clipAndIntegrate for the leaf box with one or
    two half-spaces; returns (volume, centroid) or (0, None). if
    faces_out is a list, the surviving clipped face polygons are
    appended to it (visualization)"""
    tol = SMALL
    # corners in boundBox::points order
    corners = [
        np.array([bmin[0], bmin[1], bmin[2]]),
        np.array([bmax[0], bmin[1], bmin[2]]),
        np.array([bmax[0], bmax[1], bmin[2]]),
        np.array([bmin[0], bmax[1], bmin[2]]),
        np.array([bmin[0], bmin[1], bmax[2]]),
        np.array([bmax[0], bmin[1], bmax[2]]),
        np.array([bmax[0], bmax[1], bmax[2]]),
        np.array([bmin[0], bmax[1], bmax[2]]),
    ]

    faces = [
        [corners[fc[i]].copy() for i in range(4)] for fc in FACE_CORNERS
    ]

    planes = (hs1,) if hs2 is None else (hs1, hs2)

    for hs in planes:
        p = hs[0]
        n = hs[1]
        new_faces = []
        closing = []
        for poly in faces:
            kept = []
            npoly = len(poly)
            for i in range(npoly):
                a = poly[i]
                b = poly[(i + 1) % npoly]
                da = dot(a - p, n)
                db = dot(b - p, n)
                a_in = da >= -tol
                b_in = db >= -tol
                if a_in:
                    kept.append(a)
                if a_in != b_in:
                    t = da / (da - db)
                    x = a + t * (b - a)
                    kept.append(x)
                    closing.append(x)
            if len(kept) >= 3:
                new_faces.append(kept)

        if len(new_faces) == 0:
            return 0.0, None

        if len(closing) >= 3:
            c = np.zeros(3)
            for q in closing:
                c = c + q
            c /= len(closing)

            e1 = np.array([0.0, 1.0, 0.0]) if abs(n[0]) > 0.5 \
                else np.array([1.0, 0.0, 0.0])
            e1 = e1 - dot(e1, n) * n
            e1 = e1 / (mag(e1) + VSMALL)
            e2 = cross(n, e1)

            closing_face = [
                closing[k]
                for k in sorted(
                    range(len(closing)),
                    key=lambda i: math.atan2(
                        dot(closing[i] - c, e2), dot(closing[i] - c, e1)
                    ),
                )
            ][::-1]
            new_faces.append(closing_face)

        faces = new_faces

    if len(faces) == 0:
        return 0.0, None

    if faces_out is not None:
        for poly in faces:
            faces_out.append([p.copy() for p in poly])

    # signed tetrahedra from the box midpoint (inside the kept convex
    # region)
    anchor = np.array(
        [0.5 * (bmin[i] + bmax[i]) for i in range(3)]
    )
    v = 0.0
    m = np.zeros(3)
    for poly in faces:
        n_tri = len(poly)
        for i in range(1, n_tri - 1):
            a = poly[0]
            b = poly[i]
            c = poly[i + 1]
            vt = (1.0 / 6.0) * dot(
                a - anchor, cross(b - anchor, c - anchor)
            )
            v += vt
            m = m + vt * (0.25 * (anchor + a + b + c))

    if v > VSMALL:
        return v, m / v
    return 0.0, None


# ------------------------------------------------------------------ #
# subVolume replica                                                   #
# ------------------------------------------------------------------ #
class SubVolumeInfo:
    """ibSubVolumeInfo"""

    def __init__(self):
        self.volume_type = UNKNOWN
        self.shapes_in = None
        self.clipped_volume = -1.0
        self.clipped_centroid = None

    def reset_clip(self):
        self.clipped_volume = -1.0
        self.clipped_centroid = None


class SubVolume:
    """octree node: lazy children (treeBoundBox::subBbox octant
    convention: bit0 right/x, bit1 top/y, bit2 front/z)"""

    __slots__ = ("bmin", "bmax", "c_info", "t_info", "children",
                 "parent", "is_edge", "clipped_faces")

    def __init__(self, bmin, bmax, parent=None):
        self.bmin = np.asarray(bmin, dtype=float)
        self.bmax = np.asarray(bmax, dtype=float)
        self.c_info = SubVolumeInfo()
        self.t_info = SubVolumeInfo()
        self.children = None
        self.parent = parent
        self.is_edge = False
        # exact-leaf visualization: the planePolyClip face soup of
        # this leaf (filled only on the exact path with capture)
        self.clipped_faces = None

    def ensure_clipped_faces(self):
        if self.clipped_faces is None:
            self.clipped_faces = []
        return self.clipped_faces

    def volume(self):
        s = self.bmax - self.bmin
        return float(s[0] * s[1] * s[2])

    def midpoint(self):
        return 0.5 * (self.bmin + self.bmax)

    def contains(self, p):
        return bool((p >= self.bmin).all() and (p <= self.bmax).all())

    def sub_bbox(self, octant, mid):
        bmin = self.bmin.copy()
        bmax = self.bmax.copy()
        if octant & 0x1:
            bmin[0] = mid[0]
        else:
            bmax[0] = mid[0]
        if octant & 0x2:
            bmin[1] = mid[1]
        else:
            bmax[1] = mid[1]
        if octant & 0x4:
            bmin[2] = mid[2]
        else:
            bmax[2] = mid[2]
        return bmin, bmax

    def child_sub_volumes(self):
        if self.children is None:
            mid = self.midpoint()
            self.children = [
                SubVolume(*self.sub_bbox(oct, mid), parent=self)
                for oct in range(8)
            ]
        return self.children

    def set_as_edge(self):
        self.is_edge = True


# ------------------------------------------------------------------ #
# virtualMesh replica                                                 #
# ------------------------------------------------------------------ #
class VirtualMesh:
    def __init__(self, c_body, t_body, root_bmin, root_bmax,
                 sub_volume_v, exact_sub_volume, max_sub_volumes=5e6,
                 starting_point=None):
        self.c_geom = c_body
        self.t_geom = t_body
        self.sub_volume_v = sub_volume_v
        self.exact_sub_volume = exact_sub_volume
        self.starting_point = (
            starting_point
            if starting_point is not None
            else 0.5 * (np.asarray(root_bmin) + np.asarray(root_bmax))
        )
        self.iter_count = 0
        root_volume = float(np.prod(
            np.asarray(root_bmax, dtype=float)
            - np.asarray(root_bmin, dtype=float)
        ))
        cap = 8.0 * root_volume / max(sub_volume_v, VSMALL)
        self.iter_max = int(min(cap, max(max_sub_volumes, 1.0)))
        self.truncated = False
        self.root = SubVolume(root_bmin, root_bmax)

        self.contact_volume = 0.0
        self.contact_center = np.zeros(3)
        self.edge_points = []

        # visualization collections
        self.internal_svs = []
        self.mixed_svs = []
        self.edge_svs = []

    def warn_truncated(self, where):
        self.truncated = True
        print(
            "WARNING: virtualMesh {}: octree visit cap {} reached"
            .format(where, self.iter_max),
            file=sys.stderr,
        )

    # -- classification (ensureVolumeType) --------------------------
    def ensure_volume_type(self, sv, c_ib):
        info = sv.c_info if c_ib else sv.t_info
        geom = self.c_geom if c_ib else self.t_geom
        if info.volume_type == UNKNOWN:
            info.volume_type = geom.get_volume_type(sv, c_ib)
        return info.volume_type

    # -- exact leaf evaluation (evaluateLeafExact) ------------------
    def evaluate_leaf_exact(self, sv, c_type, t_type):
        c_info = sv.c_info
        if c_info.clipped_volume > -0.5:
            return True

        c_plane = None
        t_plane = None
        if c_type == MIXED:
            c_plane = self.c_geom.get_leaf_sub_volume_plane(sv, True)
            if c_plane is None:
                return False
        if t_type == MIXED:
            t_plane = self.t_geom.get_leaf_sub_volume_plane(sv, False)
            if t_plane is None:
                return False

        if c_plane is not None:
            hs1 = c_plane
            hs2 = t_plane
        else:
            hs1 = t_plane
            hs2 = None

        vol, cen = clip_and_integrate(
            sv.bmin, sv.bmax, hs1, hs2,
            faces_out=sv.ensure_clipped_faces(),
        )
        c_info.clipped_volume = vol
        c_info.clipped_centroid = cen
        # note: the C++ (virtualMesh.C evaluateLeafExact) clips by
        # {c-plane, t-plane} when both bodies are MIXED, else the
        # single mixed plane - the same rule as here (both-plane
        # clip verified against the C++ on the sphere-lens pair,
        # vmPrtPrtConvergenceTest, 2026-10-06)
        sv.t_info.clipped_volume = vol
        sv.t_info.clipped_centroid = cen
        return True

    # -- detection (detectFirstVolumeInContact) ----------------------
    def detect_first_contact_point(self):
        found = False
        return self.detect_first_volume_in_contact(self.root, found)

    def detect_first_volume_in_contact(self, sv, start_point_found):
        self.iter_count += 1
        if self.iter_count > self.iter_max:
            if not self.truncated:
                self.warn_truncated("detectFirstVolumeInContact")
            return False

        if sv.volume() < self.sub_volume_v:
            c_type = self.ensure_volume_type(sv, True)
            t_type = self.ensure_volume_type(sv, False)

            if c_type == OUTSIDE or t_type == OUTSIDE:
                return False
            if c_type == INSIDE and t_type == INSIDE:
                return True

            if self.exact_sub_volume and self.evaluate_leaf_exact(
                sv, c_type, t_type
            ):
                return sv.c_info.clipped_volume > VSMALL

            c_bbox = (sv.bmin.copy(), sv.bmax.copy())
            t_bbox = (sv.bmin.copy(), sv.bmax.copy())

            if sv.c_info.volume_type == MIXED:
                r = self.c_geom.limit_final_sub_volume(sv, True)
                if r is None:
                    return False
                c_bbox = r
            if sv.t_info.volume_type == MIXED:
                r = self.t_geom.limit_final_sub_volume(sv, False)
                if r is None:
                    return False
                t_bbox = r

            return boxes_overlap(c_bbox, t_bbox)

        self.ensure_volume_type(sv, True)
        self.ensure_volume_type(sv, False)

        if (
            sv.c_info.volume_type == OUTSIDE
            or sv.t_info.volume_type == OUTSIDE
        ):
            return False
        if (
            sv.c_info.volume_type == INSIDE
            and sv.t_info.volume_type == INSIDE
        ):
            return True

        children = sv.child_sub_volumes()
        if not start_point_found:
            for i, ch in enumerate(children):
                if ch.contains(self.starting_point):
                    # rotate this child to the front
                    children.insert(0, children.pop(i))
                    break

        for ch in children:
            if self.detect_first_volume_in_contact(
                ch, start_point_found
            ):
                return True
            start_point_found = True

        return False

    # -- evaluation (evaluateContact / inspectSubVolume) -------------
    def evaluate_contact(self):
        self.contact_volume = 0.0
        self.edge_points = []
        self.contact_center = np.zeros(3)

        self.inspect_sub_volume(self.root)

        if self.contact_volume < VSMALL:
            return 0.0

        self.contact_center = self.contact_center / self.contact_volume
        return self.contact_volume

    def inspect_sub_volume(self, sv):
        self.iter_count += 1
        if self.iter_count > self.iter_max:
            if not self.truncated:
                self.warn_truncated("inspectSubVolume")
            return

        if sv.volume() < self.sub_volume_v:
            c_type = self.ensure_volume_type(sv, True)
            t_type = self.ensure_volume_type(sv, False)

            if c_type == OUTSIDE or t_type == OUTSIDE:
                return

            # visualization bookkeeping
            if c_type == MIXED or t_type == MIXED:
                self.mixed_svs.append(sv)
            else:
                self.internal_svs.append(sv)

            # exact leaf evaluation (second order)
            if (
                self.exact_sub_volume
                and (c_type == MIXED or t_type == MIXED)
                and self.evaluate_leaf_exact(sv, c_type, t_type)
            ):
                if c_type != INSIDE and t_type != INSIDE:
                    # rim leaf: MIXED for both bodies
                    self.edge_points.append(sv.midpoint())
                    sv.set_as_edge()
                    self.edge_svs.append(sv)

                clipped_v = sv.c_info.clipped_volume
                if clipped_v > VSMALL:
                    self.contact_volume += clipped_v
                    self.contact_center = (
                        self.contact_center
                        + sv.c_info.clipped_centroid * clipped_v
                    )
                return

            # legacy clamp (first order)
            c_bbox = (sv.bmin.copy(), sv.bmax.copy())
            t_bbox = (sv.bmin.copy(), sv.bmax.copy())
            c_volume_type_mixed = False
            if sv.c_info.volume_type == MIXED:
                c_volume_type_mixed = True
                r = self.c_geom.limit_final_sub_volume(sv, True)
                sv.c_info.volume_type = INSIDE if r is not None \
                    else OUTSIDE
                if r is not None:
                    c_bbox = r

            if sv.t_info.volume_type == MIXED:
                if c_volume_type_mixed:
                    self.edge_points.append(sv.midpoint())
                    sv.set_as_edge()
                    self.edge_svs.append(sv)
                r = self.t_geom.limit_final_sub_volume(sv, False)
                sv.t_info.volume_type = INSIDE if r is not None \
                    else OUTSIDE
                if r is not None:
                    t_bbox = r

            if (
                sv.c_info.volume_type == INSIDE
                and sv.t_info.volume_type == INSIDE
                and boxes_overlap(c_bbox, t_bbox)
            ):
                # boundBox cross section of cBBox and tBBox
                cb_min = np.maximum(c_bbox[0], t_bbox[0])
                cb_max = np.minimum(c_bbox[1], t_bbox[1])
                span = cb_max - cb_min
                if (span > 0).all():
                    cross_v = float(np.prod(span))
                    self.contact_volume += cross_v
                    self.contact_center = (
                        self.contact_center
                        + (0.5 * (cb_min + cb_max)) * cross_v
                    )
            return

        self.ensure_volume_type(sv, True)
        self.ensure_volume_type(sv, False)

        if (
            sv.c_info.volume_type == OUTSIDE
            or sv.t_info.volume_type == OUTSIDE
        ):
            return
        if (
            sv.c_info.volume_type == INSIDE
            and sv.t_info.volume_type == INSIDE
        ):
            self.contact_volume += sv.volume()
            self.contact_center = (
                self.contact_center + sv.midpoint() * sv.volume()
            )
            self.internal_svs.append(sv)
            return

        new_edge_points = []
        children = sv.child_sub_volumes()
        for ch in children:
            # child edge points stay local (newEdgePoints in the C++)
            n_before = len(self.edge_points)
            self.inspect_sub_volume(ch)
            new_edge_points.extend(self.edge_points[n_before:])

        if len(new_edge_points) > 0:
            return

        if (
            sv.c_info.volume_type == MIXED
            and sv.t_info.volume_type == MIXED
        ):
            self.edge_points.append(sv.midpoint())
            sv.set_as_edge()
            self.edge_svs.append(sv)

    # -- normal and area (get3DcontactNormalAndSurface) --------------
    def get_3d_contact_normal_and_surface(self):
        edge_pts = self.edge_points
        area = 0.0
        normal_vec = np.zeros(3)
        t_dc = self.t_geom.dC
        normal_vector = np.zeros(3)

        sub_contact_center = np.zeros(3)
        for p in edge_pts:
            sub_contact_center = sub_contact_center + p
        sub_contact_center /= len(edge_pts)

        closest_point, normal_vector = (
            self.t_geom.get_closest_point_and_normal(sub_contact_center)
        )

        if len(edge_pts) >= 3:
            norm_ok = False
            center = np.zeros(3)
            for p in edge_pts:
                center = center + p
            center /= len(edge_pts)

            xx = xy = xz = yy = yz = zz = 0.0
            for p in edge_pts:
                sp = p - center
                if mag(sp) > 0.0:
                    sp = sp / mag(sp)
                xx += sp[0] * sp[0]
                xy += sp[0] * sp[1]
                xz += sp[0] * sp[2]
                yy += sp[1] * sp[1]
                yz += sp[1] * sp[2]
                zz += sp[2] * sp[2]

            n_pts = len(edge_pts)
            xx /= n_pts
            xy /= n_pts
            xz /= n_pts
            yy /= n_pts
            yz /= n_pts
            zz /= n_pts

            # closed-form smallest-eigenvalue eigenvector, weighted
            # cross-product accumulation exactly as in the C++
            weighted_dir = np.zeros(3)

            det_x = yy * zz - yz * yz
            axis_dir_x = np.array(
                [det_x, xz * yz - xy * zz, xy * yz - xz * yy]
            )
            weight_x = det_x * det_x
            if dot(weighted_dir, axis_dir_x) < 0.0:
                weight_x = -weight_x
            weighted_dir = weighted_dir + axis_dir_x * weight_x

            det_y = xx * zz - xz * xz
            axis_dir_y = np.array(
                [xz * yz - xy * zz, det_y, xy * xz - yz * xx]
            )
            weight_y = det_y * det_y
            if dot(weighted_dir, axis_dir_y) < 0.0:
                weight_y = -weight_y
            weighted_dir = weighted_dir + axis_dir_y * weight_y

            det_z = xx * yy - xy * xy
            axis_dir_z = np.array(
                [xy * yz - xz * yy, xy * xz - yz * xx, det_z]
            )
            weight_z = det_z * det_z
            if dot(weighted_dir, axis_dir_z) < 0.0:
                weight_z = -weight_z
            weighted_dir = weighted_dir + axis_dir_z * weight_z

            if mag(weighted_dir) > SMALL:
                norm_ok = True
                normal_vec = weighted_dir / (mag(weighted_dir) + SMALL)
            if not norm_ok:
                normal_vec = normal_vector

            # best fitting plane through the contact center
            plane_n = normed(normal_vec)
            plane_p = self.contact_center

            def nearest_point(pt):
                # plane::nearestPoint
                return pt - plane_n * dot(pt - plane_p, plane_n)

            pts_in_plane = [nearest_point(p) for p in edge_pts]

            q1 = np.array([1.0, 0.0, 0.0])
            q2 = np.array([0.0, 1.0, 0.0])
            if abs(dot(q1, plane_n)) > abs(dot(q2, plane_n)):
                q1 = q2

            u = cross(plane_n, q1)
            v = cross(plane_n, u)

            # 4 x 15 angular sectors sweep around the contact center
            help_list = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0]
            comm_cells_in_sections = []

            for i in range(4):
                u_step = help_list[i + 1] - help_list[i]
                v_step = help_list[i + 2] - help_list[i + 1]
                for j in range(15):
                    # sector between the two boundary rays
                    udir1 = u * (
                        help_list[i] + u_step * j / 15.0
                    ) + v * (help_list[i + 1] + v_step * j / 15.0)
                    udir2 = u * (
                        help_list[i]
                        + u_step * (j + 1) / 15.0
                    ) + v * (
                        help_list[i + 1] + v_step * (j + 1) / 15.0
                    )
                    points_in_section = []
                    for q in pts_in_plane:
                        r = q - self.contact_center
                        # plane::sideOfPlane: NORMAL if >= 0, FLIP if
                        # < 0; keep uSide == NORMAL and vSide == FLIP
                        if dot(r, udir1) >= 0.0 and dot(r, udir2) < 0.0:
                            points_in_section.append(q)

                    if len(points_in_section) > 0:
                        max_pt = None
                        distance = 0.0
                        for pt in points_in_section:
                            m2 = mag(pt - self.contact_center)
                            if m2 > distance:
                                distance = m2
                                max_pt = pt
                        comm_cells_in_sections.append(max_pt)

            # contact area from the triangle fan
            for i in range(len(comm_cells_in_sections) - 1):
                ac = comm_cells_in_sections[i] - self.contact_center
                bc = (
                    comm_cells_in_sections[i + 1] - self.contact_center
                )
                area += mag(cross(ac, bc)) / 2.0

            if len(comm_cells_in_sections) > 2:
                ac = (
                    comm_cells_in_sections[-1] - self.contact_center
                )
                bc = comm_cells_in_sections[0] - self.contact_center
                area += mag(cross(ac, bc)) / 2.0

        if mag(normal_vec) == 0.0:
            normal_vec = normal_vector

        if dot(normal_vector, normal_vec) < 0:
            normal_vec = normal_vec * (-1.0)

        return area, normal_vec


def boxes_overlap(a, b):
    """boundBox::overlaps (closed)"""
    return bool(
        (a[1] >= b[0]).all() and (a[0] <= b[1]).all()
    )


# ------------------------------------------------------------------ #
# legacy ascii vtk writers                                             #
# ------------------------------------------------------------------ #
def write_vtk_boxes(path, boxes, label_name="label"):
    """hexahedron per box, in the writeCoverVtk style"""
    points = []
    polys = []
    data = []
    for bi, (bmin, bmax) in enumerate(boxes):
        base = len(points)
        # boundBox::points() order
        for corner in (
            (bmin[0], bmin[1], bmin[2]),
            (bmax[0], bmin[1], bmin[2]),
            (bmax[0], bmax[1], bmin[2]),
            (bmin[0], bmax[1], bmin[2]),
            (bmin[0], bmin[1], bmax[2]),
            (bmax[0], bmin[1], bmax[2]),
            (bmax[0], bmax[1], bmax[2]),
            (bmin[0], bmax[1], bmax[2]),
        ):
            points.append(corner)
        for fc in FACE_CORNERS:
            # one quad per box face; the cell scalar carries the box id
            polys.append(
                (bi, "4 {} {} {} {}".format(
                    base + fc[0], base + fc[1],
                    base + fc[2], base + fc[3],
                ))
            )

    with open(path, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("sub-volumes\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")
        f.write("POINTS {} double\n".format(len(points)))
        for p in points:
            f.write("{:.10e} {:.10e} {:.10e}\n".format(*p))
        f.write("POLYGONS {} {}\n".format(
            len(polys), len(polys) * 5
        ))
        for _, line in polys:
            f.write(line + "\n")
        if len(points) > 0:
            f.write("POINT_DATA {}\n".format(len(points)))
            f.write(
                "SCALARS {} int 1\nLOOKUP_TABLE default\n"
                .format(label_name)
            )
            for _ in points:
                f.write("0\n")
        if len(polys) > 0:
            f.write("CELL_DATA {}\n".format(len(polys)))
            f.write(
                "SCALARS {} int 1\nLOOKUP_TABLE default\n"
                .format(label_name)
            )
            for bi, _ in polys:
                f.write("{}\n".format(bi))


def write_vtk_points(path, pts, labels=None, label_name="label"):
    with open(path, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("points\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")
        f.write("POINTS {} double\n".format(len(pts)))
        for p in pts:
            f.write("{:.10e} {:.10e} {:.10e}\n".format(*p))
        if len(pts) > 0:
            f.write("POINT_DATA {}\n".format(len(pts)))
            f.write(
                "SCALARS {} int 1\nLOOKUP_TABLE default\n"
                .format(label_name)
            )
            if labels is None:
                for _ in pts:
                    f.write("0\n")
            else:
                for lb in labels:
                    f.write("{}\n".format(lb))


def write_vtk_normal(path, center, normal, length=None):
    """contact normal as an oriented line from the contact center"""
    if length is None:
        length = 3.0 * max(mag(normal), VSMALL)
    p0 = np.asarray(center, dtype=float)
    p1 = p0 + length * normed(np.asarray(normal, dtype=float))
    # geometry (LINES) must precede POINT_DATA in legacy vtk
    with open(path, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("contact normal\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")
        f.write("POINTS 2 double\n")
        f.write("{:.10e} {:.10e} {:.10e}\n".format(*p0))
        f.write("{:.10e} {:.10e} {:.10e}\n".format(*p1))
        f.write("LINES 1 3\n")
        f.write("2 0 1\n")


def write_vtk_polygon(path, poly):
    """closed polygon (list of 3d points) as a vtk polygon"""
    if poly is None or len(poly) < 3:
        return
    with open(path, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("contact plane\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")
        f.write("POINTS {} double\n".format(len(poly)))
        for p in poly:
            f.write("{:.10e} {:.10e} {:.10e}\n".format(*p))
        f.write("POLYGONS 1 {}\n".format(len(poly) + 1))
        f.write(
            "{} {}\n".format(
                len(poly), " ".join(str(i) for i in range(len(poly)))
            )
        )


def write_vtk_face_soup(path, entries, label_name="label"):
    """one polygon per clipped face; the cell scalar carries the leaf
    id so each clipped sub-volume polyhedron is selectable"""
    # entries: list of (sv, face-polygon) pairs
    pts = []
    polys = []  # (leaf_id, n, point ids)

    for sv, poly in entries:
        base = len(pts)
        for p in poly:
            pts.append(p)
        polys.append((id(sv), len(poly), list(range(base, len(pts)))))

    with open(path, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("clipped sub-volumes\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")
        f.write("POINTS {} double\n".format(len(pts)))
        for p in pts:
            f.write("{:.10e} {:.10e} {:.10e}\n".format(*p))
        f.write("POLYGONS {} {}\n".format(
            len(polys),
            sum(1 + n for _, n, _ in polys),
        ))
        for _, n, ids in polys:
            f.write(
                "{} {}\n".format(n, " ".join(str(i) for i in ids))
            )
        if len(polys) > 0:
            f.write("CELL_DATA {}\n".format(len(polys)))
            f.write(
                "SCALARS {} int 1\nLOOKUP_TABLE default\n"
                .format(label_name)
            )
            # unique leaf ids starting from 0
            leaf_ids = {}
            for sv_id, _, _ in polys:
                if sv_id not in leaf_ids:
                    leaf_ids[sv_id] = len(leaf_ids)
            for sv_id, _, _ in polys:
                f.write("{}\n".format(leaf_ids[sv_id]))


# ------------------------------------------------------------------ #
# driver                                                               #
# ------------------------------------------------------------------ #
# user-defined pair: two stl files, optional rigid transform per body
TESTS = Path(__file__).resolve().parents[2] / "virtualMesh"

BODIES = {
    "c": {
        # impellerSphere (body 0) at the recorded contact instant
        # (bodiesInfo/0.02/stlFiles). the log reference for this pair:
        #   center (0.01142363 0.04643099 -0.0031502753)
        #   normal (0.5039906 -0.85969249 0.081967934)
        #   volume 4.2637073e-13, area 4.7867579e-06
        "stl": TESTS / "impellerAndSphere_contact" / "case"
        / "bodiesInfo" / "0.02" / "stlFiles" / "0.stl",
        "transform": None,
    },
    "t": {
        # icoSphere (body 1) at the same instant
        "stl": TESTS / "impellerAndSphere_contact" / "case"
        / "bodiesInfo" / "0.02" / "stlFiles" / "1.stl",
        "transform": None,
    },
}

# virtualMesh dict settings (openHFDIBDEM.C:454)
VM_SETTINGS = {
    "level": 2,                 # virtualMeshLevel
    "charCellSize": 0.0012,      # charCellSize
    "maxSubVolumes": 5e6,        # maxSubVolumes
    "exactSubVolume": True,      # False = legacy first order
}


def build_pair_bbox(c_body, t_body):
    """prtContactInfo::getContacts_ArbShape: c-body bbox clipped to
    the t-body bbox, per-coordinate"""
    cmin, cmax = c_body.bounds()
    tmin, tmax = t_body.bounds()

    bmin = np.empty(3)
    bmax = np.empty(3)
    for i in range(3):
        bmin[i] = cmin[i]
        bmax[i] = cmax[i]
        if bmin[i] < tmin[i]:
            bmin[i] = tmin[i]
        if bmax[i] > tmax[i]:
            bmax[i] = tmax[i]
    return bmin, bmax


def run(settings, bodies_cfg, out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    c_body = StlBody(
        bodies_cfg["c"]["stl"], "c", bodies_cfg["c"].get("transform")
    )
    t_body = StlBody(
        bodies_cfg["t"]["stl"], "t", bodies_cfg["t"].get("transform")
    )

    # prtContact.C:77: ranCellVol = charCellSize^3
    level_of_division = int(2 ** settings["level"])
    char_cell = settings["charCellSize"]
    sub_volume_len = char_cell / level_of_division
    sub_volume_v = sub_volume_len ** 3

    bmin, bmax = build_pair_bbox(c_body, t_body)
    if (bmax <= bmin).any():
        print("no pair-bbox overlap - nothing to do")
        return
    print("pair bbox: [{} .. {}]".format(bmin, bmax))
    print(
        "sub-volume: len {:.6e} vol {:.6e} (level {}, ccs {})".format(
            sub_volume_len, sub_volume_v,
            settings["level"], char_cell,
        )
    )

    vm = VirtualMesh(
        c_body,
        t_body,
        bmin,
        bmax,
        sub_volume_v,
        settings["exactSubVolume"],
        settings["maxSubVolumes"],
    )

    detected = vm.detect_first_contact_point()
    print("contact detected: {}".format(detected))

    if not detected:
        print("no contact - nothing to write")
        return

    vol = vm.evaluate_contact()
    print("sub-volumes visited: {}".format(vm.iter_count))

    # prtContact.C:354: the <=4 edge-point gate
    edge_gate = len(vm.edge_points) <= 4
    if edge_gate:
        print(
            "edge-point gate: {} edge points (<= 4) -> volume zeroed"
            .format(len(vm.edge_points))
        )
        vol = 0.0

    print("contact volume: {:.6e}".format(vol))
    print("contact center: {}".format(vm.contact_center))

    if vol > 0:
        area, normal = vm.get_3d_contact_normal_and_surface()
        print("contact area:   {:.6e}".format(area))
        print("contact normal: {}".format(normal))
    else:
        area, normal = 0.0, np.zeros(3)

    # ---- vtk output -------------------------------------------------
    write_vtk_boxes(
        out / "internalSubVolumes.vtk",
        [(sv.bmin, sv.bmax) for sv in vm.internal_svs],
        "internal",
    )
    write_vtk_boxes(
        out / "mixedSubVolumes.vtk",
        [(sv.bmin, sv.bmax) for sv in vm.mixed_svs],
        "mixed",
    )
    # exact path: the planePolyClip face soups of the MIXED leaves -
    # the actual clipped sub-volume polyhedra (vs the raw octant
    # boxes in mixedSubVolumes.vtk)
    if settings["exactSubVolume"]:
        clipped = [
            (sv, poly)
            for sv in vm.mixed_svs
            if sv.clipped_faces
            for poly in sv.clipped_faces
        ]
        write_vtk_face_soup(
            out / "clippedSubVolumes.vtk", clipped, "clipped"
        )
    write_vtk_points(
        out / "contactCenter.vtk", [vm.contact_center], None, "center"
    )
    if vol > 0 and mag(normal) > 0:
        # arrow scaled to the pair bbox
        diag = mag(bmax - bmin)
        write_vtk_normal(
            out / "contactNormal.vtk",
            vm.contact_center,
            normal,
            length=0.5 * diag,
        )
    # contact plane: the projected edge points in the best-fit plane
    # are the swept hull; reconstruct them as the sector hull
    if vol > 0 and len(vm.edge_points) >= 3:
        plane_pts = project_to_best_fit_plane(
            vm.edge_points, vm.contact_center, normal
        )
        hull = sector_hull(plane_pts, vm.contact_center, normal)
        write_vtk_polygon(out / "contactPlane.vtk", hull)
    else:
        write_vtk_points(
            out / "edgePoints.vtk",
            vm.edge_points,
            None,
            "edge",
        )

    print("vtk output written to {}".format(out))


def project_to_best_fit_plane(pts, center, normal):
    n = normed(normal)
    return [p - n * dot(p - center, n) for p in pts]


def sector_hull(pts_in_plane, center, normal):
    """outermost point per angular sector (the swept hull of the
    C++ area sweep); reuses the 4x15 sector rule"""
    n = normed(normal)
    q1 = np.array([1.0, 0.0, 0.0])
    q2 = np.array([0.0, 1.0, 0.0])
    if abs(dot(q1, n)) > abs(dot(q2, n)):
        q1 = q2
    u = cross(n, q1)
    v = cross(n, u)

    help_list = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0]
    hull = []
    for i in range(4):
        u_step = help_list[i + 1] - help_list[i]
        v_step = help_list[i + 2] - help_list[i + 1]
        for j in range(15):
            udir1 = u * (help_list[i] + u_step * j / 15.0) + v * (
                help_list[i + 1] + v_step * j / 15.0
            )
            udir2 = u * (
                help_list[i] + u_step * (j + 1) / 15.0
            ) + v * (
                help_list[i + 1] + v_step * (j + 1) / 15.0
            )
            in_sec = []
            for q in pts_in_plane:
                r = q - center
                # same NORMAL/FLIP rule as the area sweep
                if dot(r, udir1) >= 0.0 and dot(r, udir2) < 0.0:
                    in_sec.append(q)
            if in_sec:
                best = max(
                    in_sec, key=lambda q: mag(q - center)
                )
                hull.append(best)
    return hull


def main():
    run(VM_SETTINGS, BODIES, TESTS / "vmReplicaOut")


if __name__ == "__main__":
    main()
