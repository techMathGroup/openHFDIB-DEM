/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  | |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___/  |_____|_|  |_|
      | |                     H ybrid F ictitious D omain - I mmersed B oundary
      |_|                                        and D iscrete E lement M ethod
-------------------------------------------------------------------------------
License

    openHFDIB-DEM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (Version 3) as published
    by the Free Software Foundation.

    openHFDIB-DEM is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

Contributors
    Martin Isoz (2019-*)
\*---------------------------------------------------------------------------*/

#include "planePolyClip.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void planePolyClip::volumeAndCentroid
(
    const boundBox& box,
    const halfSpace& hs1,
    const bool hasSecond,
    const halfSpace& hs2,
    scalar& volume,
    vector& centroid
)
{
    clipAndIntegrate(box, hs1, hasSecond, hs2, volume, centroid, nullptr);
}


void planePolyClip::volumeCentroidAndFace
(
    const boundBox& box,
    const halfSpace& hs1,
    const bool hasSecond,
    const halfSpace& hs2,
    scalar& volume,
    vector& centroid,
    DynamicList<point>& facePolygon
)
{
    clipAndIntegrate
    (
        box,
        hs1,
        hasSecond,
        hs2,
        volume,
        centroid,
        &facePolygon
    );
}


void planePolyClip::clipAndIntegrate
(
    const boundBox& box,
    const halfSpace& hs1,
    const bool hasSecond,
    const halfSpace& hs2,
    scalar& volume,
    vector& centroid,
    DynamicList<point>* faceOut
)
{
    volume = 0;
    centroid = vector::zero;
    if (faceOut)
    {
        faceOut->clear();
    }

    // Note (MI): subvolumes are always hexes

    // vertex order must match boundBox::points() (hex model):
    // 0=(min,min,min) 1=(max,min,min) 2=(max,max,min) 3=(min,max,min)
    // 4=(min,min,max) 5=(max,min,max) 6=(max,max,max) 7=(min,max,max);
    // face corner lists below keep the counter-clockwise winding seen
    // from outside the box, so face normals point away from the centre
    // and a Newell normal of each clipped face tells whether the cut
    // polygon orientation must be flipped before the closing face is
    // appended
    // Note (MI): Newell normal stems from Newell's algorithm
    //            https://en.wikipedia.org/wiki/Newell%27s_algorithm
    static const label faceCorners[6][4] =
    {
        {0, 3, 2, 1},   // zmin, outward -z
        {4, 5, 6, 7},   // zmax, outward +z
        {0, 1, 5, 4},   // ymin, outward -y
        {2, 3, 7, 6},   // ymax, outward +y
        {1, 2, 6, 5},   // xmax, outward +x
        {0, 4, 7, 3}    // xmin, outward -x
    };

    // start from the six box faces as polygons; sign tolerance: a
    // vertex exactly on a plane counts as kept, cut points are added
    // with the interpolation parameter bounded to [0, 1]
    const scalar tol(SMALL);
    const pointField corners(box.points());

    // face soup: list of closed polygons (outward winding)
    DynamicList<DynamicList<point>> faces;
    for (const auto& fc : faceCorners)
    {
        DynamicList<point> poly;
        for (label i = 0; i < 4; ++i)
        {
            poly.append(corners[fc[i]]);
        }
        faces.append(poly);
    }

    // clip the soup by each requested half-space
    for (label pass = 0; pass < (hasSecond ? 2 : 1); ++pass)
    {
        const halfSpace& hs(pass == 0 ? hs1 : hs2);

        DynamicList<DynamicList<point>> newFaces;
        DynamicList<point> closing;

        for (const auto& poly : faces)
        {
            DynamicList<point> kept;

            const label n(poly.size());
            for (label i = 0; i < n; ++i)
            {
                const point& a(poly[i]);
                const point& b(poly[(i + 1) % n]);

                const scalar da((a - hs.p) & hs.n);
                const scalar db((b - hs.p) & hs.n);

                const bool aIn(da >= -tol);
                const bool bIn(db >= -tol);

                if (aIn)
                {
                    kept.append(a);
                }

                if (aIn != bIn)
                {
                    // edge crosses the plane: add the crossing point to
                    // both the clipped face and (once) the closing ring
                    const scalar t(da/(da - db));
                    const point x(a + t*(b - a));
                    kept.append(x);
                    closing.append(x);
                }
            }

            if (kept.size() >= 3)
            {
                newFaces.append(kept);
            }
        }

        if (newFaces.size() == 0)
        {
            // everything cut away
            return;
        }

        if (closing.size() >= 3)
        {
            // order the closing ring by angle around its centroid, then
            // orient it against the plane normal so the face winds
            // outward (away from the kept material)
            vector c(vector::zero);
            forAll(closing, i)
            {
                c += closing[i];
            }
            c /= scalar(closing.size());

            vector e1(mag(hs.n.x()) > 0.5 ? vector(0, 1, 0) : vector(1, 0, 0));
            e1 = e1 - (e1 & hs.n)*hs.n;
            e1 /= mag(e1) + VSMALL;
            const vector e2(hs.n ^ e1);

            SortableList<scalar> angle(closing.size());
            forAll(closing, i)
            {
                const vector r(closing[i] - c);
                angle[i] = Foam::atan2(r & e2, r & e1);
            }
            angle.sort();

            DynamicList<point> closingFace;
            for (label k = closing.size() - 1; k >= 0; --k)
            {
                closingFace.append(closing[angle.indices()[k]]);
            }

            // the reversed-sorted ring winds with its geometric normal
            // along -n, i.e. outward from the kept material; append it
            // to close the surface. duplicate crossings (an edge shared
            // by two faces appears twice in the ring) are harmless:
            // they create zero-area triangles in the fan
            newFaces.append(closingFace);

            if (pass == 0 && faceOut)
            {
                // the requested face polygon: the kept region meets
                // the hs1 plane along this ring; later half-spaces
                // may clip it further (handled below the loop)
                faceOut->transfer(closingFace);
            }
        }

        faces = newFaces;

        if (faceOut && faceOut->size() > 0 && pass == 0 && hasSecond)
        {
            // hs2 may cut the hs1 face polygon: re-clip the captured
            // ring by the remaining half-space (kept side, same rule
            // as the face soup)
            DynamicList<point> keptFace;
            const label n(faceOut->size());
            for (label i = 0; i < n; ++i)
            {
                const point& a((*faceOut)[i]);
                const point& b((*faceOut)[(i + 1) % n]);

                const scalar da((a - hs2.p) & hs2.n);
                const scalar db((b - hs2.p) & hs2.n);

                if (da >= -tol)
                {
                    keptFace.append(a);
                }
                if ((da >= -tol) != (db >= -tol))
                {
                    const scalar t(da/(da - db));
                    keptFace.append(a + t*(b - a));
                }
            }
            faceOut->transfer(keptFace);
        }
    }

    if (faces.size() == 0)
    {
        return;
    }

    // integrate: signed tetrahedra from an anchor inside the region
    // (the box midpoint - always inside or on the boundary of the kept
    // convex region) over each face triangulation. the signed volume
    // det[a-p, b-p, c-p]/6 is positive for outward face winding, so no
    // tet is ever negative and the fan tiles the region exactly
    scalar v(0);
    vector m(vector::zero);

    const point anchor(box.midpoint());

    for (const auto& poly : faces)
    {
        const label n(poly.size());
        for (label i = 1; i < n - 1; ++i)
        {
            const point& a(poly[0]);
            const point& b(poly[i]);
            const point& c(poly[i + 1]);

            // signed volume x 6 of tet (anchor, a, b, c); positive when
            // (a, b, c) winds outward seen from the anchor
            const scalar vt
            (
                (1.0/6.0)*((a - anchor) & ((b - anchor) ^ (c - anchor)))
            );
            v += vt;
            m += vt*(0.25*(anchor + a + b + c));
        }
    }

    if (v > VSMALL)
    {
        volume = v;
        centroid = m/v;
    }
    else
    {
        volume = 0;
        centroid = vector::zero;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
