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

InNamespace
    Foam

Class
    Foam::planePolyClip

Description
    exact partial-volume credit for one virtualMesh leaf (subVolume): 
    clips an axis-aligned box by up to two half-spaces and integrates
    the volume and centroid of the kept region.

    - usable for MIXED subvolumes
    - build using Sutherland-Hodgman polygon clipping

    Technical Note:
    each half-space is {p, n}: the kept side is (x - p) & n >= 0, so n
    must point INTO the kept material (the body interior for virtualMesh
    contact evaluation).

SourceFiles
    planePolyClip.C

Contributors
    Martin Isoz (2019-*)
\*---------------------------------------------------------------------------*/

#ifndef planePolyClip_H
#define planePolyClip_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Struct planePolyClip Declaration
\*---------------------------------------------------------------------------*/

struct planePolyClip
{
    // one half-space: kept side is (x - p) & n >= 0, n unit, points into
    // the kept material
    struct halfSpace
    {
        point p;
        vector n;

        halfSpace(const point& p, const vector& n)
        :
        p(p),
        n(n)
        {}
    };

    // clip box by the half-spaces and integrate the kept region.
    // pass hasSecond = false for the single-plane variant. volume is 0
    // when the region is empty or degenerate; the result is exact for
    // any convex region (box cuts), no axis rounding involved
    static void volumeAndCentroid
    (
        const boundBox& box,
        const halfSpace& hs1,
        const bool hasSecond,
        const halfSpace& hs2,
        scalar& volume,
        vector& centroid
    );

    // as volumeAndCentroid, plus the polygon in which the kept
    // region meets the plane of hs1 (the closing face of the first
    // clip). the polygon is empty when hs1 does not cut the box or
    // the region does not survive the clip. callers that need the
    // face in a specific plane (e.g. the wetted wall face of a
    // contact leaf) must pass that half-space as hs1
    static void volumeCentroidAndFace
    (
        const boundBox& box,
        const halfSpace& hs1,
        const bool hasSecond,
        const halfSpace& hs2,
        scalar& volume,
        vector& centroid,
        DynamicList<point>& facePolygon
    );

    private:

    // shared clip + integration core; faceOut (null when no face is
    // requested) receives the closing polygon of the hs1 clip
    static void clipAndIntegrate
    (
        const boundBox& box,
        const halfSpace& hs1,
        const bool hasSecond,
        const halfSpace& hs2,
        scalar& volume,
        vector& centroid,
        DynamicList<point>* faceOut
    );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
