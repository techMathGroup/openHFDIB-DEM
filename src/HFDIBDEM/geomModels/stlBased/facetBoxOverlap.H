/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| "_ \ / _ \ "_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
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

Description
    exact overlap of one triangular facet with an axis-aligned box.

    facetBoxOverlapArea clips the triangle by the six box planes
    (Sutherland-Hodgman) and returns the area of the surviving
    polygon. unlike a bbox-volume weight it stays positive for
    facets lying flat on a lattice plane (zero-thickness bbox), so
    it is a safe ranking weight for the dominant-facet pick of
    stlBased::getShapesSurfacePlane. the area is 0 when the facet
    misses the box or only touches it in a point or line.

SourceFiles
    facetBoxOverlap.C

Contributors
    Martin Isoz (2019-*)
\*---------------------------------------------------------------------------*/

#ifndef facetBoxOverlap_H
#define facetBoxOverlap_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- area of the triangle (a, b, c) clipped to the box [bMin, bMax]
scalar facetBoxOverlapArea
(
    const point& a,
    const point& b,
    const point& c,
    const point& bMin,
    const point& bMax
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
