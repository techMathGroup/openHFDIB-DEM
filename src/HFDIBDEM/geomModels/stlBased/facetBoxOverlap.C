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

Contributors
    Martin Isoz (2019-*)
\*---------------------------------------------------------------------------*/

#include "facetBoxOverlap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::facetBoxOverlapArea
(
    const point& a,
    const point& b,
    const point& c,
    const point& bMin,
    const point& bMax
)
{
    // clipping a triangle by 6 planes adds at most 6 vertices
    point bufA[9];
    point bufB[9];
    point* in(bufA);
    point* out(bufB);
    label nIn(3);
    in[0] = a;
    in[1] = b;
    in[2] = c;

    for (direction d = 0; d < 3 && nIn > 0; ++d)
    {
        for (int side = 0; side < 2 && nIn > 0; ++side)
        {
            const scalar lim(side == 0 ? bMin[d] : bMax[d]);
            const scalar sgn(side == 0 ? 1.0 : -1.0);
            label nOut(0);

            for (label i = 0; i < nIn; ++i)
            {
                const point& p0(in[i]);
                const point& p1(in[(i + 1) % nIn]);
                const scalar d0(sgn*(p0[d] - lim));
                const scalar d1(sgn*(p1[d] - lim));

                if (d0 >= 0)
                {
                    out[nOut++] = p0;
                }

                if ((d0 >= 0) != (d1 >= 0))
                {
                    out[nOut++] = p0 + (d0/(d0 - d1))*(p1 - p0);
                }
            }

            std::swap(in, out);
            nIn = nOut;
        }
    }

    vector aSum(Zero);
    for (label i = 1; i + 1 < nIn; ++i)
    {
        aSum += (in[i] - in[0]) ^ (in[i + 1] - in[0]);
    }

    return 0.5*mag(aSum);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
