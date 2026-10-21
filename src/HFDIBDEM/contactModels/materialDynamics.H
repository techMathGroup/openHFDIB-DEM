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

    You should have received a copy of the GNU General Public License
    along with openHFDIB-DEM. If not, see <http://www.gnu.org/licenses/>.

InNamespace
    Foam

Description
    Dynamic material parameters of the surface-loading contact law: density,
    compressive flow stress, strength-cap factor and the Hugoniot of the
    material.

    Deliberately free of OpenFOAM dependencies so the law that consumes it
    (surfaceLoadingLaw.H) can be compiled and tested on its own; see
    tests/surfaceLoading/.

    Entry names in DEM/materials/<name>:
        rho          density [kg/m3]                    (required for a pair)
        Yc           compressive flow stress [Pa]       (required for a pair)
        capFactor    strength-cap factor, p = capFactor*Yc (default 2.7)
        hugoniot { c0 s c0High sHigh uMelt }            (required for a pair)

    The Hugoniot is D = c0 + s*u_p with an optional second (melt) segment
    D = c0High + sHigh*u_p for u_p >= uMelt. A missing second segment means
    the first one applies at all particle velocities.

Contributors
    Martin Isoz (2019-*)
\*---------------------------------------------------------------------------*/

#ifndef materialDynamics_H
#define materialDynamics_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Struct materialDynamics Declaration
\*---------------------------------------------------------------------------*/

struct materialDynamics
{
    // density [kg/m3]; < 0 = not given
    double rho = -1.0;

    // compressive flow stress [Pa]; < 0 = not given
    double Yc = -1.0;

    // strength-cap factor, p_str = capFactor*Yc
    double capFactor = 2.7;

    // Hugoniot D = c0 + s*u_p, optional melt segment above uMelt
    bool hasHugoniot = false;
    double c0 = 0.0;      // [m/s]
    double s = 0.0;       // [-]
    double c0High = 0.0;  // [m/s]
    double sHigh = 0.0;   // [-]
    double uMelt = 0.0;   // [m/s]

    bool twoSegment() const
    {
        return hasHugoniot && c0High > 0.0 && uMelt > 0.0;
    }

    bool usable() const
    {
        return rho > 0.0 && Yc > 0.0 && hasHugoniot && c0 > 0.0 && s > 0.0;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
