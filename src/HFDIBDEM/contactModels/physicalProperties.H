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
    class for storage of physical properties

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef physicalProperties_H
#define physicalProperties_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Struct physicalProperties Declaration
\*---------------------------------------------------------------------------*/

struct physicalProperties
{
    physicalProperties()
    :
        aY_(0),
        aG_(0),
        aMu_(0),
        maxAdhN_(0),
        curAdhN_(0),
        reduceM_(0),
        reduceBeta_(0)
    {}

    physicalProperties
    (
        const scalar aY,
        const scalar aG,
        const scalar aMu,
        const scalar maxAdhN,
        const scalar curAdhN,
        const scalar reduceM,
        const scalar reduceBeta
    ):
        aY_(aY),
        aG_(aG),
        aMu_(aMu),
        maxAdhN_(maxAdhN),
        curAdhN_(curAdhN),
        reduceM_(reduceM),
        reduceBeta_(reduceBeta)

    {}

    // Mean model parameters
    scalar aY_ = 0;         // Contact Young modulus
    scalar aG_ = 0;         // Contact shear modulus
    scalar aMu_ = 0;        // Contact frictional coef
    scalar maxAdhN_ = 0;    // Contact adhesive coef
    scalar curAdhN_ = 0;    // Current contact adhesive coef
    scalar reduceM_ = 0;    // Reduced Mass of the contact pair
    scalar reduceBeta_ = 0;  // Reduced damping coef of the contact pair
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
