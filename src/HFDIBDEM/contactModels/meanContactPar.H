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

    openHFDIB-DEM is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

Description
    class for storage of physical properties

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
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
        reduceBeta_(reduceBeta),
    {}

    // Mean model parame,ters
    scalar aY_ = 0;         // Contact Young modulus
    scalar aG_ = 0;         // Contact shear modulus
    scalar aMu_ = 0;        // Contact frictional coef
    scalar maxAdhN_ = 0;    // Contact adhesive coef
    scalar curAdhN_ = 0;    // Current contact adhesive coef
    scalar reduceM_ = 0;    // Reduced Mass of the contact pair
    scalar reduceBeta_ = 0;    // Reduced dissipation coeffitient
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
