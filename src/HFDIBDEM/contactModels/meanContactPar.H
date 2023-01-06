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
    class for storage of particle-wall contact variables

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef meanContactPar_H
#define meanContactPar_H

#include "dictionary.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class meanContactPar Declaration
\*---------------------------------------------------------------------------*/

class meanContactPar
{
    public:

    // Mean model parameters
    scalar aY_;         // Contact Young modulus
    scalar aG_;         // Contact shear modulus
    scalar aGammaN_;    // Contact normal gamma
    scalar aGammat_;    // Contact tangential gamma
    scalar aMu_;        // Contact frictional coef
    scalar maxAdhN_;    // Contact adhesive coef

    // Constructors ------------------------------------------------
    meanContactPar()
    {
        aY_ = 0;
        aG_ = 0;
        aGammaN_ = 0;
        aGammat_ = 0;
        aMu_ = 0;
        maxAdhN_ = 0;
    }

    meanContactPar
    (
        scalar aY,
        scalar aG,
        scalar aGammaN,
        scalar aGammat,
        scalar aMu,
        scalar maxAdhN
    )
    {
        aY_ = aY;
        aG_ = aG;
        aGammaN_ = aGammaN;
        aGammat_ = aGammat;
        aMu_ = aMu;
        maxAdhN_ = maxAdhN;
    }

    // Destructors -------------------------------------------------
    ~meanContactPar(){}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
