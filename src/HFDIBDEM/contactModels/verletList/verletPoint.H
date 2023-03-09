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
    class for Verlet point

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef verletPoint_H
#define verletPoint_H

#include "dictionary.H"
#include "IDLList.H"
#include "fvCFD.H"
#include "geomModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class verletPoint Declaration
\*---------------------------------------------------------------------------*/

class verletPoint
:
    public IDLList<verletPoint>::link
{
    private:

    label bodyID_;

    point& pointLoc_;
    const bool isMin_;

    public:

    // Constructors ------------------------------------------------
    verletPoint
    (
        label bodyID,
        point& pointLoc,
        const bool isMin
    )
    :
    bodyID_(bodyID),
    pointLoc_(pointLoc),
    isMin_(isMin)
    {};

    // Destructors -------------------------------------------------
    virtual ~verletPoint(){};

    bool isMin()
    {
        return isMin_;
    }

    point& getPoint()
    {
        return pointLoc_;
    }

    label getBodyId()
    {
        return bodyID_;
    }

    autoPtr<verletPoint> clone() const
    {
        return autoPtr<verletPoint>(new verletPoint(*this));
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //