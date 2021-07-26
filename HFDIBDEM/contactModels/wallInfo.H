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
    class for storage of contact variables for wall

SourceFiles
    wallInfo.C

Contributors
    Martin Isoz (2019-*), Martin Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef wallInfo_H
#define wallInfo_H

#include "dictionary.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class wallInfo Declaration
\*---------------------------------------------------------------------------*/

class wallInfo
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

    protected:

    //spring stiffness for normal contact force
    scalar  kN_;

    //spring stiffness for tangential contact force
    scalar  kt_;

    //damping coefficient for normal contact force
    scalar gammaN_;

    //damping coefficient for tangential contact force
    scalar gammat_;

    //friction coefficient
    scalar mu_;

    //adhesion coefficient for normal contact force
    scalar adhN_;

    //adhesion coefficient for normal contact force
    scalar adhEqui_;

    public:

    // Constructors ------------------------------------------------
    wallInfo
    (
        scalar kN,
        scalar kt,
        scalar gammaN,
        scalar gammat,
        scalar mu,
        scalar adhN,
        scalar adhEqui
    );

    // Destructors -------------------------------------------------
    ~wallInfo();

    const scalar& getkN()   const
    {
        return kN_;
    }

    const scalar& getkt()   const
    {
        return kt_;
    }

    const scalar& getgammaN()   const
    {
        return gammaN_;
    }

    const scalar& getgammat()   const
    {
        return gammat_;
    }

    const scalar& getmu()   const
    {
        return mu_;
    }

    const scalar& getAdhN()   const
    {
        return adhN_;
    }

    const scalar& getAdhEqui()   const
    {
        return adhEqui_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
