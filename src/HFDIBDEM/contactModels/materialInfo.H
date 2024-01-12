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
    materialInfo.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef materialInfo_H
#define materialInfo_H

#include "dictionary.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class materialInfo Declaration
\*---------------------------------------------------------------------------*/

class materialInfo
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

    protected:

        string material_;

        // Young modulus
        scalar  Y_;

        // poisson ration
        scalar  nu_;

        // friction coefficient
        scalar mu_;

        // adhesion coefficient for normal contact force
        scalar adhN_;

        //coeffitient of restitution for normal force
        scalar eps_;

    public:

    // Constructors ------------------------------------------------
    materialInfo
    (
        string material,
        scalar Y,
        scalar nu,
        scalar mu,
        scalar adhN,
        scalar eps
    );

    materialInfo
    (
        const materialInfo& mi
    );

    // Destructors -------------------------------------------------
    ~materialInfo();

    const string& getMaterial() const
    {
        return material_;
    }

    const scalar& getY()   const
    {
        return Y_;
    }

    const scalar& getNu()   const
    {
        return nu_;
    }

    const scalar& getMu()   const
    {
        return mu_;
    }

    const scalar& getAdhN()   const
    {
        return adhN_;
    }

    const scalar& getEps()   const
    {
        return eps_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
