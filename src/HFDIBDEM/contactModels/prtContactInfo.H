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
    class for storage of contact variables

SourceFiles
    prtContactInfo.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef prtContactInfo_H
#define prtContactInfo_H

#include "dictionary.H"
#include "fvCFD.H"

#include "geomModel.H"
#include "ibContactVars.H"
#include "ibContactClass.H"
#include "prtContactVars.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class prtContactInfo Declaration
\*---------------------------------------------------------------------------*/

class prtContactInfo
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

    protected:

    ibContactClass& cIbContactClass_;
    ibContactVars& cContactVars_;

    ibContactClass& tIbContactClass_;
    ibContactVars& tContactVars_;

    Tuple2<label,label> contactPair_;

    prtContactVars prtCntVars_;

    // Mean model parameters
    scalar aY_;         // Contact Young modulus
    scalar aG_;         // Contact shear modulus
    scalar aGammaN_;    // Contact normal gamma
    scalar aGammat_;    // Contact tangential gamma
    scalar aMu_;        // Contact frictional coef
    scalar maxAdhN_;    // Contact adhesive coef
    scalar curAdhN_;    // Current contact adhesive coef

    scalar reduceM_;    // Reduced mass

    vector cLVec_;       // Direction to contact center for body c
    vector tLVec_;       // Direction to contact center for body t

    vector cVeli_;       // Velocity at contact center for body c
    vector tVeli_;       // Velocity at contact center for body t

    scalar Vn_;          // Relative velocity in the normal direction
    scalar Lc_;          // Characteristic length

    vector FtPrev_;     // Tangential force from previous time-step

    Tuple2<forces,forces> outForce_; // Out variable for force

    public:

    // Constructors ------------------------------------------------
    prtContactInfo
    (
        ibContactClass& cClass,
        ibContactVars& cVars,
        ibContactClass& tClass,
        ibContactVars& tVars,
        HashTable<scalar,string,Hash<string>>& matInterAdh
    );

    // Destructors -------------------------------------------------
    ~prtContactInfo();

    ibContactClass& getcClass()
    {
        return cIbContactClass_;
    }

    ibContactClass& gettClass()
    {
        return tIbContactClass_;
    }

    Tuple2<label,label>& getcPair()
    {
        return contactPair_;
    }

    prtContactVars& getprtCntVars()
    {
        return prtCntVars_;
    }

    scalar getMu()
    {
        return aMu_;
    }

    vector getcLVec()
    {
        return cLVec_;
    }

    vector gettLVec()
    {
        return tLVec_;
    }

    vector getLVec(ibContactClass& contactClass);

    vector getVeli(ibContactVars& cVars, vector& lVec);

    void evalVariables();

    vector getFNe();

    vector getFA();

    vector getFNd();

    vector getFt(scalar deltaT);

    Tuple2<forces,forces>& getOutForce()
    {
        return outForce_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
