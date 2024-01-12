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
    prtSubContactInfo.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef prtSubContactInfo_H
#define prtSubContactInfo_H

#include "dictionary.H"
#include "fvCFD.H"

#include "geomModel.H"
#include "ibContactClass.H"
#include "ibContactVars.H"
#include "prtContactVars.H"
#include "virtualMesh.H"
#include "physicalProperties.H"
#include "forces.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class prtSubContactInfo Declaration
\*---------------------------------------------------------------------------*/

class prtSubContactInfo
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

    Tuple2<label,label> contactPair_;

    physicalProperties physicalProperties_;

    std::shared_ptr<virtualMeshInfo> vmInfo_;

    prtContactVars prtCntVars_;

    vector cLVec_;       // Direction to contact center for body c
    vector tLVec_;       // Direction to contact center for body t

    vector cVeli_;       // Velocity at contact center for body c
    vector tVeli_;       // Velocity at contact center for body t

    scalar Vn_;          // Relative velocity in the normal direction
    scalar Lc_;          // Characteristic length

    Tuple2<forces,forces> outForce_; // Out variable for force

    vector FtPrev_ = vector::zero;     // Tangential force from previous time-step

    bool contactResolved_ = false;

    vector getVeli(ibContactVars& cVars, vector& lVec);

    public:

    // Constructors ------------------------------------------------
    prtSubContactInfo
    (
        const Tuple2<label,label>& contactPair,
        const physicalProperties& physicalProperties
    );

    // Destructors -------------------------------------------------
    ~prtSubContactInfo();

    const Tuple2<label, label>& getCPair()
    {
        return contactPair_;
    }

    physicalProperties& getPhysicalProperties()
    {
        return physicalProperties_;
    }

    prtContactVars& getprtCntVars()
    {
        return prtCntVars_;
    }

    void evalVariables(
        ibContactClass& cCoM,
        ibContactClass& tCoM,
        ibContactVars& cVars,
        ibContactVars& tVars
    );

    vector getcLVec()
    {
        return cLVec_;
    }

    vector gettLVec()
    {
        return tLVec_;
    }

    vector getFNe();

    vector getFA();

    vector getFNd();

    vector getFt(scalar deltaT);

    Tuple2<forces,forces>& getOutForce()
    {
        return outForce_;
    }

    void setVMInfo(boundBox& bBox, scalar subVolumeV);

    void setVMInfo(const virtualMeshInfo& vmInfo);

    std::shared_ptr<virtualMeshInfo>& getVMInfo();

    void clearOutForces()
    {
        outForce_.first().F = vector::zero;
        outForce_.first().T = vector::zero;
        outForce_.second().F = vector::zero;
        outForce_.second().T = vector::zero;
    }

    bool contactResolved()
    {
        bool reduceBool = contactResolved_;
        reduce(reduceBool,orOp<bool>());

        return reduceBool;
    }

    void setResolvedContact(bool resolved)
    {
        contactResolved_ = resolved;
    }

    scalar& getLc()
    {
        return Lc_;
    }

    void syncData();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
