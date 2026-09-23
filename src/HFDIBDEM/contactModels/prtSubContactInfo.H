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
    class for storage of contact variables

SourceFiles
    prtSubContactInfo.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
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
#include "demTimeStepInfo.H"


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

    vector& getcLVec()
    {
        return cLVec_;
    }

    vector& gettLVec()
    {
        return tLVec_;
    }

    vector getFNe();

    vector getFA();

    vector getFNd();

    // tangential force of the current sub-step; the coulomb
    // ceiling FtCeil (= mu*mag(FN) of this sub-contact) caps the
    // APPLIED force and the stored state FtPrev_ alike, so the
    // spring stops stretching at the sliding limit
    vector getFt(scalar deltaT, scalar FtCeil = GREAT);

    Tuple2<forces,forces>& getOutForce()
    {
        return outForce_;
    }

    void setVMInfo(boundBox& bBox, scalar subVolumeV, scalar emptyScale = 1);

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

    // stability step of this contact from the cached contact state
    // (the contact class of the adaptive DEM time-step estimator).
    // the governing stiffness is taken from the implemented force
    // laws with the tangential branch scaled by tangentialFactor
    // (1 = full conservatism). damping sits inside the stepDEM
    // safety factor.
    scalar getPairDtCrit(scalar tangentialFactor = 1.0)
    {
        return demTimeStepInfo::pairDtCrit
        (
            demTimeStepInfo::effK
            (
                physicalProperties_.aY_,
                prtCntVars_.contactVolume_,
                physicalProperties_.aG_,
                prtCntVars_.contactArea_,
                Lc_,
                tangentialFactor
            ),
            physicalProperties_.reduceM_
        );
    }

    // rotational stability step: the tangential force on the cached
    // lever arms (cLVec_ for body c, tLVec_ for body t - the full
    // length is the conservative rPerp)
    scalar getPairDtCritRot
    (
        const scalar IeffC,
        const scalar IeffT,
        scalar tangentialFactor = 1.0,
        scalar rotationalFactor = 1.0
    )
    {
        const scalar kT
        (
            rotationalFactor*tangentialFactor
                *demTimeStepInfo::contactKT
                (
                    physicalProperties_.aG_,
                    prtCntVars_.contactArea_,
                    Lc_
                )
        );

        scalar dtPair(getPairDtCrit(tangentialFactor));

        dtPair = min
        (
            dtPair,
            demTimeStepInfo::rotationalDtCrit(kT, mag(cLVec_), IeffC)
        );
        dtPair = min
        (
            dtPair,
            demTimeStepInfo::rotationalDtCrit(kT, mag(tLVec_), IeffT)
        );

        return dtPair;
    }

    void syncData();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
