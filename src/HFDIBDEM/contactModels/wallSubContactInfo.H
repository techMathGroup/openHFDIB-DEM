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
    wallSubContactInfo.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef wallSubContactInfo_H
#define wallSubContactInfo_H

#include "dictionary.H"
#include "fvCFD.H"

#include "geomModel.H"
#include "ibContactVars.H"
#include "ibContactClass.H"
#include "wallContactVars.H"
#include "materialInfo.H"
#include "physicalProperties.H"
#include "forces.H"
#include "demTimeStepInfo.H"
#include "spectatorMesh.H"
#include "virtualMeshWall.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class wallSubContactInfo Declaration
\*---------------------------------------------------------------------------*/

class wallSubContactInfo
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;
    // list of collision pathes 
    List<string> contactPatches_;

    //bool list for internal SM elements
    List<Tuple2<point,boundBox>> internalBBData_;

    // Mean model parameters
    HashTable<physicalProperties,string,Hash<string>> wallMeanPars_;

    boundBox BB_;

    label bodyId_;

    wallContactVars wallCntvar_;

    // per-patch contact variables (the inclusion-exclusion
    // attribution of a multi-patch sub-contact, or the aggregate
    // as one entry otherwise); index-aligned with
    // contactPatches_
    List<wallContactVars> patchVars_;

    // per-patch tangential springs restored (by patch name)
    // before the patch variables are built; consumed by the
    // patchVars fill, which happens inside solveWallContact
    // (wallContactInfo::restoreWallFtHistory runs earlier in
    // the sub-step)
    HashTable<vector, string, Hash<string>> patchFtPrev_;

    // creating member of a cluster sub-contact (null on plain
    // bodies); set at detection so the evaluation runs the
    // member exclusively
    std::shared_ptr<geomModel> creatorModel_ = nullptr;

    scalar reduceM_ = 0;    // Reduced mass

    forces outForce_; // Out variable for force

    List<autoPtr<virtualMeshWallInfo>> vmWInfoList_; 

    List<autoPtr<virtualMeshWallInfo>> vmPlaneInfoList_;

    bool contactResolved_ = false;

    public:

    // Constructors ------------------------------------------------
    wallSubContactInfo
    (
    List<Tuple2<point,boundBox>> contactBBData,
    List<Tuple2<point,boundBox>> planeBBData,
    List<string> contactPatches,
    List<Tuple2<point,boundBox>> internalBBData,
    HashTable<physicalProperties,string,Hash<string>> wallMeanPars,
    boundBox BB,
    label bodyId
    );

    // Destructors -------------------------------------------------
    ~wallSubContactInfo();

    label getVMContactSize()
    {
        return vmWInfoList_.size();
    }
    
    label getVMPlaneSize()
    {
        return vmPlaneInfoList_.size();
    }

    autoPtr<virtualMeshWallInfo>& getVMContactInfo
    (
        label ID
    );

    autoPtr<virtualMeshWallInfo>& getVMPlaneInfo
    (
        label ID
    );

    wallContactVars& getWallCntVars()
    {
        return wallCntvar_;
    }

    // per-patch force variables of every sub-contact; empty when
    // the volume gate rejected the sub-contact this sub-step
    List<wallContactVars>& getPatchCntVars()
    {
        return patchVars_;
    }

    // staged per-patch springs from the tangential history
    HashTable<vector, string, Hash<string>>& getPatchFtPrev()
    {
        return patchFtPrev_;
    }

    // creating member of a cluster sub-contact (null on plain
    // bodies); the cluster evaluation runs this member only
    std::shared_ptr<geomModel>& getCreatorModel()
    {
        return creatorModel_;
    }

    void setCreatorModel(std::shared_ptr<geomModel> creator)
    {
        creatorModel_ = std::move(creator);
    }

    scalar& getReduceM()
    {
        return reduceM_;
    }

    const boundBox& getBB() const
    {
        return BB_;
    }

    scalar getMu(wallContactVars& wallCntvar)
    {
        return wallCntvar.getMeanCntPar().aMu_;
    }

    vector getLVec
    (
        wallContactVars& wallCntvar,
        ibContactClass ibCClass
    );

    vector getVeli
    (
        wallContactVars& wallCntvar,
        ibContactVars& cVars
    );

    void evalVariables
    (
        wallContactVars& wallCntvar,
        ibContactClass& ibCClass,
        ibContactVars& cVars
    );

    vector getFNe(wallContactVars& wallCntvar);

    vector getFA(wallContactVars& wallCntvar);

    vector getFNd(wallContactVars& wallCntvar);

    // tangential force of the current sub-step; the coulomb
    // ceiling FtCeil (= mu*mag(FN) of this sub-contact) caps the
    // APPLIED force and the stored state FtPrev_ alike, so the
    // spring stops stretching at the sliding limit
    vector getFt
    (
        wallContactVars& wallCntvar,
        scalar deltaT,
        scalar FtCeil = GREAT
    );

    // stability step of this wall contact from the cached contact
    // state (the contact class of the adaptive DEM time-step
    // estimator). the wall mass is infinite, so the reduced mass
    // is the full body mass M0 - same value reduceM_ now carries
    // for the force laws. tangentialFactor scales the tangential
    // branch (1 = full conservatism); damping sits inside the
    // stepDEM safety factor.
    scalar getPairDtCrit
    (
        wallContactVars& wallCntvar,
        ibContactVars& cVars,
        scalar tangentialFactor = 1.0
    )
    {
        return demTimeStepInfo::pairDtCrit
        (
            demTimeStepInfo::effK
            (
                wallCntvar.getMeanCntPar().aY_,
                wallCntvar.contactVolume_,
                wallCntvar.getMeanCntPar().aG_,
                wallCntvar.contactArea_,
                wallCntvar.Lc_,
                tangentialFactor
            ),
            cVars.M0_
        );
    }

    // rotational stability step: the tangential
    // force on the cached lever arm lVec_ acts on the body's
    // lambda_min inertia (the wall contributes no rotational mode).
    scalar getPairDtCritRot
    (
        wallContactVars& wallCntvar,
        ibContactVars& cVars,
        scalar IeffC,
        scalar tangentialFactor = 1.0,
        scalar rotationalFactor = 1.0
    )
    {
        const scalar kT
        (
            rotationalFactor*tangentialFactor
                *demTimeStepInfo::contactKT
                (
                    wallCntvar.getMeanCntPar().aG_,
                    wallCntvar.contactArea_,
                    wallCntvar.Lc_
                )
        );

        return min
        (
            getPairDtCrit(wallCntvar, cVars, tangentialFactor),
            demTimeStepInfo::rotationalDtCrit(kT, mag(wallCntvar.lVec_), IeffC)
        );
    }

    forces& getOutForce()
    {
        return outForce_;
    }

    void setResolvedContact(bool resolved)
    {
        contactResolved_ = resolved;
    }

    bool getContactResolved()
    {
        return contactResolved_;
    }

    bool isContactPatch(string patchName)
    {
        return wallMeanPars_.found(patchName);
    }

    HashTable<physicalProperties,string,Hash<string>>& getWallMeanPars()
    {
        return wallMeanPars_;
    }

    List<autoPtr<virtualMeshWallInfo>>& getVMContactInfo()
    {
        return vmWInfoList_;
    }

    List<autoPtr<virtualMeshWallInfo>>& getVMPlaneInfo()
    {
        return vmPlaneInfoList_;
    }

    const List<Tuple2<point,boundBox>>& getInternalElements()
    {
        return internalBBData_;
    }

    const List<string>& getContactPatches()
    {
        return contactPatches_;
    }

    void syncData();

    void syncContactResolve();

    label& getBodyId()
    {
        return bodyId_;
    }

    boundBox& getsWCBB()
    {
        return BB_;
    }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //