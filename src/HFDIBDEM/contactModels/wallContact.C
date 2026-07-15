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

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "wallContact.H"
#include "wallMatInfo.H"

#include "virtualMeshLevel.H"
#include "wallPlaneInfo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace contactModel
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//---------------------------------------------------------------------------//
bool detectWallContact(
    const fvMesh&   mesh,
    ibContactClass& ibClass,
    wallContactInfo& wallCntInfo
)
{
    ibClass.setWallContact(false);

    if(ibClass.getGeomModel().getcType() == sphere)
    {
        return detectWallContact_Sphere
        (
            mesh,
            ibClass,
            wallCntInfo
        );
    }
    else
    if(ibClass.getGeomModel().getcType() == cluster)
    {
        return detectWallContact_Cluster
        (
            mesh,
            ibClass,
            wallCntInfo
        );
    }
    else
    {
        return detectWallContact_ArbShape
        (
            mesh,
            ibClass,
            wallCntInfo
        );
    }
}
//---------------------------------------------------------------------------/
bool detectWallContact_ArbShape(
    const fvMesh&   mesh,
    ibContactClass& ibClass,
    wallContactInfo& wallCntInfo
)
{
    bool isContact(false);
    if(wallCntInfo.detectWallContact())
    {
        // wallCntInfo.findContactAreas();
        isContact = true;
    }
    return(isContact);
}
//---------------------------------------------------------------------------//
bool detectWallContact_Sphere(
    const fvMesh&   mesh,
    ibContactClass& ibClass,
    wallContactInfo& wallCntInfo
)
{
    if (wallCntInfo.detectWallContact())
    {
        vector cCenter(wallCntInfo.getcClass().getGeomModel().getCoM());
        List<string>& contactPatches = wallCntInfo.getContactPatches();
        List<string> contactPatchesTmp;
        forAll(contactPatches, patchI)
        {
            List<vector> planeInfo = wallPlaneInfo::getWallPlaneInfo()[contactPatches[patchI]];
            plane p(planeInfo[1], planeInfo[0]);
            point nearestPoint = p.nearestPoint(cCenter);
            if(mag(cCenter - nearestPoint)-wallCntInfo.getcClass().getGeomModel().getDC()/2 > 0)
            {
                continue;
            }
            contactPatchesTmp.append(contactPatches[patchI]);
        }
        contactPatches = contactPatchesTmp;
        
        if(contactPatches.size() == 0)
        {
            return false;
        }
        wallCntInfo.getWallSCList().emplace_back(
            std::make_shared<wallSubContactInfo>(
                List<Tuple2<point,boundBox>>(),
                List<Tuple2<point,boundBox>>(),
                wallCntInfo.getContactPatches(),
                List<Tuple2<point,boundBox>>(),
                wallCntInfo.getWallMeanPars(),
                ibClass.getGeomModel().getBounds(),
                wallCntInfo.getBodyId()
            )
        );

        return true;
    }

    return false;
}
//---------------------------------------------------------------------------//
bool detectWallContact_Cluster(
    const fvMesh&   mesh,
    ibContactClass& ibClass,
    wallContactInfo& wallCntInfo
)
{
    wallCntInfo.getContactPatches().clear();

    periodicBody& cCluster = dynamic_cast<periodicBody&>(ibClass.getGeomModel());
    std::vector<std::shared_ptr<geomModel>> cBodies = cCluster.getClusterBodies();

    for(std::shared_ptr<geomModel>& cgModel : cBodies)
    {
        ibContactClass cIbClassI(
            cgModel,
            ibClass.getMatInfo().getMaterial()
        );

        ibContactVars cIbVars(
            wallCntInfo.getBodyId(),
            wallCntInfo.getcVars().Vel_,
            wallCntInfo.getcVars().omega_,
            wallCntInfo.getcVars().Axis_,
            cgModel->getM0(),
            cgModel->getM(),
            cgModel->getRhoS()
        );

        wallContactInfo tmpWallCntI(
            cIbClassI,
            cIbVars
        );

        if (detectWallContact(
            mesh,
            cIbClassI,
            tmpWallCntI
        ))
        {
            wallCntInfo.getWallSCList().insert(
                std::end(wallCntInfo.getWallSCList()),
                std::begin(tmpWallCntI.getWallSCList()),
                std::end(tmpWallCntI.getWallSCList())
            );

            wallCntInfo.getContactPatches().append(
                tmpWallCntI.getContactPatches()
            );

            cIbClassI.setWallContact(true);
            cIbClassI.inContactWithStatic(true);
        }

        if(cIbClassI.checkWallContact())
        {
            return true;
        }
    }

    return false;
}
//---------------------------------------------------------------------------//
void getWallContactVars(
    const fvMesh&   mesh,
    wallContactInfo& wallCntInfo,
    const scalar deltaT,
    wallSubContactInfo& sWC
)
{
    if (wallCntInfo.getcClass().getGeomModel().getcType() == sphere)
    {
        getWallContactVars_Sphere
        (
            mesh,
            wallCntInfo,
            deltaT,
            sWC
        );
    }
    else if(wallCntInfo.getcClass().getGeomModel().getcType() == cluster)
    {
        getWallContactVars_Cluster
        (
            mesh,
            wallCntInfo,
            deltaT,
            sWC
        );
    }
    else
    {
        getWallContactVars_ArbShape
        (
            mesh,
            wallCntInfo,
            deltaT,
            sWC
        );
    }
}
//---------------------------------------------------------------------------//
void getWallContactVars_ArbShape(
    const fvMesh&   mesh,
    wallContactInfo& wallCntInfo,
    const scalar deltaT,
    wallSubContactInfo& sCW
)
{
    scalar intersectVolume(0);
    vector contactCenter(vector::zero);
    scalar contactArea(0);
    vector contactNormal(vector::zero);

    autoPtr<DynamicVectorList> contactCenters(
        new DynamicVectorList);

    autoPtr<DynamicVectorList> contactPlaneCenters(
        new DynamicVectorList);

    autoPtr<DynamicScalarList> contactAreas(
        new DynamicScalarList);

    label vMContactInfoSize = sCW.getVMContactSize();
    label vMPlaneInfoSize  = sCW.getVMPlaneSize();
    const List<Tuple2<point,boundBox>>& sCInternalInfo = sCW.getInternalElements();
    const List<string>& contactPatches = sCW.getContactPatches();

    for(label i = 0; i< vMContactInfoSize; i++)
    {
        autoPtr<virtualMeshWallInfo> vmWInfo = sCW.getVMContactInfo(i);
        if (!vmWInfo.valid())
        {
            continue;
        }

        virtualMeshWall virtMeshWall(
            vmWInfo(),
            wallCntInfo.getcClass().getGeomModel()
        );

        if(virtMeshWall.detectFirstContactPoint())
        {
            intersectVolume += virtMeshWall.evaluateContact();
            contactCenters().append(virtMeshWall.getContactCenter());
        }
    }

    forAll(sCInternalInfo,sCII)
    {
        intersectVolume += sCInternalInfo[sCII].second().volume();
        contactCenters().append(sCInternalInfo[sCII].first());
    }

    if(intersectVolume>0)
    {
        for(label i = 0; i< vMPlaneInfoSize; i++)
        {
            autoPtr<virtualMeshWallInfo> vmWInfo = sCW.getVMPlaneInfo(i);
            if (!vmWInfo.valid())
            {
                continue;
            }

            autoPtr<virtualMeshWall> virtMeshPlane(new virtualMeshWall(
                vmWInfo(),
                wallCntInfo.getcClass().getGeomModel()
            ));

            if(virtMeshPlane->detectFirstFaceContactPoint())
            {
                scalar contactAreaLoc = (virtMeshPlane->evaluateContact()/vmWInfo->getSVVolume())*(pow(vmWInfo->getSVVolume(),2.0/3));
                contactAreas().append(contactAreaLoc);
                contactPlaneCenters().append(virtMeshPlane->getContactCenter());
            }
            else
            {
                contactAreas().append((intersectVolume/vmWInfo->getSVVolume())*(pow(vmWInfo->getSVVolume(),2.0/3)));
                contactPlaneCenters().append(virtMeshPlane->getContactCenter());
            }
        }

        forAll(contactCenters(),cC)
        {
            contactCenter += contactCenters()[cC];
        }
        // Pout << "contactCenter " << contactCenter.size() << endl;
        contactCenter /= contactCenters().size();
        // Pout << "Survived #0 " << endl;
        // Pout << "contactAreas() " << contactAreas() << endl;
        forAll(contactAreas(),cA)
        {
            contactArea += contactAreas()[cA];
        }

        if(contactArea == 0)
        {
            return;
        }
        // Pout << "contactArea " << contactArea << endl;
        // Pout << "Survived #0.5 " << endl;
        forAll(contactPatches,cP)
        {
            contactNormal -= wallPlaneInfo::getWallPlaneInfo()[contactPatches[cP]][0]*contactAreas()[cP];
        }
        // Pout << "contactNormal " << mag(contactNormal) << endl;
        contactNormal /=mag(contactNormal);
        // Pout << "contactNormal " << contactNormal << endl;
        wallCntInfo.getcClass().setWallContact(true);
        wallCntInfo.getcClass().inContactWithStatic(true);
        // Pout << "Survived #1 " << endl;
        wallContactVars& wallCntVars = sCW.getWallCntVars();
        // Pout << "Survived #2 " << endl;
        wallCntVars.contactCenter_ = contactCenter;
        wallCntVars.contactArea_   = contactArea;
        wallCntVars.contactVolume_ = intersectVolume;
        wallCntVars.contactNormal_ = contactNormal;
        // Pout << "Survived #3 " << endl;
        wallCntVars.setMeanCntPars_Plane
        (
            contactAreas(),
            contactPatches,
            wallCntInfo.getWallMeanPars()
        );
        // Pout << "Survived #4 " << endl;
    }
}
// //---------------------------------------------------------------------------//
void getWallContactVars_Sphere(
    const fvMesh&   mesh,
    wallContactInfo& wallCntInfo,
    const scalar deltaT,
    wallSubContactInfo& sCW
)
{
    vector cCenter(wallCntInfo.getcClass().getGeomModel().getCoM());

    List<string>& contactPatches = wallCntInfo.getContactPatches();

    List<wallContactVars> wallCntVarsList;
    forAll(contactPatches, patchI)
    {
        List<vector> planeInfo = wallPlaneInfo::getWallPlaneInfo()[contactPatches[patchI]];
        plane p(planeInfo[1], planeInfo[0]);
        point nearestPoint = p.nearestPoint(cCenter);

        vector nVec = (cCenter - nearestPoint)/mag(cCenter - nearestPoint);

        wallCntVarsList.append(wallContactVars());

        wallCntVarsList.last().contactCenter_ = nearestPoint;
        wallCntVarsList.last().contactNormal_ = nVec;
        wallCntVarsList.last().contactArea_ = sphereContactArea
        (
            mesh,
            wallCntInfo.getcClass(),
            wallCntInfo.getcVars(),
            p
        );

        wallCntVarsList.last().contactVolume_ = getInterVolume_Sphere
        (
            mesh,
            wallCntInfo.getcClass(),
            wallCntInfo.getcVars(),
            p
        );
    }

    wallContactVars& wallCntVars = sCW.getWallCntVars();

    List<scalar> contactAreas;

    for (int i = 0; i < wallCntVarsList.size(); ++i)
    {
        wallCntVars.contactCenter_ += wallCntVarsList[i].contactCenter_ * wallCntVarsList[i].contactVolume_;
        wallCntVars.contactNormal_ += wallCntVarsList[i].contactNormal_ * wallCntVarsList[i].contactVolume_;
        wallCntVars.contactArea_ += wallCntVarsList[i].contactArea_;
        contactAreas.append(wallCntVarsList[i].contactArea_);
        wallCntVars.contactVolume_ += wallCntVarsList[i].contactVolume_;
    }

    wallCntVars.contactCenter_ /= wallCntVars.contactVolume_;
    wallCntVars.contactNormal_ /= wallCntVars.contactVolume_;

    wallCntVars.setMeanCntPars_Plane
    (
        contactAreas,
        contactPatches,
        wallCntInfo.getWallMeanPars()
    );

    if(wallCntVars.contactVolume_ > 0)
    {
        wallCntInfo.getcClass().setWallContact(true);
        wallCntInfo.getcClass().inContactWithStatic(true);
    }
}
//---------------------------------------------------------------------------//
void getWallContactVars_Cluster(
    const fvMesh&   mesh,
    wallContactInfo& wallCntInfo,
    const scalar deltaT,
    wallSubContactInfo& sWC
)
{
    periodicBody& cCluster = dynamic_cast<periodicBody&>(wallCntInfo.getcClass().getGeomModel());
    std::vector<std::shared_ptr<geomModel>> cBodies = cCluster.getClusterBodies();

    for(std::shared_ptr<geomModel>& cgModel : cBodies)
    {
        if(!sWC.getsWCBB().overlaps(cgModel->getBounds()))
        {
            continue;
        }

        ibContactClass cIbClassI(
            cgModel,
            wallCntInfo.getcClass().getMatInfo().getMaterial()
        );

        wallContactInfo cWallCntI(
            cIbClassI,
            wallCntInfo.getcVars()
        );

        cWallCntI.getContactPatches() = wallCntInfo.getContactPatches();

        getWallContactVars(
            mesh,
            cWallCntI,
            deltaT,
            sWC
        );
    }
}
//---------------------------------------------------------------------------//
bool solveWallContact
(
    const fvMesh&   mesh,
    wallContactInfo& wallCntInfo,
    scalar deltaT,
    wallSubContactInfo& sCI
)
{
    getWallContactVars(
        mesh,
        wallCntInfo,
        deltaT,
        sCI
    );

    vector outF = vector::zero;
    vector cLVecOut = vector::zero;

    wallContactVars& wallCntVar = sCI.getWallCntVars();

    sCI.evalVariables(wallCntVar,wallCntInfo.getcClass(),wallCntInfo.getcVars());

    if(wallCntVar.contactVolume_ == 0)
    {
        sCI.getOutForce().F = vector::zero;
        sCI.getOutForce().T = vector::zero;
        return false;
    }

    InfoH << parallelDEM_Info << "-- Detected Particle-wall contact: -- body "
        << sCI.getBodyId() << endl;
    InfoH << parallelDEM_Info << "-- body "<< sCI.getBodyId() <<"  linear velocity:"
        << wallCntInfo.getcVars().Vel_ << " magnitude: " << mag(wallCntInfo.getcVars().Vel_) <<endl;
    InfoH << parallelDEM_Info << "-- body "<< sCI.getBodyId() <<"  angular velocity:"
        << wallCntInfo.getcVars().omega_ << " magnitude: " << mag(wallCntInfo.getcVars().omega_) <<endl;
    InfoH << parallelDEM_Info << "-- Particle-wall body "<< sCI.getBodyId() <<" contact center "
        << wallCntVar.contactCenter_ << endl;
    InfoH << parallelDEM_Info << "-- Particle-wall body "<< sCI.getBodyId() <<" contact normal "
        << wallCntVar.contactNormal_ << endl;
    InfoH << parallelDEM_Info << "-- Particle-wall body "<< sCI.getBodyId() <<" contact volume "
        << wallCntVar.contactVolume_ << endl;
    InfoH << parallelDEM_Info << "-- Particle-wall body "<< sCI.getBodyId() <<" contact area "
        << wallCntVar.contactArea_ << endl;

    vector F = sCI.getFNe(wallCntVar);
    InfoH << parallelDEM_Info << "-- Particle-wall body "<< sCI.getBodyId() <<" contact FNe " << F << endl;

    vector FNd = sCI.getFNd(wallCntVar);
    InfoH << parallelDEM_Info << "-- Particle-wall body "<< sCI.getBodyId() <<" contact FNd " << FNd << endl;

    if ((F & FNd) < 0 && mag(FNd) > mag(F))
    {
        FNd *= mag(F) / mag(FNd);
        InfoH << parallelDEM_Info << "FNd was Clipped to "<< FNd << endl;
    }

    F += FNd;
    InfoH << parallelDEM_Info << "-- Particle-wall body "<< sCI.getBodyId() <<" contact FN " << F << endl;

    vector Ft = sCI.getFt(wallCntVar, deltaT);
    InfoH << parallelDEM_Info << "-- Particle-wall body "<< sCI.getBodyId() <<" contact Ft " << Ft << endl;

    if (mag(Ft) > sCI.getMu(wallCntVar) * mag(F))
    {
        Ft *= sCI.getMu(wallCntVar) * mag(F) / mag(Ft);
    }
    InfoH << parallelDEM_Info << "-- Particle-wall body "<< sCI.getBodyId() <<" contact Ft clamped" << Ft << endl;
    F += Ft;

    vector FA = sCI.getFA(wallCntVar);
    InfoH << parallelDEM_Info << "-- Particle-wall body "<< sCI.getBodyId() <<" contact FA " << FA << endl;
    F -= FA;

    outF += F;

    cLVecOut += wallCntVar.lVec_;

    sCI.getOutForce().F = outF;
    sCI.getOutForce().T = cLVecOut ^  outF;
    return true;
}
// //---------------------------------------------------------------------------//
scalar getInterVolume_Sphere(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactVars& cVars,
    plane& cPlane
)
{
    scalar cRadius(cClass.getGeomModel().getDC() / 2);
    vector cCenter(cClass.getGeomModel().getCoM());

    scalar dist = cPlane.distance(cCenter);
    scalar xH = cRadius - dist;

    if(case3D)
    {
        return (Foam::constant::mathematical::pi*sqr(xH)*(3*cRadius - xH)) / 3;
    }
    else
    {
        boundBox meshBounds = mesh.bounds();
        scalar emptyLength = meshBounds.max()[emptyDim]
                            - meshBounds.min()[emptyDim];

        scalar cCirSeg = sqr(cRadius)*acos((dist)/cRadius)
                         - (dist)*sqrt(sqr(cRadius) - sqr(dist));
        return cCirSeg*emptyLength;
    }
}
//---------------------------------------------------------------------------//
scalar sphereContactArea
(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactVars& cVars,
    plane& cPlane
)
{
    scalar cRadius(cClass.getGeomModel().getDC() / 2);
    vector cCenter(cClass.getGeomModel().getCoM());

    scalar dist = cPlane.distance(cCenter);
    scalar contactRad = sqrt(sqr(cRadius) - sqr(dist));

    if(case3D)
    {
        return Foam::constant::mathematical::pi*sqr(contactRad);
    }
    else
    {
        boundBox meshBounds = mesh.bounds();
        scalar emptyLength = meshBounds.max()[emptyDim]
                            - meshBounds.min()[emptyDim];

        return contactRad*emptyLength;
    }
}
//---------------------------------------------------------------------------//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace contactModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
