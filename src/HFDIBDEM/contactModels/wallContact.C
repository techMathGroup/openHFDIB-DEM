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

    // if(ibClass.getGeomModel().getcType() == sphere)
    // {
    //     return detectWallContact_Sphere
    //     (
    //         mesh,
    //         ibClass
    //     );
    // }
    // else
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
        wallCntInfo.findContactAreas();
        isContact = true;
    }
    return(isContact);
}
//---------------------------------------------------------------------------//
// bool detectWallContact_Sphere(
//     const fvMesh&   mesh,
//     ibContactClass& ibClass
// )
// {
//     scalar cRadius(ibClass.getGeomModel().getDC() / 2);
//     vector cCenter(ibClass.getGeomModel().getCoM());

//     label nCells = mesh.nCells();

//     List<DynamicLabelList>& surfCells(ibClass.getSurfCells());
//     bool inContact = false;


//     // go through all surfCells and check if there is any surfCell whose face is a boundary face
//     forAll (surfCells[Pstream::myProcNo()],sCellI)
//     {
//         label cCell(surfCells[Pstream::myProcNo()][sCellI]);
//         if(cCell < nCells)
//         {
//             const labelList& cFaces = mesh.cells()[cCell];
//             forAll (cFaces,faceI)
//             {
//                 if (!mesh.isInternalFace(cFaces[faceI]))
//                 {
//                     // get reference to the patch which is in contact with IB. There is contact only if the patch is marked as a wall
//                     label facePatchId = mesh.boundaryMesh().whichPatch(cFaces[faceI]);
//                     const polyPatch& cPatch = mesh.boundaryMesh()[facePatchId];
//                     if (cPatch.type()=="wall")
//                     {
//                         vector nVec(-mesh.Sf()[cFaces[faceI]]/mag(mesh.Sf()[cFaces[faceI]]));

//                         if(!case3D)
//                         {
//                             if(mag(emptyDir - nVec) < 1e-4 || mag(emptyDir - (-1)*nVec) < 1e-4)
//                             {
//                                 continue;
//                             }
//                         }

//                         face currentFace = mesh.faces()[cFaces[faceI]];
//                         pointHit pointH = currentFace.nearestPoint(cCenter,mesh.points());
//                         if(mag(pointH.rawPoint() - cCenter) < cRadius)
//                         {
//                             inContact = true;
//                         }
//                     }
//                 }
//                 if(inContact)
//                     break;
//             }
//         }
//         if(inContact)
//             break;
//     }

//     reduce(inContact, orOp<bool>());

//     if(inContact)
//     {
//         return true;
//     }
//     return false;
// }
//---------------------------------------------------------------------------//
bool detectWallContact_Cluster(
    const fvMesh&   mesh,
    ibContactClass& ibClass,
    wallContactInfo& wallCntInfo
)
{
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
    if(wallCntInfo.getcClass().getGeomModel().getcType() == cluster)
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
    //run forAll for all surface Out Of Wall Contact elements get local contactCenters and intersectVolume            // InfoH << DEM_Info << " virtMeshWall.evaluateContact() " << virtMeshWall.evaluateContact() <<endl;
    // InfoH << DEM_Info << " -- WallContact detectVolume VM " << endl;
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

    //run forAll for all internal Out Of Wall Contact elements get local contactCenters and intersectVolume
    forAll(sCInternalInfo,sCII)
    {
        intersectVolume += sCInternalInfo[sCII].second().volume();
        contactCenters().append(sCInternalInfo[sCII].first());
    }
    //run forAll contactPlanes to get contact areas for the givenPatch
    InfoH << parallelDEM_Info << " -- VM prtId : "<< " intersectVolume : " <<  intersectVolume << endl;
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
            // InfoH << DEM_Info << " -- VM vmWInfo().BBox() " <<vmWInfo->getbBox() << endl;

            if(virtMeshPlane->detectFirstFaceContactPoint())
            {
                scalar contactAreaLoc = (virtMeshPlane->evaluateContact()/vmWInfo->getSVVolume())*(pow(vmWInfo->getSVVolume(),2.0/3));
                contactAreas().append(contactAreaLoc);
                contactPlaneCenters().append(virtMeshPlane->getContactCenter());
            }
        }

        forAll(contactCenters(),cC)
        {
            contactCenter += contactCenters()[cC];
        }
        contactCenter /= contactCenters().size();

        // InfoH << DEM_Info << " -- VM contactAreas() " << contactAreas() << endl;
        forAll(contactAreas(),cA)
        {
            contactArea += contactAreas()[cA];
        }

        forAll(contactPatches,cP)
        {
            contactNormal -= wallPlaneInfo::getWallPlaneInfo()[contactPatches[cP]][0]*contactAreas()[cP];
        }
        InfoH << parallelDEM_Info << " -- VM prtId : "<< sCW.getBodyId()<<" contactAreas() " << contactAreas() << endl;
        if (contactAreas().size() < SMALL)
        {
            InfoH << parallelDEM_Info << " -- VM CoM Of Failed Prt  : " << wallCntInfo.getcClass().getGeomModel().getCoM() << endl;
            InfoH << parallelDEM_Info << " -- VM CoM Of ContactArea : " << contactCenter << endl;
        }
        contactNormal /=mag(contactNormal);

        InfoH << parallelDEM_Info << " -- VM prtId : "<< sCW.getBodyId()<<" contactCenter_ : " << contactCenter << endl;
        InfoH << parallelDEM_Info << " -- VM prtId : "<< sCW.getBodyId()<<" contactArea_   : " << contactArea << endl;
        InfoH << parallelDEM_Info << " -- VM prtId : "<< sCW.getBodyId()<<" contactNormal_ : " << contactNormal << endl;

        wallCntInfo.getcClass().setWallContact(true);
        wallCntInfo.getcClass().inContactWithStatic(true);

        wallContactVars& wallCntVars = sCW.getWallCntVars();
        wallCntVars.contactCenter_ = contactCenter;
        wallCntVars.contactArea_   = contactArea;
        wallCntVars.contactVolume_ = intersectVolume;
        wallCntVars.contactNormal_ = contactNormal;

        wallCntVars.setMeanCntPars_Plane
        (
            contactAreas(),
            contactPatches,
            wallCntInfo.getWallMeanPars()
        );
    }
}
// //---------------------------------------------------------------------------//
// DynamicList<Tuple2<label,string>> getContactFacesArbShape
// (
//     const fvMesh&   mesh,
//     wallContactInfo& wallCntInfo,
//     label faceLabel,
//     DynamicList<Tuple2<label,string>>& checkedFaces,
//     HashTable<bool,label,Hash<label>>& pointTable,
//     DynamicVectorList& center,
//     DynamicVectorList& normal,
//     DynamicLabelList& centerPoints,
//     DynamicScalarList& area
// )
// {
//     DynamicList<Tuple2<label,string>> facesReturnList;

//     labelList faceEdges = mesh.faceEdges()[faceLabel];
//     forAll(faceEdges,edge)
//     {
//         labelList edgeFaces = mesh.edgeFaces()[faceEdges[edge]];
//         forAll(edgeFaces,faceI)
//         {
//             if (!mesh.isInternalFace(edgeFaces[faceI]))
//             {
//                 // get reference to the patch which is in contact with IB.
//                 // There is contact only if the patch is marked as a wall
//                 label facePatchId = mesh.boundaryMesh().whichPatch(edgeFaces[faceI]);
//                 const polyPatch& cPatch = mesh.boundaryMesh()[facePatchId];
//                 if (wallCntInfo.isContactPatch(cPatch.name()))
//                 {
//                     vector nVec(-mesh.Sf()[edgeFaces[faceI]]/mag(mesh.Sf()[edgeFaces[faceI]]));

//                     if(!case3D)
//                     {
//                         if(mag(emptyDir - nVec) < 1e-4 || mag(emptyDir - (-1)*nVec) < 1e-4)
//                             continue;
//                     }

//                     bool cont = false;

//                     forAll(checkedFaces,face)
//                     {
//                         if(checkedFaces[face].first() == edgeFaces[faceI])
//                         {
//                             cont = true;
//                             break;
//                         }
//                     }

//                     if(cont)
//                         continue;

//                     labelList facePoints(mesh.faces()[edgeFaces[faceI]]);//list of vertex indicies
//                     label numOfPoints = 0;

//                     forAll(facePoints,pointI)
//                     {
//                         if(pointTable.found(facePoints[pointI]))
//                         {
//                             if(pointTable[facePoints[pointI]])
//                             {
//                                 numOfPoints++;
//                                 center.last() += mesh.points()[facePoints[pointI]];
//                                 normal.last() += nVec;
//                             }
//                         }
//                         else
//                         {
//                             if(wallCntInfo.getcClass().getGeomModel().pointInside(mesh.points()[facePoints[pointI]]))
//                             {
//                                 numOfPoints++;
//                                 center.last() += mesh.points()[facePoints[pointI]];
//                                 normal.last() += nVec;
//                                 pointTable.insert(facePoints[pointI], true);
//                             }
//                             else
//                             {
//                                 pointTable.insert(facePoints[pointI], false);
//                             }
//                         }
//                     }

//                     if(numOfPoints > 0)
//                     {
//                         centerPoints.last() += numOfPoints;
//                         area.last() += mag(mesh.Sf()[edgeFaces[faceI]])
//                             *numOfPoints/facePoints.size();
//                         facesReturnList.append(Tuple2<label,string>
//                         (
//                             edgeFaces[faceI],
//                             cPatch.name()
//                         ));
//                     }
//                 }
//             }
//         }
//     }

//     return facesReturnList;
// }
//---------------------------------------------------------------------------//
// void getWallContactVars_Sphere(
//     const fvMesh&   mesh,
//     wallContactInfo& wallCntInfo,
//     const scalar deltaT,
//     label sWC
// )
// {
//     scalar cRadius(wallCntInfo.getcClass().getGeomModel().getDC() / 2);
//     vector cCenter(wallCntInfo.getcClass().getGeomModel().getCoM());

//     label nCells = mesh.nCells();

//     // Tuple: first() = face label; second() = patch name
//     List<DynamicList<Tuple2<label,string>>> contactFaces;

//     DynamicVectorList center(Pstream::nProcs(),vector::zero);
//     DynamicVectorList normal(Pstream::nProcs(),vector::zero);

//     List<DynamicLabelList>& surfCells(wallCntInfo.getcClass().getSurfCells());

//     // go through all surfCells and check if there is any surfCell whose face is a boundary face
//     forAll (surfCells[Pstream::myProcNo()],sCellI)
//     {
//         label cCell(surfCells[Pstream::myProcNo()][sCellI]);
//         if(cCell < nCells)
//         {
//             const labelList& cFaces = mesh.cells()[cCell];
//             forAll (cFaces,faceI)
//             {
//                 if (!mesh.isInternalFace(cFaces[faceI]))
//                 {
//                     // get reference to the patch which is in contact with IB. There is contact only if the patch is marked as a wall
//                     label facePatchId = mesh.boundaryMesh().whichPatch(cFaces[faceI]);
//                     const polyPatch& cPatch = mesh.boundaryMesh()[facePatchId];
//                     if (wallCntInfo.isContactPatch(cPatch.name()))
//                     {
//                         vector nVec(-mesh.Sf()[cFaces[faceI]]/mag(mesh.Sf()[cFaces[faceI]]));

//                         if(!case3D)
//                         {
//                             if(mag(emptyDir - nVec) < 1e-4 || mag(emptyDir - (-1)*nVec) < 1e-4)
//                             {
//                                 continue;
//                             }
//                         }

//                         bool cont = false;
//                         forAll(contactFaces,list)
//                         {
//                             forAll(contactFaces[list],face)
//                             {
//                                 if(contactFaces[list][face].first() == cFaces[faceI])
//                                 {
//                                     cont = true;
//                                     break;
//                                 }
//                             }
//                             if(cont)
//                                 break;
//                         }
//                         if(cont)
//                             continue;

//                         face currentFace = mesh.faces()[cFaces[faceI]];
//                         pointHit pointH = currentFace.nearestPoint(cCenter,mesh.points());
//                         if(mag(pointH.rawPoint() - cCenter) < cRadius)
//                         {
//                             DynamicList<Tuple2<label,string>> newContactFaces(
//                                 1,
//                                 Tuple2<label,string>(cFaces[faceI], cPatch.name())
//                             );

//                             labelList nextToCheck(1,cFaces[faceI]);

//                             while (nextToCheck.size() > 0)
//                             {
//                                 DynamicLabelList auxToCheck;

//                                 forAll (nextToCheck,faceToCheck)
//                                 {
//                                     DynamicList<Tuple2<label,string>> newfaces =
//                                     getContactFacesSphere(
//                                                 mesh,
//                                                 wallCntInfo,
//                                                 nextToCheck[faceToCheck],
//                                                 newContactFaces
//                                     );

//                                     forAll(newfaces,face)
//                                     {
//                                         newContactFaces.append(newfaces[face]);
//                                         auxToCheck.append(newfaces[face].first());
//                                     }
//                                 }
//                                 nextToCheck = auxToCheck;
//                             }

//                             contactFaces.append(newContactFaces);
//                         }
//                     }
//                 }
//             }
//         }
//     }

//     bool inContact = false;

//     forAll(contactFaces,list)
//     {
//         bool firstFace = true;
//         point bPoint = vector::zero;
//         forAll(contactFaces[list],faceI)
//         {
//             face currentFace = mesh.faces()[contactFaces[list][faceI].first()];
//             pointHit pointH = currentFace.nearestPoint(cCenter,mesh.points());
//             if(mag(pointH.rawPoint() - cCenter) < mag(bPoint - cCenter) || firstFace)
//             {
//                 bPoint = pointH.rawPoint();
//                 firstFace = false;
//             }
//         }

//         vector nVec = (cCenter - bPoint)/mag(cCenter - bPoint);

//         wallContactVars& wallCntVars(wallCntInfo.getWallCntVars(
//             bPoint,
//             nVec,
//             deltaT
//         ));

//         wallCntVars.contactCenter_ = bPoint;
//         wallCntVars.contactNormal_ = nVec;
//         wallCntVars.contactArea_ = sphereContactArea
//         (
//             mesh,
//             wallCntInfo.getcClass(),
//             wallCntInfo.getcVars(),
//             nVec,
//             bPoint
//         );
//         wallCntVars.contactVolume_ = getInterVolume_Sphere
//         (
//             mesh,
//             wallCntInfo.getcClass(),
//             wallCntInfo.getcVars(),
//             nVec,
//             bPoint
//         );

//         wallCntVars.setMeanCntPars
//         (
//             mesh,
//             contactFaces[list],
//             wallCntInfo.getWallMeanPars()
//         );

//         if(wallCntVars.contactVolume_ > 0)
//         {
//             inContact = true;
//         }
//     }

//     reduce(inContact, orOp<bool>());

//     if(inContact)
//     {
//         wallCntInfo.getcClass().setWallContact(true);
//         wallCntInfo.getcClass().inContactWithStatic(true);
//     }
// }
// //---------------------------------------------------------------------------//
// DynamicList<Tuple2<label,string>> getContactFacesSphere
// (
//     const fvMesh&   mesh,
//     wallContactInfo& wallCntInfo,
//     label faceLabel,
//     DynamicList<Tuple2<label,string>>& checkedFaces
// )
// {
//     scalar cRadius(wallCntInfo.getcClass().getGeomModel().getDC() / 2);
//     vector cCenter(wallCntInfo.getcClass().getGeomModel().getCoM());

//     DynamicList<Tuple2<label,string>> facesReturnList;

//     labelList faceEdges = mesh.faceEdges()[faceLabel];
//     forAll(faceEdges,edge)
//     {
//         labelList edgeFaces = mesh.edgeFaces()[faceEdges[edge]];
//         forAll(edgeFaces,faceI)
//         {
//             if (!mesh.isInternalFace(edgeFaces[faceI]))
//             {
//                 // get reference to the patch which is in contact with IB.
//                 // There is contact only if the patch is marked as a wall
//                 label facePatchId = mesh.boundaryMesh().whichPatch(edgeFaces[faceI]);
//                 const polyPatch& cPatch = mesh.boundaryMesh()[facePatchId];
//                 if (wallCntInfo.isContactPatch(cPatch.name()))
//                 {
//                     vector nVec(-mesh.Sf()[edgeFaces[faceI]]/mag(mesh.Sf()[edgeFaces[faceI]]));

//                     if(!case3D)
//                     {
//                         if(mag(emptyDir - nVec) < 1e-4 || mag(emptyDir - (-1)*nVec) < 1e-4)
//                         {
//                             continue;
//                         }
//                     }

//                     bool cont = false;

//                     forAll(checkedFaces,face)
//                     {
//                         if(checkedFaces[face].first() == edgeFaces[faceI])
//                         {
//                             cont = true;
//                             break;
//                         }
//                     }

//                     if(cont)
//                     {
//                         continue;
//                     }

//                     face currentFace = mesh.faces()[edgeFaces[faceI]];
//                     pointHit pointH = currentFace.nearestPoint(cCenter,mesh.points());
//                     if(mag(pointH.rawPoint() - cCenter) < cRadius)
//                     {
//                         facesReturnList.append(Tuple2<label,string>(
//                             edgeFaces[faceI],
//                             cPatch.name()
//                         ));
//                     }
//                 }
//             }
//         }
//     }

//     return facesReturnList;
// }
// //---------------------------------------------------------------------------//
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

        getWallContactVars(
            mesh,
            cWallCntI,
            deltaT,
            sWC
        );
    }
}
// //---------------------------------------------------------------------------//
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
    InfoH << parallelDEM_Info << "-- Particle-wall id "<< sCI.getBodyId() <<" contact center "
        << wallCntVar.contactCenter_ << endl;
    InfoH << parallelDEM_Info << "-- Particle-wall id "<< sCI.getBodyId() <<" contact normal "
        << wallCntVar.contactNormal_ << endl;
    InfoH << parallelDEM_Info << "-- Particle-wall id "<< sCI.getBodyId() <<" contact volume "
        << wallCntVar.contactVolume_ << endl;
    InfoH << parallelDEM_Info << "-- Particle-wall id "<< sCI.getBodyId() <<" contact area "
        << wallCntVar.contactArea_ << endl;

    vector F = sCI.getFNe(wallCntVar);
    InfoH << parallelDEM_Info << "-- Particle-wall id "<< sCI.getBodyId() <<" contact FNe " << F << endl;

    vector FNd = sCI.getFNd(wallCntVar);
    InfoH << parallelDEM_Info << "-- Particle-wall id "<< sCI.getBodyId() <<" contact FNd " << FNd << endl;

    if ((F & FNd) < 0 && mag(FNd) > mag(F))
    {
        FNd *= mag(F) / mag(FNd);
        InfoH << parallelDEM_Info << "FNd was Clipped to "<< FNd << endl;
    }

    F += FNd;
    InfoH << parallelDEM_Info << "-- Particle-wall id "<< sCI.getBodyId() <<" contact FN " << F << endl;

    vector Ft = sCI.getFt(wallCntVar, deltaT);
    InfoH << parallelDEM_Info << "-- Particle-wall id "<< sCI.getBodyId() <<" contact Ft " << Ft << endl;

    if (mag(Ft) > sCI.getMu(wallCntVar) * mag(F))
    {
        Ft *= sCI.getMu(wallCntVar) * mag(F) / mag(Ft);
    }
    InfoH << parallelDEM_Info << "-- Particle-wall id "<< sCI.getBodyId() <<" contact Ft clamped" << Ft << endl;
    F += Ft;

    vector FA = sCI.getFA(wallCntVar);
    InfoH << parallelDEM_Info << "-- Particle-wall id "<< sCI.getBodyId() <<" contact FA " << FA << endl;
    F -= FA;

    outF += F;

    cLVecOut += wallCntVar.lVec_;

    sCI.getOutForce().F = outF;
    sCI.getOutForce().T = cLVecOut ^  outF;
    return true;
}
// //---------------------------------------------------------------------------//
// scalar getInterVolume_ArbShape(
//     const fvMesh&   mesh,
//     ibContactClass& cClass,
//     ibContactVars& cVars,
//     vector nVec,
//     vector center
// )
// {
//     return (cVars.M0_-cVars.M_)/(cVars.rhoS_.value() + SMALL);
// }
// //---------------------------------------------------------------------------//
// scalar getInterVolume_Sphere(
//     const fvMesh&   mesh,
//     ibContactClass& cClass,
//     ibContactVars& cVars,
//     vector nVec,
//     vector center
// )
// {
//     scalar cRadius(cClass.getGeomModel().getDC() / 2);
//     vector cCenter(cClass.getGeomModel().getCoM());

//     plane contPlane(center, nVec);
//     scalar dist = contPlane.distance(cCenter);
//     scalar xH = cRadius - dist;

//     if(case3D)
//     {
//         return (Foam::constant::mathematical::pi*sqr(xH)*(3*cRadius - xH)) / 3;
//     }
//     else
//     {
//         boundBox meshBounds = mesh.bounds();
//         scalar emptyLength = meshBounds.max()[emptyDim]
//                             - meshBounds.min()[emptyDim];

//         scalar cCirSeg = sqr(cRadius)*acos((dist)/cRadius)
//                          - (dist)*sqrt(sqr(cRadius) - sqr(dist));
//         return cCirSeg*emptyLength;
//     }
// }
// //---------------------------------------------------------------------------//
// scalar sphereContactArea
// (
//     const fvMesh&   mesh,
//     ibContactClass& cClass,
//     ibContactVars& cVars,
//     vector nVec,
//     vector center
// )
// {
//     scalar cRadius(cClass.getGeomModel().getDC() / 2);
//     vector cCenter(cClass.getGeomModel().getCoM());

//     plane contPlane(center, nVec);

//     scalar dist = contPlane.distance(cCenter);
//     scalar contactRad = sqrt(sqr(cRadius) - sqr(dist));

//     if(case3D)
//     {
//         return Foam::constant::mathematical::pi*sqr(contactRad);
//     }
//     else
//     {
//         boundBox meshBounds = mesh.bounds();
//         scalar emptyLength = meshBounds.max()[emptyDim]
//                             - meshBounds.min()[emptyDim];

//         return contactRad*emptyLength;
//     }
// }
// //---------------------------------------------------------------------------//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace contactModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
