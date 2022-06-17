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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace contactModel
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//---------------------------------------------------------------------------//
void detectWallContact(
    const fvMesh&   mesh,
    wallContactInfo& wallCntInfo
)
{
    wallCntInfo.getcClass().setWallContact(false);

    if(wallCntInfo.getcClass().getGeomModel().getcType() == sphere)
    {
        detectWallContact<sphere>
        (
            mesh,
            wallCntInfo
        );
    }
    else if(wallCntInfo.getcClass().getGeomModel().getcType() == cluster)
    {
        detectWallContact<cluster>
        (
            mesh,
            wallCntInfo
        );
    }
    else
    {
        detectWallContact<arbShape>(
            mesh,
            wallCntInfo
        );
    }
}
//---------------------------------------------------------------------------//
template <contactType cT>
void detectWallContact(
    const fvMesh&   mesh,
    wallContactInfo& wallCntInfo
)
{
    bool inContact = false;

    label nCells = mesh.nCells();

    List<DynamicLabelList> contactFaces;
    List<DynamicLabelList>& surfCells(wallCntInfo.getcClass().getSurfCells());

    // go through all surfCells and check if there is any surfCell whose face is a boundary face
    forAll (surfCells[Pstream::myProcNo()],sCellI)
    {
        label cCell(surfCells[Pstream::myProcNo()][sCellI]);
        if(cCell < nCells)
        {
            const labelList& cFaces = mesh.cells()[cCell];
            forAll (cFaces,faceI)
            {
                if (!mesh.isInternalFace(cFaces[faceI]))
                {
                    // get reference to the patch which is in contact with IB. There is contact only if the patch is marked as a wall
                    label facePatchId = mesh.boundaryMesh().whichPatch(cFaces[faceI]);
                    const polyPatch& cPatch = mesh.boundaryMesh()[facePatchId];
                    if (cPatch.type()=="wall")
                    {
                        vector nVec(-mesh.Sf()[cFaces[faceI]]/mag(mesh.Sf()[cFaces[faceI]]));

                        if(!case3D)
                        {
                            if(mag(emptyDir - nVec) < 1e-4 || mag(emptyDir - (-1)*nVec) < 1e-4)
                            {
                                continue;
                            }
                        }

                        labelList facePoints(mesh.faces()[cFaces[faceI]]);//list of vertex indicies

                        forAll(facePoints,pointI)
                        {
                            if(wallCntInfo.getcClass().getGeomModel().pointInside(mesh.points()[facePoints[pointI]]))
                            {
                                inContact = true;
                                break;
                            }
                        }
                    }
                }
                if(inContact)
                    break;
            }
        }
        if(inContact)
            break;
    }

    reduce(inContact, orOp<bool>());

    if(inContact)
    {
        wallCntInfo.getcClass().setWallContact(true);
        wallCntInfo.getcClass().inContactWithStatic(true);
    }
}
//---------------------------------------------------------------------------//
template <>
void detectWallContact <sphere>(
    const fvMesh&   mesh,
    wallContactInfo& wallCntInfo
)
{
    scalar cRadius(wallCntInfo.getcClass().getGeomModel().getDC() / 2);
    vector cCenter(wallCntInfo.getcClass().getGeomModel().getCoM());

    label nCells = mesh.nCells();

    List<DynamicLabelList>& surfCells(wallCntInfo.getcClass().getSurfCells());
    bool inContact = false;

    // go through all surfCells and check if there is any surfCell whose face is a boundary face
    forAll (surfCells[Pstream::myProcNo()],sCellI)
    {
        label cCell(surfCells[Pstream::myProcNo()][sCellI]);
        if(cCell < nCells)
        {
            const labelList& cFaces = mesh.cells()[cCell];
            forAll (cFaces,faceI)
            {
                if (!mesh.isInternalFace(cFaces[faceI]))
                {
                    // get reference to the patch which is in contact with IB. There is contact only if the patch is marked as a wall
                    label facePatchId = mesh.boundaryMesh().whichPatch(cFaces[faceI]);
                    const polyPatch& cPatch = mesh.boundaryMesh()[facePatchId];
                    if (cPatch.type()=="wall")
                    {
                        vector nVec(-mesh.Sf()[cFaces[faceI]]/mag(mesh.Sf()[cFaces[faceI]]));

                        if(!case3D)
                        {
                            if(mag(emptyDir - nVec) < 1e-4 || mag(emptyDir - (-1)*nVec) < 1e-4)
                            {
                                continue;
                            }
                        }

                        face currentFace = mesh.faces()[cFaces[faceI]];
                        pointHit pointH = currentFace.nearestPoint(cCenter,mesh.points());
                        if(mag(pointH.rawPoint() - cCenter) < cRadius)
                        {
                            inContact = true;
                        }
                    }
                }
                if(inContact)
                    break;
            }
        }
        if(inContact)
            break;
    }

    reduce(inContact, orOp<bool>());

    if(inContact)
    {
        wallCntInfo.getcClass().setWallContact(true);
        wallCntInfo.getcClass().inContactWithStatic(true);
    }
}
//---------------------------------------------------------------------------//
template <>
void detectWallContact <cluster>(
    const fvMesh&   mesh,
    wallContactInfo& wallCntInfo
)
{
    clusterBody& cCluster = dynamic_cast<clusterBody&>(wallCntInfo.getcClass().getGeomModel());
    PtrList<geomModel>& cBodies = cCluster.getClusterBodies();

    forAll(cBodies, cIbI)
    {
        autoPtr<geomModel> cGeomModel(cBodies[cIbI].getGeomModel());
        autoPtr<ibContactClass> cIbClassI(new ibContactClass(
            cGeomModel,
            wallCntInfo.getcClass().getMatInfo()
        ));
        autoPtr<wallContactInfo> cWallCntI(new wallContactInfo(
            cIbClassI(),
            wallCntInfo.getcVars(),
            wallCntInfo.getWInfos(),
            wallCntInfo.getMatInterAdh()
        ));

        detectWallContact(
            mesh,
            cWallCntI()
        );

        if(cWallCntI().getcClass().checkWallContact())
        {
            wallCntInfo.getcClass().setWallContact(true);
        }

        if(cWallCntI().getcClass().checkInContactWithStatic())
        {
            wallCntInfo.getcClass().inContactWithStatic(true);
        }
    }
}
//---------------------------------------------------------------------------//
void getWallContactVars(
    const fvMesh&   mesh,
    wallContactInfo& wallCntInfo,
    const scalar deltaT
)
{
    if(wallCntInfo.getcClass().getGeomModel().getcType() == sphere)
    {
        getWallContactVars<sphere>
        (
            mesh,
            wallCntInfo,
            deltaT
        );
        wallCntInfo.clearWallCntVars();
    }
    else if(wallCntInfo.getcClass().getGeomModel().getcType() == cluster)
    {
        getWallContactVars<cluster>
        (
            mesh,
            wallCntInfo,
            deltaT
        );
        wallCntInfo.clearWallCntVars();
    }
    else
    {
        getWallContactVars<arbShape>(
            mesh,
            wallCntInfo,
            deltaT
        );
        wallCntInfo.clearWallCntVars();
    }
}
//---------------------------------------------------------------------------//
template <contactType cT>
void getWallContactVars(
    const fvMesh&   mesh,
    wallContactInfo& wallCntInfo,
    const scalar deltaT
)
{
    bool inContact(false);

    DynamicVectorList center;
    DynamicVectorList normal;
    DynamicLabelList centerPoints;
    DynamicScalarList area;
    HashTable<bool,label,Hash<label>> pointTable;

    label nCells = mesh.nCells();

    // Tuple: first() = face label; second() = patch name
    List<DynamicList<Tuple2<label,string>>> contactFaces;
    List<DynamicLabelList>& surfCells(wallCntInfo.getcClass().getSurfCells());

    // go through all surfCells and check if there is any surfCell whose face is a boundary face
    forAll (surfCells[Pstream::myProcNo()],sCellI)
    {
        label cCell(surfCells[Pstream::myProcNo()][sCellI]);
        if(cCell < nCells)
        {
            const labelList& cFaces = mesh.cells()[cCell];
            forAll (cFaces,faceI)
            {
                if (!mesh.isInternalFace(cFaces[faceI]))
                {
                    // get reference to the patch which is in contact with IB. There is contact only if the patch is marked as a wall
                    label facePatchId = mesh.boundaryMesh().whichPatch(cFaces[faceI]);
                    const polyPatch& cPatch = mesh.boundaryMesh()[facePatchId];
                    if (wallCntInfo.isContactPatch(cPatch.name()))
                    {
                        vector nVec(-mesh.Sf()[cFaces[faceI]]/mag(mesh.Sf()[cFaces[faceI]]));

                        if(!case3D)
                        {
                            if(mag(emptyDir - nVec) < 1e-4 || mag(emptyDir - (-1)*nVec) < 1e-4)
                            {
                                continue;
                            }
                        }

                        bool cont = false;
                        forAll(contactFaces,list)
                        {
                            forAll(contactFaces[list],face)
                            {
                                if(contactFaces[list][face].first() == cFaces[faceI])
                                {
                                    cont = true;
                                    break;
                                }
                            }
                            if(cont)
                                break;
                        }
                        if(cont)
                            continue;

                        labelList facePoints(mesh.faces()[cFaces[faceI]]);//list of vertex indicies
                        label numOfPoints = 0;

                        forAll(facePoints,pointI)
                        {
                            if(pointTable.found(facePoints[pointI]))
                            {
                                if(pointTable[facePoints[pointI]])
                                {
                                    numOfPoints++;
                                    inContact = true;
                                }
                            }
                            else
                            {
                                if(wallCntInfo.getcClass().getGeomModel().pointInside(mesh.points()[facePoints[pointI]]))
                                {
                                    numOfPoints++;
                                    inContact = true;
                                    pointTable.insert(facePoints[pointI], true);
                                }
                                else
                                {
                                    pointTable.insert(facePoints[pointI], false);
                                }
                            }
                        }

                        if(numOfPoints > 0)
                        {
                            center.append(vector::zero);
                            normal.append(vector::zero);
                            centerPoints.append(0);
                            area.append(0);

                            DynamicList<Tuple2<label,string>> newContactFaces
                            (
                                1,
                                Tuple2<label,string>(cFaces[faceI], cPatch.name())
                            );

                            labelList nextToCheck(1,cFaces[faceI]);

                            while (nextToCheck.size() > 0)
                            {
                                DynamicLabelList auxToCheck;

                                forAll (nextToCheck,faceToCheck)
                                {
                                    DynamicList<Tuple2<label,string>> newfaces =
                                    getContactFacesArbShape
                                    (
                                        mesh,
                                        wallCntInfo,
                                        nextToCheck[faceToCheck],
                                        newContactFaces,
                                        pointTable,
                                        center,
                                        normal,
                                        centerPoints,
                                        area
                                    );

                                    forAll(newfaces,face)
                                    {
                                        newContactFaces.append(newfaces[face]);
                                        auxToCheck.append(newfaces[face].first());
                                    }
                                }
                                nextToCheck = auxToCheck;
                            }

                            contactFaces.append(newContactFaces);
                        }
                    }
                }
            }
        }
    }

    reduce(inContact, orOp<bool>());

    if(inContact)
    {
        wallCntInfo.getcClass().setWallContact(true);
        wallCntInfo.getcClass().inContactWithStatic(true);

        forAll(center, cntI)
        {
            center[cntI] /= centerPoints[cntI];
            normal[cntI] /= centerPoints[cntI];
            normal[cntI] /= mag(normal[cntI]);

            wallContactVars& wallCntVars(wallCntInfo.getWallCntVars(
                center[cntI],
                normal[cntI],
                deltaT
            ));

            wallCntVars.contactArea_ = area[cntI];
            wallCntVars.contactCenter_ = center[cntI];
            wallCntVars.contactVolume_ = getInterVolume<cT>
            (
                mesh,
                wallCntInfo.getcClass(),
                wallCntInfo.getcVars(),
                normal[cntI],
                center[cntI]
            );
            wallCntVars.contactNormal_ = normal[cntI];

            wallCntVars.setMeanCntPars
            (
                mesh,
                contactFaces[cntI],
                wallCntInfo.getWallMeanPars()
            );
        }
    }
}
//---------------------------------------------------------------------------//
DynamicList<Tuple2<label,string>> getContactFacesArbShape
(
    const fvMesh&   mesh,
    wallContactInfo& wallCntInfo,
    label faceLabel,
    DynamicList<Tuple2<label,string>>& checkedFaces,
    HashTable<bool,label,Hash<label>>& pointTable,
    DynamicVectorList& center,
    DynamicVectorList& normal,
    DynamicLabelList& centerPoints,
    DynamicScalarList& area
)
{
    DynamicList<Tuple2<label,string>> facesReturnList;

    labelList faceEdges = mesh.faceEdges()[faceLabel];
    forAll(faceEdges,edge)
    {
        labelList edgeFaces = mesh.edgeFaces()[faceEdges[edge]];
        forAll(edgeFaces,faceI)
        {
            if (!mesh.isInternalFace(edgeFaces[faceI]))
            {
                // get reference to the patch which is in contact with IB.
                // There is contact only if the patch is marked as a wall
                label facePatchId = mesh.boundaryMesh().whichPatch(edgeFaces[faceI]);
                const polyPatch& cPatch = mesh.boundaryMesh()[facePatchId];
                if (wallCntInfo.isContactPatch(cPatch.name()))
                {
                    vector nVec(-mesh.Sf()[edgeFaces[faceI]]/mag(mesh.Sf()[edgeFaces[faceI]]));

                    if(!case3D)
                    {
                        if(mag(emptyDir - nVec) < 1e-4 || mag(emptyDir - (-1)*nVec) < 1e-4)
                            continue;
                    }

                    bool cont = false;

                    forAll(checkedFaces,face)
                    {
                        if(checkedFaces[face].first() == edgeFaces[faceI])
                        {
                            cont = true;
                            break;
                        }
                    }

                    if(cont)
                        continue;

                    labelList facePoints(mesh.faces()[edgeFaces[faceI]]);//list of vertex indicies
                    label numOfPoints = 0;

                    forAll(facePoints,pointI)
                    {
                        if(pointTable.found(facePoints[pointI]))
                        {
                            if(pointTable[facePoints[pointI]])
                            {
                                numOfPoints++;
                                center.last() += mesh.points()[facePoints[pointI]];
                                normal.last() += nVec;
                            }
                        }
                        else
                        {
                            if(wallCntInfo.getcClass().getGeomModel().pointInside(mesh.points()[facePoints[pointI]]))
                            {
                                numOfPoints++;
                                center.last() += mesh.points()[facePoints[pointI]];
                                normal.last() += nVec;
                                pointTable.insert(facePoints[pointI], true);
                            }
                            else
                            {
                                pointTable.insert(facePoints[pointI], false);
                            }
                        }
                    }

                    if(numOfPoints > 0)
                    {
                        centerPoints.last() += numOfPoints;
                        area.last() += mag(mesh.Sf()[edgeFaces[faceI]])
                            *numOfPoints/facePoints.size();
                        facesReturnList.append(Tuple2<label,string>
                        (
                            edgeFaces[faceI],
                            cPatch.name()
                        ));
                    }
                }
            }
        }
    }

    return facesReturnList;
}
//---------------------------------------------------------------------------//
template <>
void getWallContactVars <sphere>(
    const fvMesh&   mesh,
    wallContactInfo& wallCntInfo,
    const scalar deltaT
)
{
    scalar cRadius(wallCntInfo.getcClass().getGeomModel().getDC() / 2);
    vector cCenter(wallCntInfo.getcClass().getGeomModel().getCoM());

    label nCells = mesh.nCells();

    // Tuple: first() = face label; second() = patch name
    List<DynamicList<Tuple2<label,string>>> contactFaces;

    DynamicVectorList center(Pstream::nProcs(),vector::zero);
    DynamicVectorList normal(Pstream::nProcs(),vector::zero);

    List<DynamicLabelList>& surfCells(wallCntInfo.getcClass().getSurfCells());

    // go through all surfCells and check if there is any surfCell whose face is a boundary face
    forAll (surfCells[Pstream::myProcNo()],sCellI)
    {
        label cCell(surfCells[Pstream::myProcNo()][sCellI]);
        if(cCell < nCells)
        {
            const labelList& cFaces = mesh.cells()[cCell];
            forAll (cFaces,faceI)
            {
                if (!mesh.isInternalFace(cFaces[faceI]))
                {
                    // get reference to the patch which is in contact with IB. There is contact only if the patch is marked as a wall
                    label facePatchId = mesh.boundaryMesh().whichPatch(cFaces[faceI]);
                    const polyPatch& cPatch = mesh.boundaryMesh()[facePatchId];
                    if (wallCntInfo.isContactPatch(cPatch.name()))
                    {
                        vector nVec(-mesh.Sf()[cFaces[faceI]]/mag(mesh.Sf()[cFaces[faceI]]));

                        if(!case3D)
                        {
                            if(mag(emptyDir - nVec) < 1e-4 || mag(emptyDir - (-1)*nVec) < 1e-4)
                            {
                                continue;
                            }
                        }

                        bool cont = false;
                        forAll(contactFaces,list)
                        {
                            forAll(contactFaces[list],face)
                            {
                                if(contactFaces[list][face].first() == cFaces[faceI])
                                {
                                    cont = true;
                                    break;
                                }
                            }
                            if(cont)
                                break;
                        }
                        if(cont)
                            continue;

                        face currentFace = mesh.faces()[cFaces[faceI]];
                        pointHit pointH = currentFace.nearestPoint(cCenter,mesh.points());
                        if(mag(pointH.rawPoint() - cCenter) < cRadius)
                        {
                            DynamicList<Tuple2<label,string>> newContactFaces(
                                1,
                                Tuple2<label,string>(cFaces[faceI], cPatch.name())
                            );

                            labelList nextToCheck(1,cFaces[faceI]);

                            while (nextToCheck.size() > 0)
                            {
                                DynamicLabelList auxToCheck;

                                forAll (nextToCheck,faceToCheck)
                                {
                                    DynamicList<Tuple2<label,string>> newfaces =
                                    getContactFacesSphere(
                                                mesh,
                                                wallCntInfo,
                                                nextToCheck[faceToCheck],
                                                newContactFaces
                                    );

                                    forAll(newfaces,face)
                                    {
                                        newContactFaces.append(newfaces[face]);
                                        auxToCheck.append(newfaces[face].first());
                                    }
                                }
                                nextToCheck = auxToCheck;
                            }

                            contactFaces.append(newContactFaces);
                        }
                    }
                }
            }
        }
    }

    bool inContact = false;

    forAll(contactFaces,list)
    {
        bool firstFace = true;
        point bPoint = vector::zero;
        forAll(contactFaces[list],faceI)
        {
            face currentFace = mesh.faces()[contactFaces[list][faceI].first()];
            pointHit pointH = currentFace.nearestPoint(cCenter,mesh.points());
            if(mag(pointH.rawPoint() - cCenter) < mag(bPoint - cCenter) || firstFace)
            {
                bPoint = pointH.rawPoint();
                firstFace = false;
            }
        }

        vector nVec = (cCenter - bPoint)/mag(cCenter - bPoint);

        wallContactVars& wallCntVars(wallCntInfo.getWallCntVars(
            bPoint,
            nVec,
            deltaT
        ));

        wallCntVars.contactCenter_ = bPoint;
        wallCntVars.contactNormal_ = nVec;
        wallCntVars.contactArea_ = sphereContactArea
        (
            mesh,
            wallCntInfo.getcClass(),
            wallCntInfo.getcVars(),
            nVec,
            bPoint
        );
        wallCntVars.contactVolume_ = getInterVolume<sphere>
        (
            mesh,
            wallCntInfo.getcClass(),
            wallCntInfo.getcVars(),
            nVec,
            bPoint
        );

        wallCntVars.setMeanCntPars
        (
            mesh,
            contactFaces[list],
            wallCntInfo.getWallMeanPars()
        );

        if(wallCntVars.contactVolume_ > 0)
        {
            inContact = true;
        }
    }

    reduce(inContact, orOp<bool>());

    if(inContact)
    {
        wallCntInfo.getcClass().setWallContact(true);
        wallCntInfo.getcClass().inContactWithStatic(true);
    }
}
//---------------------------------------------------------------------------//
DynamicList<Tuple2<label,string>> getContactFacesSphere
(
    const fvMesh&   mesh,
    wallContactInfo& wallCntInfo,
    label faceLabel,
    DynamicList<Tuple2<label,string>>& checkedFaces
)
{
    scalar cRadius(wallCntInfo.getcClass().getGeomModel().getDC() / 2);
    vector cCenter(wallCntInfo.getcClass().getGeomModel().getCoM());

    DynamicList<Tuple2<label,string>> facesReturnList;

    labelList faceEdges = mesh.faceEdges()[faceLabel];
    forAll(faceEdges,edge)
    {
        labelList edgeFaces = mesh.edgeFaces()[faceEdges[edge]];
        forAll(edgeFaces,faceI)
        {
            if (!mesh.isInternalFace(edgeFaces[faceI]))
            {
                // get reference to the patch which is in contact with IB.
                // There is contact only if the patch is marked as a wall
                label facePatchId = mesh.boundaryMesh().whichPatch(edgeFaces[faceI]);
                const polyPatch& cPatch = mesh.boundaryMesh()[facePatchId];
                if (wallCntInfo.isContactPatch(cPatch.name()))
                {
                    vector nVec(-mesh.Sf()[edgeFaces[faceI]]/mag(mesh.Sf()[edgeFaces[faceI]]));

                    if(!case3D)
                    {
                        if(mag(emptyDir - nVec) < 1e-4 || mag(emptyDir - (-1)*nVec) < 1e-4)
                        {
                            continue;
                        }
                    }

                    bool cont = false;

                    forAll(checkedFaces,face)
                    {
                        if(checkedFaces[face].first() == edgeFaces[faceI])
                        {
                            cont = true;
                            break;
                        }
                    }

                    if(cont)
                    {
                        continue;
                    }

                    face currentFace = mesh.faces()[edgeFaces[faceI]];
                    pointHit pointH = currentFace.nearestPoint(cCenter,mesh.points());
                    if(mag(pointH.rawPoint() - cCenter) < cRadius)
                    {
                        facesReturnList.append(Tuple2<label,string>(
                            edgeFaces[faceI],
                            cPatch.name()
                        ));
                    }
                }
            }
        }
    }

    return facesReturnList;
}
//---------------------------------------------------------------------------//
template <>
void getWallContactVars <cluster>(
    const fvMesh&   mesh,
    wallContactInfo& wallCntInfo,
    const scalar deltaT
)
{
    clusterBody& cCluster = dynamic_cast<clusterBody&>(wallCntInfo.getcClass().getGeomModel());
    PtrList<geomModel>& cBodies = cCluster.getClusterBodies();

    forAll(cBodies, cIbI)
    {
        autoPtr<geomModel> cGeomModel(cBodies[cIbI].getGeomModel());
        autoPtr<ibContactClass> cIbClassI(new ibContactClass(
            cGeomModel,
            wallCntInfo.getcClass().getMatInfo()
        ));
        autoPtr<wallContactInfo> cWallCntI(new wallContactInfo(
            cIbClassI(),
            wallCntInfo.getcVars(),
            wallCntInfo.getWInfos(),
            wallCntInfo.getMatInterAdh()
        ));

        getWallContactVars(
            mesh,
            cWallCntI(),
            deltaT
        );

        DynamicLabelList& wallCntHashesI(cWallCntI().getWallCntHashes());
        forAll(wallCntHashesI,contVar)
        {
            wallContactVars& wallCntVarI(
                cWallCntI().getWallcVarsTable()[wallCntHashesI[contVar]]
            );

            label newHash = wallCntInfo.getWallCntHashes().size();
            wallCntInfo.getWallcVarsTable().insert(newHash, wallCntVarI);
            wallCntInfo.getWallCntHashes().append(newHash);
            wallCntInfo.getCurUsedHashes().append(wallCntInfo.getWallCntHashes().last());
        }
    }
}
//---------------------------------------------------------------------------//
void solveWallContact
(
    const fvMesh&   mesh,
    wallContactInfo& wallCntInfo,
    scalar deltaT
)
{
    getWallContactVars(
        mesh,
        wallCntInfo,
        deltaT
    );

    vector outF = vector::zero;
    vector cLVecOut = vector::zero;
    DynamicLabelList& wallCntHashes(wallCntInfo.getWallCntHashes());
    label numOfCntVars = wallCntHashes.size();
    forAll(wallCntHashes,contVar)
    {
        wallContactVars& wallCntVar(
            wallCntInfo.getWallcVarsTable()[wallCntHashes[contVar]]
        );

        wallCntInfo.evalVariables(wallCntVar);

        if(wallCntVar.contactVolume_ == 0)
        {
            continue;
        }

        InfoH << DEM_Info << "-- Detected Particle-wall contact: -- body "
            << wallCntInfo.getBodyId() << endl;
        InfoH << DEM_Info << "-- Particle-wall contact center "
            << wallCntVar.contactCenter_ << endl;
        InfoH << DEM_Info << "-- Particle-wall contact normal "
            << wallCntVar.contactNormal_ << endl;
        InfoH << DEM_Info << "-- Particle-wall contact volume "
            << wallCntVar.contactVolume_ << endl;
        InfoH << DEM_Info << "-- Particle-wall contact area "
            << wallCntVar.contactArea_ << endl;

        vector F = wallCntInfo.getFNe(wallCntVar);
        InfoH << DEM_Info << "-- Particle-wall contact FNe " << F << endl;

        vector FNd = wallCntInfo.getFNd(wallCntVar);
        InfoH << DEM_Info << "-- Particle-wall contact FNd " << FNd << endl;

        F += FNd;
        InfoH << DEM_Info << "-- Particle-wall contact FN " << F << endl;

        vector Ft = wallCntInfo.getFt(wallCntVar, deltaT);
        InfoH << DEM_Info << "-- Particle-wall contact Ft " << Ft << endl;

        if (mag(Ft) > wallCntInfo.getMu(wallCntVar) * mag(F))
        {
            Ft *= wallCntInfo.getMu(wallCntVar) * mag(F) / mag(Ft);
        }
        InfoH << DEM_Info << "-- Particle-wall contact Ft clamped" << Ft << endl;
        F += Ft;        

        vector FA = wallCntInfo.getFA(wallCntVar);
        InfoH << DEM_Info << "-- Particle-wall contact FA " << FA << endl;
        F -= FA;

        outF += F + Ft;
        cLVecOut += wallCntVar.lVec_;
    }

    reduce(numOfCntVars, sumOp<label>());
    if(numOfCntVars == 0)
    {
        wallCntInfo.getOutForce().F = vector::zero;
        wallCntInfo.getOutForce().T = vector::zero;
        return;
    }

    reduce(outF, sumOp<vector>());
    reduce(cLVecOut, sumOp<vector>());
    cLVecOut /= numOfCntVars;
    if(mag(outF) == 0)
    {
        wallCntInfo.getOutForce().F = vector::zero;
        wallCntInfo.getOutForce().T = vector::zero;
        return;
    }

    InfoH << DEM_Info << "-- Resolved Particle-wall contact: -- body "
        << wallCntInfo.getBodyId() << endl;

    wallCntInfo.getOutForce().F = outF;
    wallCntInfo.getOutForce().T = cLVecOut ^  outF;
}
//---------------------------------------------------------------------------//
template <contactType cT>
scalar
getInterVolume(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactVars& cVars,
    vector nVec,
    vector center
)
{
    return (cVars.M0_-cVars.M_)/(cVars.rhoS_.value() + SMALL);
}
//---------------------------------------------------------------------------//
template <>
scalar
getInterVolume <sphere>(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactVars& cVars,
    vector nVec,
    vector center
)
{
    scalar cRadius(cClass.getGeomModel().getDC() / 2);
    vector cCenter(cClass.getGeomModel().getCoM());

    plane contPlane(center, nVec);
    scalar dist = contPlane.distance(cCenter);
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
    vector nVec,
    vector center
)
{
    scalar cRadius(cClass.getGeomModel().getDC() / 2);
    vector cCenter(cClass.getGeomModel().getCoM());

    plane contPlane(center, nVec);

    scalar dist = contPlane.distance(cCenter);
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
