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
    Martin Isoz (2019-*), Martin Šourek (2019-*),
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
    const dynamicFvMesh&   mesh,
    contactInfo& cInfo,
    contactVars* cVars,
    vector geometricD
)
{
    const contactType cT(cInfo.getGeomModel().getcType());

    if(cT == sphere)
    {
        return detectWallContact<sphere>
        (
            mesh,
            cInfo,
            cVars,
            geometricD
        );
    }
    else
    {
        return detectWallContact<arbShape>(
            mesh,
            cInfo,
            cVars,
            geometricD
        );
    }
}
//---------------------------------------------------------------------------//
template <contactType cT>
void detectWallContact(
    const dynamicFvMesh&   mesh,
    contactInfo& cInfo,
    contactVars* cVars,
    vector geometricD
)
{
    bool inContact(false);
    cInfo.clearContactInfo();

    vector center = vector::zero;
    vector normal = vector::zero;
    label centerPoints = 0;

    scalar area = 0;

    bool case3D = true;
    vector emptydir = vector::zero;

    forAll (geometricD, direction)
    {
        if (geometricD[direction] == -1)
        {
            case3D = false;
            emptydir[direction] = 1;
            break;
        }
    }

    // go through all surfCells and check if there is any surfCell whose face is a boundary face
    forAll (cInfo.getSurfCells()[Pstream::myProcNo()],sCellI)
    {
        label cCell(cInfo.getSurfCells()[Pstream::myProcNo()][sCellI]);

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
                        if(mag(emptydir - nVec) < 1e-4)
                        {
                            continue;
                        }
                    }

                    labelList facePoints(mesh.faces()[cFaces[faceI]]);//list of vertex indicies
                    label numOfPoints = 0;

                    forAll(facePoints,pointI)
                    {
                        if(cInfo.getGeomModel().pointInside(mesh.points()[facePoints[pointI]]))
                        {
                            numOfPoints++;
                            center += mesh.points()[facePoints[pointI]];
                            normal += nVec;
                            inContact = true;
                        }
                    }
                    centerPoints += numOfPoints;
                    area += mag(mesh.Sf()[cFaces[faceI]]) * numOfPoints / facePoints.size();
                }
            }
        }
    }

    reduce(inContact, orOp<bool>());

    if(inContact)
    {
        cInfo.setWallContact(true);
        forAll(cInfo.getIntCells()[Pstream::myProcNo()],iCellI)
        {
            label cCell(cInfo.getIntCells()[Pstream::myProcNo()][iCellI]);

            const labelList& cFaces = mesh.cells()[cCell];

            forAll (cFaces,faceI)
            {
                if (!mesh.isInternalFace(cFaces[faceI]))
                {
                    // get reference to the patch which is in contact with IB. There is contact only if the patch is marked as a wall
                    label facePatchId(-1);
                    facePatchId = mesh.boundaryMesh().whichPatch(cFaces[faceI]);
                    const polyPatch& cPatch = mesh.boundaryMesh()[facePatchId];
                    if (cPatch.type()=="wall")
                    {
                        vector nVec(-mesh.Sf()[cFaces[faceI]]/mag(mesh.Sf()[cFaces[faceI]]));

                        if(!case3D)
                        {
                            if(mag(emptydir - nVec) < 1e-4)
                            {
                                continue;
                            }
                        }

                        labelList facePoints(mesh.faces()[cFaces[faceI]]);//list of vertex indicies
                        label numOfPoints = 0;

                        forAll(facePoints,pointI)
                        {
                            if(cInfo.getGeomModel().pointInside(mesh.points()[facePoints[pointI]]))
                            {
                                numOfPoints++;
                                center += mesh.points()[facePoints[pointI]];
                                normal += nVec;
                            }
                        }
                        centerPoints += numOfPoints;
                        area += mag(mesh.Sf()[cFaces[faceI]]) * numOfPoints / facePoints.size();
                    }
                }
            }
        }
        cInfo.inContactWithStatic(true);

        reduce(area, sumOp<scalar>());

        reduce(center, sumOp<vector>());
        reduce(normal, sumOp<vector>());
        reduce(centerPoints, sumOp<label>());
        center /= centerPoints;
        normal /= centerPoints;
        normal /= mag(normal);

        wallContactVars wallContactVar;

        wallContactVar.contactArea_ = area;
        wallContactVar.contactCenter_ = center;
        wallContactVar.contactVolume_ = getInterVolume<cT>(mesh,cInfo,cVars,geometricD,normal,center);
        wallContactVar.contactNormal_ = normal;
        cInfo.getWallContactVar().append(wallContactVar);
    }
    else
    {
        wallContactVars wallContactVar;

        wallContactVar.contactArea_ = 0;
        wallContactVar.contactCenter_ = vector::zero;
        wallContactVar.contactVolume_ = 0;
        wallContactVar.contactNormal_ = vector::zero;

        cInfo.getWallContactVar().append(wallContactVar);
    }
}
//---------------------------------------------------------------------------//
template <>
void detectWallContact <sphere>(
    const dynamicFvMesh&   mesh,
    contactInfo& cInfo,
    contactVars* cVars,
    vector geometricD
)
{
    cInfo.clearContactInfo();

    scalar cRadius(cInfo.getGeomModel().getDC() / 2);
    vector cCenter(cInfo.getGeomModel().getCoM());

    label nCells = mesh.nCells();

    bool case3D = true;
    vector emptydir = vector::zero;

    forAll (geometricD, direction)
    {
        if (geometricD[direction] == -1)
        {
            case3D = false;
            emptydir[direction] = 1;
            break;
        }
    }

    List<DynamicLabelList> contactFaces;

    DynamicVectorList center(Pstream::nProcs(),vector::zero);
    DynamicVectorList normal(Pstream::nProcs(),vector::zero);

    // go through all surfCells and check if there is any surfCell whose face is a boundary face
    forAll (cInfo.getSurfCells()[Pstream::myProcNo()],sCellI)
    {
        label cCell(cInfo.getSurfCells()[Pstream::myProcNo()][sCellI]);
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
                            if(mag(emptydir - nVec) < 1e-4 || mag(emptydir - (-1)*nVec) < 1e-4)
                            {
                                continue;
                            }
                        }

                        bool cont = false;

                        forAll(contactFaces,list)
                        {
                            forAll(contactFaces[list],face)
                            {
                                if(contactFaces[list][face] == cFaces[faceI])
                                {
                                    cont = true;
                                    break;
                                }
                            }

                            if(cont)
                            {
                                break;
                            }
                        }

                        if(cont)
                        {
                            continue;
                        }

                        face currentFace = mesh.faces()[cFaces[faceI]];
                        pointHit pointH = currentFace.nearestPoint(cCenter,mesh.points());
                        if(mag(pointH.rawPoint() - cCenter) < cRadius)
                        {
                            DynamicLabelList newContactFaces(1,cFaces[faceI]);
                            labelList nextToCheck(1,cFaces[faceI]);

                            while (nextToCheck.size() > 0)
                            {
                                DynamicLabelList auxToCheck;

                                forAll (nextToCheck,faceToCheck)
                                {
                                    labelList newfaces = getContactFaces(
                                                mesh,
                                                cInfo,
                                                nextToCheck[faceToCheck],
                                                newContactFaces,
                                                case3D,
                                                emptydir
                                    );

                                    forAll(newfaces,face)
                                    {
                                        newContactFaces.append(newfaces[face]);
                                        auxToCheck.append(newfaces[face]);
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
        point bPoint = vector::zero;
        forAll(contactFaces[list],faceI)
        {
            face currentFace = mesh.faces()[contactFaces[list][faceI]];
            pointHit pointH = currentFace.nearestPoint(cCenter,mesh.points());
            if(mag(pointH.rawPoint() - cCenter) < mag(bPoint - cCenter))
            {
                bPoint = pointH.rawPoint();
            }
        }

        vector nVec = (cCenter - bPoint)/mag(cCenter - bPoint);

        wallContactVars wallContactVar;

        wallContactVar.contactCenter_ = bPoint;
        wallContactVar.contactNormal_ = nVec;
        wallContactVar.contactArea_ = sphereContactArea(mesh,cInfo,cVars,geometricD,nVec,bPoint);
        wallContactVar.contactVolume_ = getInterVolume<sphere>(mesh,cInfo,cVars,geometricD,nVec,bPoint);

        cInfo.getWallContactVar().append(wallContactVar);

        if(wallContactVar.contactVolume_ > 0)
        {
            inContact = true;
        }
    }

    reduce(inContact, orOp<bool>());

    if(inContact)
    {
        cInfo.setWallContact(true);
        cInfo.inContactWithStatic(true);
    }
}
//---------------------------------------------------------------------------//
DynamicLabelList getContactFaces
(
    const dynamicFvMesh&   mesh,
    contactInfo& cInfo,
    label faceLabel,
    DynamicLabelList& checkedFaces,
    bool case3D,
    vector& emptydir
)
{
    scalar cRadius(cInfo.getGeomModel().getDC() / 2);
    vector cCenter(cInfo.getGeomModel().getCoM());

    DynamicLabelList facesReturnList;

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
                if (cPatch.type()=="wall")
                {
                    vector nVec(-mesh.Sf()[edgeFaces[faceI]]/mag(mesh.Sf()[edgeFaces[faceI]]));

                    if(!case3D)
                    {
                        if(mag(emptydir - nVec) < 1e-4 || mag(emptydir - (-1)*nVec) < 1e-4)
                        {
                            continue;
                        }
                    }

                    bool cont = false;

                    forAll(checkedFaces,face)
                    {
                        if(checkedFaces[face] == edgeFaces[faceI])
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
                        facesReturnList.append(edgeFaces[faceI]);
                    }
                }
            }
        }
    }

    return facesReturnList;
}
//---------------------------------------------------------------------------//
void solveWallContact
(
    wallInfo& wInfo,
    contactInfo& cInfo,
    contactVars* cVars,
    scalar deltaT,
    Tuple2<vector,vector>& outVars
)
{
    // compute mean model parameters
    scalar aY = 1/((1 - sqr(cInfo.getNu()))/cInfo.getY()
                + (1 - sqr(wInfo.getNu()))/wInfo.getY());                // Contact Young modulus
    scalar aG = 1/(2*(2 - cInfo.getNu())*(1 + cInfo.getNu())/cInfo.getY()
                + 2*(2 - wInfo.getNu())*(1 + wInfo.getNu())/wInfo.getY()); // Contact shear modulus
    scalar aGammaN = aY*(cInfo.getGamma()*wInfo.getGamma())
                    /(cInfo.getGamma()+wInfo.getGamma()+SMALL);         // Contact normal gamma
    scalar aGammat = aG*(cInfo.getGamma()*wInfo.getGamma())
                    /(cInfo.getGamma()+wInfo.getGamma()+SMALL);         // Contact tangential gamma
    scalar aMu = (cInfo.getmu()+wInfo.getmu())/2;                       // Contact frictional coef
    scalar aadhN = (cInfo.getAdhN()+wInfo.getAdhN())/2;                 // Contact adhesive coef

    vector FnormOut = vector::zero;
    vector FtanOut = vector::zero;
    vector cLVecOut = vector::zero;
    bool FtLastFinded(false);

    forAll(cInfo.getWallContactVar(),contVar)
    {

        scalar intersectedVolume = cInfo.getWallContactVar()[contVar].contactVolume_;
        scalar overallContactArea = cInfo.getWallContactVar()[contVar].contactArea_;
        vector nVec = cInfo.getWallContactVar()[contVar].contactNormal_;
        vector cLVec = cInfo.getWallContactVar()[contVar].contactCenter_
                       - cInfo.getGeomModel().getCoM();

        if(intersectedVolume == 0)
        {
            continue;
        }

        Info << "-- Detected Particle-wall contact: -- body " << cVars->bodyId_ << endl;
        Info << "-- Particle-wall contact center " << cInfo.getWallContactVar()[contVar].contactCenter_ << endl;
        Info << "-- Particle-wall contact normal " << cInfo.getWallContactVar()[contVar].contactNormal_ << endl;
        Info << "-- Particle-wall contact volume " << cInfo.getWallContactVar()[contVar].contactVolume_ << endl;
        Info << "-- Particle-wall contact area "   << cInfo.getWallContactVar()[contVar].contactArea_ << endl;

        scalar Lc(4*mag(cLVec)*mag(cLVec)/(mag(cLVec)+mag(cLVec)));
        scalar reduceM(cVars->M0_*cVars->M0_/(cVars->M0_+cVars->M0_));

        vector planarVec       =  cLVec - cVars->Axis_*(cLVec & cVars->Axis_);
        vector rotDir(planarVec^cVars->Axis_);

        vector cVel(-rotDir*cVars->omega_ + cVars->Vel_);
        vector wVel(vector::zero);
        scalar Vn(-(cVel-wVel) & nVec);

        vector FN = (aY*intersectedVolume/(Lc+SMALL)
            + aGammaN*sqrt(aY*reduceM/pow(Lc+SMALL,3))
            *(Vn*overallContactArea))*nVec;

        Info << "-- Particle-wall contact FN " << FN << endl;

        // if the IB was in contact in previous DEM time step, find the information about tangential force and assigne it
        vector FtLast(vector::zero);
        forAll (cInfo.getHistoryFt(),Fti)
        {
            if (cInfo.getHistoryFt()[Fti].first() == -1)
            {
                FtLastFinded = true;
                FtLast = cInfo.getHistoryFt()[Fti].second().second();
                break;
            }
        }

        // project last Ft to new tangential direction
        vector FtLastP(FtLast - (FtLast & nVec) * nVec);
        // scale the projected vector to remain the magnitude
        vector FtLastr(mag(FtLast) * (FtLastP/(mag(FtLastP)+SMALL)));

        vector cVeliNorm((cVel & nVec)*nVec);
        vector Vt = cVel-cVeliNorm-wVel;
        vector Ft = FtLastr - aG*Lc*Vt*deltaT - aGammat*sqrt(aG*reduceM*Lc)*Vt;

        Info << "-- Particle-wall contact Ft " << Ft << endl;
        if (mag(Ft) > aMu * mag(FN))
        {
            Ft *= aMu * mag(FN) / mag(Ft);
        }
        Info << "-- Particle-wall contact Ft clamped" << Ft << endl;

        scalar pi = Foam::constant::mathematical::pi;
        vector FA((sqrt(8*pi*aY*aadhN*intersectedVolume)) * nVec);

        Info << "-- Particle-wall FA " << FA << endl;
        FN -= FA;

        FnormOut += FN;
        FtanOut += Ft;
        cLVecOut += cLVec;
    }

    reduce(FnormOut, sumOp<vector>());
    reduce(FtanOut, sumOp<vector>());
    reduce(cLVecOut, sumOp<vector>());

    if(mag(FnormOut) == 0)
    {
        return;
    }

    // update or add the history of tangential force
    if (FtLastFinded)
    {
        forAll (cInfo.getHistoryFt(),Fti)
        {
            if (cInfo.getHistoryFt()[Fti].first() == -1)
            {
                Tuple2<label,vector> help(1,FtanOut);
                cInfo.getHistoryFt()[Fti].second() = help;
                break;
            }
        }
    }
    else
    {
        Tuple2<label,vector> help(1,FtanOut);
        Tuple2<label,Tuple2<label,vector>> help2(-1, help);

        cInfo.getHistoryFt().append(help2);
    }

    outVars.first() = FnormOut + FtanOut;
    outVars.second() = cLVecOut ^ (FnormOut + FtanOut);
}
//---------------------------------------------------------------------------//
template <contactType cT>
scalar
getInterVolume(
    const dynamicFvMesh&   mesh,
    contactInfo& cInfo,
    contactVars* cVars,
    vector geometricD,
    vector nVec,
    vector center
)
{
    return (cVars->M0_-cVars->M_)/(cVars->rhoS_.value() + SMALL);
}
//---------------------------------------------------------------------------//
template <>
scalar
getInterVolume <sphere>(
    const dynamicFvMesh&   mesh,
    contactInfo& cInfo,
    contactVars* cVars,
    vector geometricD,
    vector nVec,
    vector center
)
{
    scalar cRadius(cInfo.getGeomModel().getDC() / 2);
    vector cCenter(cInfo.getGeomModel().getCoM());

    plane contPlane(center, nVec);
    scalar dist = contPlane.distance(cCenter);
    scalar xH = cRadius - dist;

    bool case3D = true;
    label emptyDim = 0;
    // check if the case is 3D
    forAll (geometricD, direction)
    {
        if (geometricD[direction] == -1)
        {
            case3D = false;
            emptyDim = direction;
            break;
        }
    }

    if(case3D)
    {
        scalar pi = Foam::constant::mathematical::pi;
        return (pi*sqr(xH)*(3*cRadius - xH)) / 3;
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
    const dynamicFvMesh&   mesh,
    contactInfo& cInfo,
    contactVars* cVars,
    vector geometricD,
    vector nVec,
    vector center
)
{
    scalar cRadius(cInfo.getGeomModel().getDC() / 2);
    vector cCenter(cInfo.getGeomModel().getCoM());

    plane contPlane(center, nVec);

    scalar pi = Foam::constant::mathematical::pi;
    scalar dist = contPlane.distance(cCenter);
    scalar contactRad = sqrt(sqr(cRadius) - sqr(dist));

    bool case3D = true;
    label emptyDim = 0;
    // check if the case is 3D
    forAll (geometricD, direction)
    {
        if (geometricD[direction] == -1)
        {
            case3D = false;
            emptyDim = direction;
            break;
        }
    }

    if(case3D)
    {
        return pi*sqr(contactRad);
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
