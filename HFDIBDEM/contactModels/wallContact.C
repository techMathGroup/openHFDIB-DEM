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
                    labelList facePoints(mesh.faces()[cFaces[faceI]]);//list of vertex indicies
                    label numOfPoints = 0;

                    forAll(facePoints,pointI)
                    {
                        if(cInfo.getGeomModel().pointInside(mesh.points()[facePoints[pointI]]))
                        {
                            numOfPoints++;
                            center += mesh.points()[facePoints[pointI]];
                            normal += -mesh.Sf()[cFaces[faceI]]/mag(mesh.Sf()[cFaces[faceI]]);
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
                        labelList facePoints(mesh.faces()[cFaces[faceI]]);//list of vertex indicies
                        label numOfPoints = 0;

                        forAll(facePoints,pointI)
                        {
                            if(cInfo.getGeomModel().pointInside(mesh.points()[facePoints[pointI]]))
                            {
                                numOfPoints++;
                                center += mesh.points()[facePoints[pointI]];
                                normal += -mesh.Sf()[cFaces[faceI]]/mag(mesh.Sf()[cFaces[faceI]]);
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

    while(true)
    {
        bool inContact(false);
        DynamicVectorList center(Pstream::nProcs(),vector::zero);
        DynamicVectorList normal(Pstream::nProcs(),vector::zero);

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
                        plane contPlane(mesh.Cf()[cFaces[faceI]], nVec);
                        scalar dist = contPlane.distance(cCenter);
                        scalar xH = cRadius - dist;
                        if(xH > 0)
                        {
                            bool newVec(true);
                            forAll(cInfo.getWallContactVar(),contVar)
                            {
                                if(cInfo.getWallContactVar()[contVar].contactNormal_ == nVec)
                                {
                                    newVec = false;
                                    break;
                                }
                            }

                            if(newVec)
                            {
                                center[Pstream::myProcNo()] = contPlane.nearestPoint(cCenter);
                                normal[Pstream::myProcNo()] = nVec;
                                inContact = true;
                                break;
                            }
                        }
                    }
                }

                if(inContact)
                {
                    break;
                }
            }
        }

        reduce(inContact, orOp<bool>());

        if(inContact)
        {
            cInfo.setWallContact(true);

            Pstream::gatherList(center, 0);
            Pstream::scatter(center, 0);

            Pstream::gatherList(normal, 0);
            Pstream::scatter(normal, 0);

            Info << center << endl;
            Info << normal << endl;

            label contProc(0);

            forAll(normal,ci)
            {
                if(normal[ci] != vector::zero)
                {
                    contProc = ci;
                    break;
                }
            }

            wallContactVars wallContactVar;

            wallContactVar.contactCenter_ = center[contProc];
            wallContactVar.contactNormal_ = normal[contProc];
            wallContactVar.contactArea_ = sphereContactArea(mesh,cInfo,cVars,geometricD,normal[contProc],center[contProc]);
            wallContactVar.contactVolume_ = getInterVolume<sphere>(mesh,cInfo,cVars,geometricD,normal[contProc],center[contProc]);

            cInfo.getWallContactVar().append(wallContactVar);
        }
        else
        {
            Info << "Num of contacts: " << cInfo.getWallContactVar().size() << endl;
            break;
        }
    }
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
    scalar aKN((cInfo.getkN()*wInfo.getkN())/(cInfo.getkN()+wInfo.getkN()+SMALL));
    scalar aGammaN((cInfo.getgammaN()*wInfo.getgammaN())/(cInfo.getgammaN()+wInfo.getgammaN()+SMALL));
    scalar aKt((cInfo.getkt()*wInfo.getkt())/(cInfo.getkt()+wInfo.getkt()+SMALL));
    scalar aGammat((cInfo.getgammat()*wInfo.getgammat())/(cInfo.getgammat()+wInfo.getgammat()+SMALL));
    scalar amu((cInfo.getmu()*wInfo.getmu())/(cInfo.getmu()+wInfo.getmu()+SMALL));
    scalar aadhN((cInfo.getAdhN()*wInfo.getAdhN())/(cInfo.getAdhN()+wInfo.getAdhN()+SMALL));
    scalar aadhEqui(0.5*(cInfo.getAdhEqui()+wInfo.getAdhEqui()));

    vector FnormOut = vector::zero;
    vector FtanOut = vector::zero;
    vector cLVecOut = vector::zero;
    bool FtLastFinded(false);

    forAll(cInfo.getWallContactVar(),contVar)
    {
        Info << "-// Wall contact num: " << (contVar + 1) << endl;

        scalar intersectedVolume = cInfo.getWallContactVar()[contVar].contactVolume_;
        scalar overallContactArea = cInfo.getWallContactVar()[contVar].contactArea_;
        vector nVec = cInfo.getWallContactVar()[contVar].contactNormal_;
        vector cLVec = cInfo.getWallContactVar()[contVar].contactCenter_
                       - cInfo.getGeomModel().getCoM();
        Info << "-// Volume: " << intersectedVolume << endl;
        Info << "-// Area: " << overallContactArea << endl;
        Info << "-// nVec: " << nVec << endl;

        if(intersectedVolume == 0)
        {
            return;
        }

        scalar Lc(4*mag(cLVec)*mag(cLVec)/(mag(cLVec)+mag(cLVec)));
        scalar reduceM(cVars->M0_*cVars->M0_/(cVars->M0_+cVars->M0_));

        vector planarVec       =  cLVec - cVars->Axis_*(cLVec & cVars->Axis_);
        vector rotDir(planarVec^cVars->Axis_);

        vector cVel(-rotDir*cVars->omega_ + cVars->Vel_);
        vector wVel(vector::zero);
        scalar VnF(-(cVel-wVel) & nVec);

        vector FN = (aKN*intersectedVolume/(Lc+SMALL)
            + aGammaN*sqrt(aKN*reduceM/pow(Lc+SMALL,3))
            *(VnF*overallContactArea))*nVec;

        Info << "-// VnF: " << VnF << endl;
        Info << "-// Fn: " << FN << endl;
        Info << "-// Fn norm: " << aKN*intersectedVolume/(Lc+SMALL) << endl;
        Info << "-// Fn dump: " << (aGammaN*sqrt(aKN*reduceM/pow(Lc+SMALL,3))*(VnF*overallContactArea)) << endl;

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
        vector Ft = FtLastr - aKt*Lc*Vt*deltaT - aGammat*sqrt(aKt*reduceM*Lc)*Vt;

        Info << "-// Vt: " << Vt << endl;
        Info << "-// AngVel: " << cVars->omega_ << endl;
        Info << "-// Ft: " << Ft << endl;

        if (mag(Ft) > amu * mag(FN))
        {
            Ft *= amu * mag(FN) / mag(Ft);
        }

        Info << "-// Ft small: " << Ft << endl;

        scalar FAc(aadhN*overallContactArea);
        scalar FAeq(aKN*((aadhEqui*cVars->M0_)/(cVars->rhoS_.value() + SMALL))/(Lc+SMALL));
        scalar partMul((cVars->M0_-cVars->M_)/(cVars->M0_+SMALL)/(aadhEqui+SMALL));
        if(partMul > 1)
        {
            partMul = 1;
        }
        vector FA((FAeq * partMul  + FAc * (1-partMul)) * nVec);
        FN -= FA;

        FnormOut += FN;
        FtanOut += Ft;
        cLVecOut += cLVec;
    }

    Info << "FnormOut: " << FnormOut << endl;
    Info << "FtanOut: " << FtanOut << endl;
    Info << "cLVecOut: " << cLVecOut << endl;

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

    Info << "cRadius " << cRadius << " cCenter " << cCenter << endl;
    plane contPlane(center, nVec);
    scalar dist = contPlane.distance(cCenter);
    Info << "dist: " << dist << endl;
    scalar xH = cRadius - dist;
    Info << "xH " << xH << endl;

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

        Info << "emptyLength" << emptyLength << endl;
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
