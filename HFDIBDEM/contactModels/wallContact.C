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
    contactInfo& cInfo
)
{
    bool inContact(false);
    cInfo.clearContactInfo();

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

                    forAll(facePoints,pointI)
                    {
                        if(cInfo.getGeomModel().pointInside(mesh.points()[facePoints[pointI]]))
                        {
                            cInfo.wallContactFaces()[Pstream::myProcNo()].append(cFaces[faceI]);
                            inContact = true;
                            break;
                        }
                    }
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

                        forAll(facePoints,pointI)
                        {
                            if(cInfo.getGeomModel().pointInside(mesh.points()[facePoints[pointI]]))
                            {
                                cInfo.wallContactFaces()[Pstream::myProcNo()].append(cFaces[faceI]);
                                break;
                            }
                        }
                    }
                }
            }
        }
        cInfo.inContactWithStatic(true);
        cInfo.distributeWallContactFaces();
    }
}
//---------------------------------------------------------------------------//
void solveWallContact
(
    const dynamicFvMesh&   mesh,
    wallInfo& wInfo,
    contactInfo& cInfo,
    contactVars* cVars,
    vector geometricD,
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

    label nContactFaces(cInfo.wallContactFaces()[Pstream::myProcNo()].size());
    // create placeholders for forces
    vector FN(vector::zero);
    vector Ft(vector::zero);
    vector cLVec(vector::zero);
    vector nVecF(vector::zero);
    scalar overallContactArea(0);

    // if the IB was in contact in previous DEM time step, find the information about tangential force and assigne it
    vector FtLast(vector::zero);
    bool FtLastFinded(false);
    forAll (cInfo.getHistoryhistoryFt(),Fti)
    {
        if (cInfo.getHistoryhistoryFt()[Fti].first() == -1)
        {
            FtLastFinded = true;
            FtLast = cInfo.getHistoryhistoryFt()[Fti].second().second();
            break;
        }
    }

    DynamicLabelList contactCells;
    vector cVel(vector::zero);
    scalar reduceM(cVars->M0_*cVars->M0_/(cVars->M0_+cVars->M0_));

    // loop over all faces in contact with walls
    forAll (cInfo.wallContactFaces()[Pstream::myProcNo()],faceI)
    {
        label cFace(cInfo.wallContactFaces()[Pstream::myProcNo()][faceI]);

        contactCells.append(mesh.faceOwner()[cFace]);
        // get the local wall velocity
        vector wVel(vector::zero);

        // compute normal to movement and relative velocitydeltaTDEM
        vector nVec(-mesh.Sf()[cFace]/mag(mesh.Sf()[cFace]));

        nVecF += nVec * mag(mesh.Sf()[cFace]);
        overallContactArea += mag(mesh.Sf()[cFace]);

        // project last Ft to new tangential direction
        vector FtLastP(FtLast - (FtLast & nVec) * nVec);
        // scale the projected vector to remain the magnitude
        vector FtLastr(mag(FtLast) * (FtLastP/(mag(FtLastP)+SMALL)));

        vector ri(mesh.Cf()[cFace] - cInfo.getGeomModel().getCoM());
        // evaluate tangential velocity
        vector planarVec       =  ri - cVars->Axis_*(ri & cVars->Axis_);
        vector rotDir(planarVec^cVars->Axis_);

        vector cVeli(-rotDir*cVars->omega_ + cVars->Vel_);
        cVel += cVeli;
        vector cVeliNorm((cVeli & nVec)*nVec);
        cVeli -= cVeliNorm;
        vector Vt(cVeli-wVel);
        // compute tangential force
        vector Fti(FtLastr - aKt*Vt*deltaT);

        scalar Lci(4*mag(ri)*mag(ri)/(mag(ri)+mag(ri)));

        vector Ftdi(- aGammat*sqrt(aKN*reduceM*Lci)*Vt);
        Ft += (Fti + Ftdi);
        cLVec += ri;
    }


    reduce(nContactFaces, sumOp<label>());
    reduce(overallContactArea, sumOp<scalar>());
    reduce(Ft, sumOp<vector>());
    reduce(cVel, sumOp<vector>());
    reduce(cLVec, sumOp<vector>());
    reduce(nVecF, sumOp<vector>());

    vector wVel(vector::zero);
    cLVec /= (nContactFaces+SMALL);
    cVel /= (nContactFaces+SMALL);
    nVecF /= (overallContactArea+SMALL);
    scalar VnF(-(cVel-wVel) & nVecF);
    Ft /= (nContactFaces+SMALL);

    Info << "-// Ft: " << Ft << endl;

    scalar intersectedVolume;
    if(cInfo.getGeomModel().getcType() == sphere)
    {
        intersectedVolume = getInterVolume<sphere>(mesh,cInfo,cVars,geometricD);
        overallContactArea = sphereContactArea(mesh,cInfo,cVars,geometricD);
    }
    else
    {
        intersectedVolume = getInterVolume<arbShape>(mesh,cInfo,cVars,geometricD);
    }
    reduce(intersectedVolume,maxOp<scalar>());
    reduce(overallContactArea, maxOp<scalar>());

    Info << "-// Volume: " << intersectedVolume << endl;
    Info << "-// Area: " << overallContactArea << endl;

    scalar Lc(4*mag(cLVec)*mag(cLVec)/(mag(cLVec)+mag(cLVec)));

    FN = (aKN*intersectedVolume/(Lc+SMALL) + aGammaN*sqrt(aKN*reduceM/pow(Lc+SMALL,3))*(VnF*overallContactArea))*nVecF;

    if (mag(Ft) > amu * mag(FN))
    {
        Ft *= amu * mag(FN) / mag(Ft);
    }

    Info << "-// Ft: " << Ft << endl;

    scalar FAc(aadhN*overallContactArea);
    scalar FAeq(aKN*((aadhEqui*cVars->M0_)/(cVars->rhoS_.value() + SMALL))/(Lc+SMALL));
    scalar partMul((cVars->M0_-cVars->M_)/(cVars->M0_+SMALL)/(aadhEqui+SMALL));
    if(partMul > 1)
    {
        partMul = 1;
    }
    vector FA((FAeq * partMul  + FAc * (1-partMul)) * nVecF);
    FN -= FA;

    // update or add the history of tangential force
    if (FtLastFinded)
    {
        forAll (cInfo.getHistoryhistoryFt(),Fti)
        {
            if (cInfo.getHistoryhistoryFt()[Fti].first() == -1)
            {
                Tuple2<label,vector> help(1,Ft);
                cInfo.getHistoryhistoryFt()[Fti].second() = help;
                break;
            }
        }
    }
    else
    {
        Tuple2<label,vector> help(1,Ft);
        Tuple2<label,Tuple2<label,vector>> help2(-1, help);

        cInfo.getHistoryhistoryFt().append(help2);
    }

    outVars.first() = FN + Ft;
    outVars.second() = cLVec ^ (FN + Ft);
}
//---------------------------------------------------------------------------//
template <contactType cT>
scalar
getInterVolume(
    const dynamicFvMesh&   mesh,
    contactInfo& cInfo,
    contactVars* cVars,
    vector geometricD
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
    vector geometricD
)
{
    scalar cRadius(cInfo.getGeomModel().getDC() / 2);
    vector cCenter(cInfo.getGeomModel().getCoM());
    if(cInfo.wallContactFaces()[Pstream::myProcNo()].size() == 0)
    {
        return 0;
    }

    label cFace(cInfo.wallContactFaces()[Pstream::myProcNo()][0]);

    Info << "cRadius " << cRadius << " cCenter " << cCenter << endl;
    Info << "contactWallsSize: " << cInfo.wallContactFaces()[Pstream::myProcNo()].size() << endl;
    vector nVec(-mesh.Sf()[cFace]/mag(mesh.Sf()[cFace]));
    Info  << "nVec " << nVec << endl;
    plane contPlane(mesh.Cf()[cFace], nVec);
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
    vector geometricD
)
{
    scalar cRadius(cInfo.getGeomModel().getDC() / 2);
    vector cCenter(cInfo.getGeomModel().getCoM());
    if(cInfo.wallContactFaces()[Pstream::myProcNo()].size() == 0)
    {
        return 0;
    }

    label cFace(cInfo.wallContactFaces()[Pstream::myProcNo()][0]);
    vector nVec(-mesh.Sf()[cFace]/mag(mesh.Sf()[cFace]));
    plane contPlane(mesh.Cf()[cFace], nVec);

    scalar pi = Foam::constant::mathematical::pi;
    scalar dist = contPlane.distance(cCenter);
    scalar contactRad = sqr(cRadius) - sqr(dist);

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
