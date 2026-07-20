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
#include "nonConvexBody.H"

using namespace Foam;

//---------------------------------------------------------------------------//
// auxiliary for createImmersedBody
labelList nonConvexBody::getBBoxCellsByOctTree
(
    label cellToCheck,
    bool& insideBB,
    vector& bBoxMin,
    vector& bBoxMax,
    List<DynamicLabelList>& bBoxCells,
    Field<label>& octreeField
)
{
    labelList retList;

    if (octreeField[cellToCheck] ==0)
    {
        octreeField[cellToCheck] = 1;
        vector cCenter = mesh_.C()[cellToCheck];
        label   partCheck(0);
        forAll (bBoxMin,vecI)
        {
            if (cCenter[vecI] >= bBoxMin[vecI] && cCenter[vecI] <= bBoxMax[vecI])
            {
                partCheck++;
            }
        }
        bool cellInside = (partCheck == 3) ? true : false;
        if (cellInside)
        {
            bBoxCells[Pstream::myProcNo()].append(cellToCheck);
            insideBB = true;
        }
        if (not insideBB or cellInside)
        {
            retList = mesh_.cellCells()[cellToCheck];
        }
    }
    return retList;
}
//---------------------------------------------------------------------------//
// create immersed body for convex body
void nonConvexBody::createImmersedBody
(
    volScalarField& body,
    Field<label>& octreeField,
    List<labelList>& cellPoints
)
{
    // clear old list contents
    intCells_[Pstream::myProcNo()].clear();
    surfCells_[Pstream::myProcNo()].clear();

    // reduce computational domain to the body bounding box
    scalar inflFact(2*sqrt(mesh_.magSf()[0]));

    boundBox ibBound(getBounds());

    vector expMinBBox = ibBound.min() - vector::one*inflFact;
    vector expMaxBBox = ibBound.max() + vector::one*inflFact;

    octreeField *= 0;
    List<DynamicLabelList> bBoxCells(Pstream::nProcs());

    // prepare for sending
    autoPtr<List<DynamicLabelList>> neighboursToSend(
        new List<DynamicLabelList>(Pstream::nProcs()));

    // get number of empty directions
    label nEmpty(0);
    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& cPatch = mesh_.boundaryMesh()[patchI];
        if (cPatch.type() == "empty")
        {
            nEmpty += 1;
        }
    }

    bool isInsideBB(false);
    labelList nextToCheck(1,0);
    label iterCount(0);label iterMax(mesh_.nCells());

    label nextSize = nextToCheck.size();
    reduce(nextSize, maxOp<label>());
    while ((nextSize > 0 or not isInsideBB) && iterCount < iterMax)
    {
        iterCount++;
        DynamicLabelList auxToCheck;

        forAll (nextToCheck,cellToCheck)
        {
            auxToCheck.append(
                getBBoxCellsByOctTree(
                    nextToCheck[cellToCheck],
                    isInsideBB,
                    expMinBBox,
                    expMaxBBox,
                    bBoxCells,
                    octreeField
                )
            );

            if (!isInsideBB)
            {
                continue;
            }

            // get cell label
            label cCell = nextToCheck[cellToCheck];

            // if number of neighbours equals number of faces, skip processors 
            label nProcFaces = mesh_.cells()[cCell].size() - mesh_.cellCells()[cCell].size();
            nProcFaces -= nEmpty;
            if(nProcFaces == 0)
            {
                continue;
            }

            // add processor neighbors
            DynamicList<labelList> procFaces;
            forAll(mesh_.cells()[cCell], fI)
            {
                // get face label
                label faceI = mesh_.cells()[cCell][fI];

                // check if face is a processor face
                if(!mesh_.isInternalFace(faceI))
                {
                    // get the patch the face belongs to
                    label facePatchI(mesh_.boundaryMesh().whichPatch(faceI));
                    const polyPatch& cPatch = mesh_.boundaryMesh()[facePatchI];

                    // check if it is a processor boundary
                    if (cPatch.type() == "processor")
                    {
                        // get the processor patch
                        const processorPolyPatch& procPatch
                            = refCast<const processorPolyPatch>(cPatch);

                        // get the neighboring processor id
                        label iProc = (Pstream::myProcNo() == procPatch.myProcNo())
                            ? procPatch.neighbProcNo() : procPatch.myProcNo();

                        // get local face value
                        label iFace = cPatch.whichFace(faceI);

                        // save to send
                        procFaces.append({iProc, iFace});
                        neighboursToSend()[iProc].append(iFace);
                    }
                }
            }
        }

        // send processor neighbors to add to next to check 
        PstreamBuffers pBufsIFaces(Pstream::commsTypes::nonBlocking);
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            if(proci != Pstream::myProcNo())
            {
                UOPstream sendIFaces(proci, pBufsIFaces);
                sendIFaces << neighboursToSend()[proci];
                neighboursToSend()[proci].clear();
            }
        }

        pBufsIFaces.finishedSends();

        // recieve and add to aux to check
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            if (proci != Pstream::myProcNo())
            {
                UIPstream recvIFaces(proci, pBufsIFaces);
                DynamicList<label> recIFaces (recvIFaces);

                // find cells for faces
                forAll(recIFaces, rFace)
                {
                    // get the cell label
                    label faceI = recIFaces[rFace]; // local face labels

                    // find the respective cell 
                    forAll(mesh_.boundaryMesh(), patchI)
                    {
                        if (isA<processorPolyPatch>(mesh_.boundaryMesh()[patchI]))
                        {
                            const processorPolyPatch& procPatch
                                = refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchI]);

                            // get the neighboring processor id
                            label iProc = (Pstream::myProcNo() == procPatch.myProcNo())
                                ? procPatch.neighbProcNo() : procPatch.myProcNo();

                            if (iProc == proci)
                            {
                                // get the cell label
                                label rCellI = mesh_.boundaryMesh()[patchI].faceCells()[faceI];
                                auxToCheck.append(
                                    getBBoxCellsByOctTree(
                                        rCellI,
                                        isInsideBB,
                                        expMinBBox,
                                        expMaxBBox,
                                        bBoxCells,
                                        octreeField
                                    )
                                );
                            }
                        }
                    }
                }
            }
        }

        // clear buffer
        pBufsIFaces.clear();

        // next iteration
        nextToCheck = auxToCheck;

        // check if all processors finished 
        nextSize = nextToCheck.size();
        reduce(nextSize, maxOp<label>());
    }

    // get cell centers inside the body bounding box
    const pointField& cp = mesh_.C();
    const pointField fCp = filterField(cp,bBoxCells[Pstream::myProcNo()]);
    boolList fCentersInside = triSurfSearch_().calcInside(fCp);

    // find the processor with most of this IB inside
    ibPartialVolume_[Pstream::myProcNo()] = 0;

    vector sDSpan(4.0*(mesh_.bounds().max()-mesh_.bounds().min()));

    // first loop, construction of body field and identification of
    // the number of inside and surface cells
    const pointField& pp = mesh_.points();
    forAll (bBoxCells[Pstream::myProcNo()],bCellI)                      //go only through bBox
    {
        label cellI(bBoxCells[Pstream::myProcNo()][bCellI]);

        // check if partially or completely inside
        const pointField vertexPoints(pp,cellPoints[cellI]);
        boolList vertexesInside = triSurfSearch_().calcInside( vertexPoints );
        bool centerInside(fCentersInside[bCellI]);
        scalar rVInSize(0.5/vertexesInside.size());
        // Note: weight of a single vertex in the cell

        scalar cBody(0);
        forAll (vertexesInside, verIn)
        {
            if (vertexesInside[verIn]==true)
            {
                cBody  += rVInSize;                                     //fraction of cell covered
            }
        }

        if (centerInside)
        {
            cBody+=0.5;
        }
        if (cBody > thrSurf_)
        {
            if (cBody > (1.0-thrSurf_))
            {
                intCells_[Pstream::myProcNo()].append(cellI);
                cellToStartInCreateIB_ = cellI;
            }
            else if (cBody  <= (1.0-thrSurf_))
            {
                surfCells_[Pstream::myProcNo()].append(cellI);
                if (sdBasedLambda_)
                {
                    pointIndexHit pointHit(
                        triSurfSearch_().nearest(
                            mesh_.C()[cellI],
                            sDSpan
                        )
                    );
                    scalar signedDist(0.0);
                    if (pointHit.hit())
                    {
                        signedDist = mag(pointHit.hitPoint()-cp[cellI]);
                    }
                    else
                    {
                        InfoH << iB_Info
                            << "Missed point in signedDist computation !!"
                            << endl;
                    }
                    if (centerInside)
                    {
                        cBody = 0.5*(Foam::tanh(intSpan_*signedDist/Foam::pow(mesh_.V()[cellI],0.333))+1.0);
                    }
                    else
                    {
                        cBody = 0.5*(-1.0*Foam::tanh(intSpan_*signedDist/Foam::pow(mesh_.V()[cellI],0.333))+1.0);
                    }
                }
            }
            ibPartialVolume_[Pstream::myProcNo()] += 1;
        }
        body[cellI]+= cBody;
        // clip the body field values
        body[cellI] = min(max(0.0,body[cellI]),1.0);
    }

    // gather partial volume from other processors
    Pstream::gatherList(ibPartialVolume_, 0);
    Pstream::broadcast(ibPartialVolume_, 0);                            //OF.com: scatter -> broadcast

    for (label i = 0; i < ibPartialVolume_.size(); i++)
    {
        if (ibPartialVolume_[i] == max(ibPartialVolume_))
        {
        //Set owner of the IB which will move this IB
            owner_ = i;
            break;
        }
    }
}
//---------------------------------------------------------------------------//
