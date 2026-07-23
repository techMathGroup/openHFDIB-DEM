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
#include "convexBody.H"

using namespace Foam;

//---------------------------------------------------------------------------//
// create immersed body for convex body
void convexBody::createImmersedBody
(
    volScalarField& body,
    Field<label>& octreeField,
    List<labelList>& cellPoints
)
{
    // clear old list contents
    intCells_[Pstream::myProcNo()].clear();
    surfCells_[Pstream::myProcNo()].clear();
    haloCells_[Pstream::myProcNo()].clear();

    // find the processor with most of this IB inside
    ibPartialVolume_[Pstream::myProcNo()] = 0;

    label nextSize = 1;
    if(!isBBoxInMesh())
    {
        nextSize = 0;
    }

    label cellInIB = getCellInBody(octreeField);

    if(cellInIB == -1)
    {
        nextSize = 0;
    }

    // get the list of cell centroids
    const pointField& cp = mesh_.C();

    autoPtr<DynamicLabelList> nextToCheck(
        new DynamicLabelList(nextSize,cellInIB));
    autoPtr<DynamicLabelList> auxToCheck(
        new DynamicLabelList);
    autoPtr<List<DynamicLabelList>> neighboursToSend(
        new List<DynamicLabelList>(Pstream::nProcs()));

    label tableSize = 128;
    if(cachedNeighbours_.valid() && getRefineBuffers() <= 0)
    {
        tableSize = cachedNeighbours_().toc().size()*1.5;
    }
    else
    {
        cachedNeighbours_.reset(new HashTable<labelList, label, Hash<label>>);
    }

    // initialize cache for processor neighbours
    if(!procNeighbours_.valid())
    {
        procNeighbours_.reset(new HashTable<DynamicList<labelList>, label, Hash<label>>);
    }

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

    HashTable<bool, label, Hash<label>> cellInside(tableSize);

    label iterCount(0);label iterMax(mesh_.nCells());
    reduce(nextSize, maxOp<label>());
    while (nextSize > 0 and iterCount++ < iterMax)
    {
        // clear neighbors found in previous iteration
        auxToCheck().clear();

        // loop over neighbors to check found in previous iteration
        forAll (nextToCheck(),cellToCheck)
        {
            // get cell label
            label cCell = nextToCheck()[cellToCheck];

            // continue if it was not visited
            if (!cellInside.found(cCell))
            {
                iterCount++;

                // if inside body add neighbors to check
                if(pointInside(cp[cCell]))
                {
                    cellInside.set(cCell, true);

                    if(cachedNeighbours_.valid() && cachedNeighbours_().found(cCell))
                    {
                        auxToCheck().append(cachedNeighbours_()[cCell]);
                    }
                    else
                    {
                        const labelList& neigh = mesh_.cellCells(cCell);
                        cachedNeighbours_().insert(cCell, neigh);
                        auxToCheck().append(neigh);
                    }

                    // if number of neighbours equals number of faces, skip processors 
                    label nProcFaces = mesh_.cells()[cCell].size() - mesh_.cellCells()[cCell].size();
                    nProcFaces -= nEmpty;
                    if(nProcFaces == 0)
                    {
                        continue;
                    }

                    // check processor neighbours
                    if(procNeighbours_.valid() && procNeighbours_().found(cCell))
                    {
                        DynamicList<labelList>& procFaces = procNeighbours_()[cCell];
                        forAll(procFaces, pI)
                        {
                            neighboursToSend()[procFaces[pI][0]].append(procFaces[pI][1]);
                        }
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
                    procNeighbours_().insert(cCell, procFaces);
                }
                else
                {
                    cellInside.set(cCell, false);
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
                                auxToCheck().append(rCellI);
                            }
                        }
                    }
                }
            }
        }

        // clear buffer
        pBufsIFaces.clear();

        // clean up and prep for next iter
        autoPtr<DynamicLabelList> helpPtr(nextToCheck.ptr()); // removed const Type pointer
        nextToCheck.reset(auxToCheck.ptr()); // issue set -> reset compiler warning
        auxToCheck = std::move(helpPtr); // added std::move

        // check if all processors finished 
        nextSize = nextToCheck().size();
        reduce(nextSize, maxOp<label>());
    }

    // clear cached neighbours of cells that have not center inside body
    DynamicLabelList keyToErase;
    for(auto it = cachedNeighbours_().begin(); it != cachedNeighbours_().end(); ++it)
    {
        if(!cellInside.found(it.key()))
        {
            keyToErase.append(it.key());
        }
    }
    cachedNeighbours_().erase(keyToErase);

    // clear processor neighbours of cells that have not center inside body
    keyToErase.clear();
    for(auto it = procNeighbours_().begin(); it != procNeighbours_().end(); ++it)
    {
        if(!cellInside.found(it.key()))
        {
            keyToErase.append(it.key());
        }
    }
    procNeighbours_().erase(keyToErase);

    DynamicLabelList potentSurfCells = 
        getPotentSurfCells(
            body,
            cellInside,
            cellPoints
        );

    correctSurfCells
    (
        body,
        potentSurfCells,
        cellInside,
        cellPoints
    );

    if(intCells_[Pstream::myProcNo()].size() > 0)
    {
        cellToStartInCreateIB_ = min(intCells_[Pstream::myProcNo()]);
    }

    // look for halo cells
    findHaloCells(body);
}
//---------------------------------------------------------------------------//
// Find first cell with center inside the body
label convexBody::getCellInBody
(
    Field<label>& octreeField
)
{
    // octreeField *= 0;
    labelHashSet checkedCells;
    // get the list of cell centroids
    const pointField& cp = mesh_.C();

    if(cellToStartInCreateIB_ >= octreeField.size())
        cellToStartInCreateIB_ = 0;

    autoPtr<DynamicLabelList> nextToCheck(
        new DynamicLabelList(1,cellToStartInCreateIB_));
    autoPtr<DynamicLabelList> auxToCheck(
        new DynamicLabelList);

    label iterCount(0);label iterMax(mesh_.nCells());

    while (nextToCheck().size() > 0 and iterCount < iterMax)
    {
        auxToCheck().clear();
        forAll (nextToCheck(),cellToCheck)
        {
            if (!checkedCells.found(nextToCheck()[cellToCheck]))
            {
                checkedCells.insert(nextToCheck()[cellToCheck]);
                iterCount++;

                if(pointInside(cp[nextToCheck()[cellToCheck]]))
                {
                    return nextToCheck()[cellToCheck];
                }
                else
                {
                    auxToCheck().append(mesh_.cellCells()[nextToCheck()[cellToCheck]]);
                }
            }
        }
        autoPtr<DynamicLabelList> helpPtr(nextToCheck.ptr());
        nextToCheck.reset(auxToCheck.ptr());
        auxToCheck = std::move(helpPtr);
    }
    return -1;
}
//---------------------------------------------------------------------------//
