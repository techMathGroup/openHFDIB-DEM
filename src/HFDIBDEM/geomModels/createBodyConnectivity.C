/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  | |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___| |_____|_|  |_|
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

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/
#include "geomModel.H"

using namespace Foam;

#include "meshSearch.H"

//---------------------------------------------------------------------------//
// vertex-based lambda evaluation with shared-vertex cache
// (extracted from correctSurfCells)
scalar geomModel::evaluateLambda
(
    label cCell,
    List<labelList>& cellPoints,
    HashTable<bool, label, Hash<label>>& verticesStatus,
    bool& centerInside
)
{
    scalar cBody(0);

    centerInside = pointInside(mesh_.C()[cCell]);
    if (centerInside)
    {
        cBody += 0.5;
    }

    const labelList& cVerts = cellPoints[cCell];
    scalar rVInSize(0.5/cVerts.size());
    forAll(cVerts, vertI)
    {
        if(!verticesStatus.found(cVerts[vertI]))
        {
            bool vertexInside = pointInside(mesh_.points()[cVerts[vertI]]);
            verticesStatus.insert
            (
                cVerts[vertI],
                vertexInside
            );

            if(vertexInside)
            {
                cBody += rVInSize;
            }
        }
        else if(verticesStatus[cVerts[vertI]])
        {
            cBody += rVInSize;
        }
    }

    return cBody;
}
//---------------------------------------------------------------------------//
// phase 0: resolve a seed cell for this rank
label geomModel::resolveSeed
(
    const pointField& cp,
    List<labelList>& cellPoints
)
{
    HashTable<bool, label, Hash<label>> verticesStatus(128);
    bool centerInside(false);

    // 1. fast path: cached surf seed still covered by this body
    if (surfSeed_ != -1 && surfSeed_ < mesh_.nCells()
        && evaluateLambda(surfSeed_, cellPoints, verticesStatus, centerInside)
            > SMALL)
    {
        return surfSeed_;
    }

    // 2. cold path: rank-local findCell on the body points.
    meshSearch searchEng(mesh_);

    pointField candidates(1, getCoM());
    {
        const pointField bodyPoints(getBodyPoints());
        const label maxSamples(min(bodyPoints.size(), label(128)));
        if (maxSamples > 0)
        {
            const label stride(max(bodyPoints.size()/maxSamples, label(1)));
            for (label pI = 0; pI < bodyPoints.size(); pI += stride)
            {
                candidates.append(bodyPoints[pI]);
            }
        }
    }

    forAll(candidates, pI)
    {
        label hit(searchEng.findCell(candidates[pI]));
        if (hit != -1
            && evaluateLambda(hit, cellPoints, verticesStatus, centerInside)
                > SMALL)
        {
            return hit;
        }
    }

    return -1;                                                          // -1 on ranks not holding the body -> they idle through
}
//---------------------------------------------------------------------------//
// phase 0: march from a fully-covered seed to a partial (ring) cell
label geomModel::marchToRing(label seed, const pointField& cp)
{
    // march toward the bbox minimum;
    const Vector<label> geomDirs(mesh_.geometricD());
    const vector validDirs
    (
        (geomDirs.x() == -1 ? 0.0 : 1.0),
        (geomDirs.y() == -1 ? 0.0 : 1.0),
        (geomDirs.z() == -1 ? 0.0 : 1.0)
    );
    const vector dir(cmptMultiply(getBounds().min() - cp[seed], validDirs));

    if (mag(dir) < SMALL)
    {
        return -1;
    }

    const vector dirN(dir/mag(dir));
    label cCell(seed);
    label nextCell(-1);

    label iterCount(0);
    const label iterMax(mesh_.nCells());

    while (iterCount++ < iterMax)
    {
        nextCell = walkOneCell(cCell, dirN);

        if (nextCell == -1 || nextCell == cCell)
        {
            // ran into a domain wall before leaving the body
            return -1;
        }

        if (!pointInside(cp[nextCell]))
        {
            // surface crossed between the two centers -> one of them
            // is a partial cell; prefer the inside one as the ring seed
            return cCell;
        }

        cCell = nextCell;
    }

    return -1;
}
//---------------------------------------------------------------------------//
// one cell step along dir (getFaceInDir pattern from lineIntInfo)
label geomModel::walkOneCell(label cCell, const vector& dir)
{
    const labelList& cFaces = mesh_.cells()[cCell];
    label bestFace(-1);
    label bestIntFace(-1);
    scalar bestDot(-2.0);
    scalar bestIntDot(-2.0);

    forAll(cFaces, fI)
    {
        label faceI = cFaces[fI];
        vector fC(mesh_.faceCentres()[faceI]);
        scalar dirDot((fC - mesh_.C()[cCell]) & dir);
        // prefer internal faces
        if (mesh_.isInternalFace(faceI))
        {
            if (dirDot > bestIntDot)
            {
                bestIntDot = dirDot;
                bestIntFace = faceI;
            }
        }
        if (dirDot > bestDot)
        {
            bestDot = dirDot;
            bestFace = faceI;
        }
    }

    if (bestIntFace != -1)
    {
        label nCell(mesh_.owner()[bestIntFace]);
        if (nCell == cCell)
        {
            nCell = mesh_.neighbour()[bestIntFace];
        }
        return nCell;
    }

    if (bestFace == -1)
    {
        return -1;
    }

    // boundary face: no walk-through
    return -1;
}
//---------------------------------------------------------------------------//
// shared processor-face exchange: send per-rank face lists, receive
// the corresponding cell labels adjacent across the processor boundary
void geomModel::exchangeFaceLabels
(
    List<DynamicLabelList>& facesToSend,
    List<DynamicLabelList>& facesRecv
)
{
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UOPstream sendFaces(proci, pBufs);
            sendFaces << facesToSend[proci];
            facesToSend[proci].clear();
        }
    }
    pBufs.finishedSends();

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recvFaces(proci, pBufs);
            DynamicLabelList recFaces(recvFaces);
            facesRecv[proci] = recFaces;
        }
    }
    pBufs.clear();

    // convert the received patch-local faces to local cell labels
    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& cPatch = mesh_.boundaryMesh()[patchI];
        if (cPatch.type() == "processor")
        {
            const processorPolyPatch& procPatch
                = refCast<const processorPolyPatch>(cPatch);

            label sProc = (Pstream::myProcNo() == procPatch.myProcNo())
                ? procPatch.neighbProcNo() : procPatch.myProcNo();

            forAll(facesRecv[sProc], faceI)
            {
                facesRecv[sProc][faceI]
                    = mesh_.faceOwner()[cPatch.start()
                    + facesRecv[sProc][faceI]];
            }
        }
    }
}
//---------------------------------------------------------------------------//
// connectivity-based immersed body creation
// returns false -> caller falls back to legacy
bool geomModel::createImmersedBodyConnectivity
(
    volScalarField& body,
    Field<label>& octreeField,
    List<labelList>& cellPoints
)
{
    // ---- bookkeeping: clear old lists ----
    // NOTE (MI): the old lambda is NOT zeroed here (legacy resetBody leaves it
    // untouched too); on audit failure everything written is wiped before
    // falling back, so the legacy re-run starts from a clean state
    surfCells_[Pstream::myProcNo()].clear();
    intCells_[Pstream::myProcNo()].clear();
    haloCells_[Pstream::myProcNo()].clear();

    octreeField *= 0;

    ibPartialVolume_[Pstream::myProcNo()] = 0;

    const pointField& cp = mesh_.C();
    const vector sDSpan(4.0*(mesh_.bounds().max()-mesh_.bounds().min()));

    // ================= phase 0: seed =================
    label seed = resolveSeed(cp, cellPoints);                           // -1 on ranks without the body

    label ringSeed(-1);
    if (seed != -1)
    {
        HashTable<bool, label, Hash<label>> seedVerts;
        bool seedCenterInside(false);
        scalar seedLambda = evaluateLambda
        (
            seed, cellPoints, seedVerts, seedCenterInside
        );
        if (seedLambda > SMALL && seedLambda < 1.0 - SMALL)
        {
            ringSeed = seed;      // already a partial (ring) cell
        }
        else if (seedLambda >= 1.0 - SMALL)
        {
            ringSeed = marchToRing(seed, cp);
        }
        // seedLambda <= SMALL: stale seed or degenerate overlap
    }

    // ================= phase 1: surf/shell discovery =================
    // no seed on any rank -> cannot start; let the caller fall back
    {
        label anySeed = (ringSeed != -1) ? 1 : 0;
        reduce(anySeed, maxOp<label>());
        if (anySeed == 0)
        {
            resetSeeds();
            InfoH << iB_Info << "connectivity-based creation found no seed, "
                << "falling back to legacy" << endl;
            return false;
        }
    }

    HashTable<bool, label, Hash<label>> verticesStatus(128);

    // frontier split by the class of the enqueuing sender:
    //  - fromSurf:   enqueued by a surf (or seed) sender -> an int-classified
    //                cell is face-adjacent to surf -> MARK_SHELL_INT
    //  - fromShell:  enqueued by a shell-int sender -> MARK_FILLED_INT,
    //                no phase-1 propagation (phase 2 fills from it)
    DynamicLabelList nextFromSurf;
    DynamicLabelList nextFromShell;
    DynamicLabelList fillFrontier;
    List<DynamicLabelList> facesFromSurf(Pstream::nProcs());
    List<DynamicLabelList> facesFromShellInt(Pstream::nProcs());
    List<DynamicLabelList> facesRecvA(Pstream::nProcs());
    List<DynamicLabelList> facesRecvB(Pstream::nProcs());

    if (ringSeed != -1)
    {
        // the seed may itself classify as int (fully covered); treat it as
        // surf-sent so that it still propagates the discovery
        nextFromSurf.append(ringSeed);
    }
    label nextSize = nextFromSurf.size() + nextFromShell.size();
    reduce(nextSize, maxOp<label>());
    label iterCount(0);
    const label iterMax(mesh_.nCells());

    while (nextSize > 0 && iterCount++ < iterMax)
    {
        DynamicLabelList auxFromSurf;
        DynamicLabelList auxFromShell;

        // ---- process the surf-sent frontier (propagating) ----
        forAll(nextFromSurf, cI)
        {
            label cCell = nextFromSurf[cI];
            if (octreeField[cCell] != MARK_NONE)
            {
                continue;
            }

            bool centerInside(false);
            scalar cBody = evaluateLambda
            (
                cCell, cellPoints, verticesStatus, centerInside
            );

            if (cBody <= SMALL)                                         // outside this body
            {
                octreeField[cCell] = MARK_OUTSIDE;
                continue;                                               // no propagation
            }

            bool isSurf = (cBody < 1.0 - SMALL);

            // --- lambda write: replicate correctSurfCells block ---
            if (isSurf && sdBasedLambda_)
            {
                point closestPoint(point::zero);
                vector normal(vector::zero);
                getClosestPointAndNormal
                (
                    mesh_.C()[cCell],
                    sDSpan,
                    closestPoint,
                    normal
                );
                const scalar signedDist = mag(closestPoint - mesh_.C()[cCell]);
                if (centerInside)
                {
                    cBody = 0.5*(Foam::tanh(intSpan_*signedDist/Foam::pow(mesh_.V()[cCell],0.333))+1.0);
                }
                else
                {
                    cBody = 0.5*(-1.0*Foam::tanh(intSpan_*signedDist/Foam::pow(mesh_.V()[cCell],0.333))+1.0);
                }
            }

            if (isSurf)
            {
                octreeField[cCell] = MARK_SURF;
                surfCells_[Pstream::myProcNo()].append(cCell);
            }
            else
            {
                octreeField[cCell] = MARK_SHELL_INT;
                intCells_[Pstream::myProcNo()].append(cCell);
            }
            ibPartialVolume_[Pstream::myProcNo()] += 1;

            if (overlayLambda_ == "add")
            {
                body[cCell] += cBody;
            }
            else if (overlayLambda_ == "max")
            {
                // take maximum of what was there and what cBody wants
                body[cCell] = max(body[cCell], cBody);
            }
            else // not sure what is wanted default behavior
            {
                body[cCell] += cBody;
            }

            // clip the body field values
            body[cCell] = min(max(0.0,body[cCell]),1.0);

            // ---- propagate: enqueue in-rank face-neighbours ----
            const labelList& cFaces = mesh_.cells()[cCell];
            forAll(cFaces, fI)
            {
                label faceI = cFaces[fI];
                if (mesh_.isInternalFace(faceI))
                {
                    label nCell(mesh_.owner()[faceI]);
                    if (nCell == cCell)
                    {
                        nCell = mesh_.neighbour()[faceI];
                    }
                    if (octreeField[nCell] == MARK_NONE)
                    {
                        if (isSurf)
                        {
                            auxFromSurf.append(nCell);
                        }
                        else
                        {
                            auxFromShell.append(nCell);
                        }
                    }
                }
                else
                {
                    label facePatchI(mesh_.boundaryMesh().whichPatch(faceI));
                    const polyPatch& cPatch = mesh_.boundaryMesh()[facePatchI];
                    if (cPatch.type() == "processor")
                    {
                        const processorPolyPatch& procPatch
                            = refCast<const processorPolyPatch>(cPatch);
                        label iProc = (Pstream::myProcNo() == procPatch.myProcNo())
                            ? procPatch.neighbProcNo() : procPatch.myProcNo();
                        label iFace = cPatch.whichFace(faceI);

                        if (isSurf)
                        {
                            facesFromSurf[iProc].append(iFace);
                        }
                        else
                        {
                            facesFromShellInt[iProc].append(iFace);
                        }
                    }
                }
            }
        }

        // ---- process the shell-sent frontier (non-propagating shell layer) ----
        // these are cells enqueued by a shell-int sender:
        // fully covered => beyond the shell layer => filled and phase 2
        // partially covered => the sender was wrong => surf and 
        // re-enter through the surf path
        forAll(nextFromShell, cI)
        {
            label cCell = nextFromShell[cI];
            if (octreeField[cCell] != MARK_NONE)
            {
                continue;
            }

            bool centerInside(false);
            scalar cBody = evaluateLambda
            (
                cCell, cellPoints, verticesStatus, centerInside
            );

            if (cBody <= SMALL)
            {
                octreeField[cCell] = MARK_OUTSIDE;
                continue;
            }

            if (cBody < 1.0 - SMALL)
            {
                octreeField[cCell] = MARK_SURF;
                surfCells_[Pstream::myProcNo()].append(cCell);
                ibPartialVolume_[Pstream::myProcNo()] += 1;

                if (sdBasedLambda_)
                {
                    point closestPoint(point::zero);
                    vector normal(vector::zero);
                    getClosestPointAndNormal
                    (
                        mesh_.C()[cCell],
                        sDSpan,
                        closestPoint,
                        normal
                    );
                    const scalar signedDist
                        = mag(closestPoint - mesh_.C()[cCell]);
                    if (centerInside)
                    {
                        cBody = 0.5*(Foam::tanh(intSpan_*signedDist/Foam::pow(mesh_.V()[cCell],0.333))+1.0);
                    }
                    else
                    {
                        cBody = 0.5*(-1.0*Foam::tanh(intSpan_*signedDist/Foam::pow(mesh_.V()[cCell],0.333))+1.0);
                    }
                }

                if (overlayLambda_ == "add")
                {
                    body[cCell] += cBody;
                }
                else if (overlayLambda_ == "max")
                {
                    body[cCell] = max(body[cCell], cBody);
                }
                else
                {
                    body[cCell] += cBody;
                }
                body[cCell] = min(max(0.0,body[cCell]),1.0);

                // re-enter through the propagating surf path
                const labelList& cFaces = mesh_.cells()[cCell];
                forAll(cFaces, fI)
                {
                    label faceI = cFaces[fI];
                    if (mesh_.isInternalFace(faceI))
                    {
                        label nCell(mesh_.owner()[faceI]);
                        if (nCell == cCell)
                        {
                            nCell = mesh_.neighbour()[faceI];
                        }
                        if (octreeField[nCell] == MARK_NONE)
                        {
                            auxFromSurf.append(nCell);
                        }
                    }
                    else
                    {
                        label facePatchI(mesh_.boundaryMesh().whichPatch(faceI));
                        const polyPatch& cPatch = mesh_.boundaryMesh()[facePatchI];
                        if (cPatch.type() == "processor")
                        {
                            const processorPolyPatch& procPatch
                                = refCast<const processorPolyPatch>(cPatch);
                            label iProc = (Pstream::myProcNo() == procPatch.myProcNo())
                                ? procPatch.neighbProcNo() : procPatch.myProcNo();
                            label iFace = cPatch.whichFace(faceI);
                            facesFromSurf[iProc].append(iFace);
                        }
                    }
                }
                continue;
            }

            octreeField[cCell] = MARK_FILLED_INT;
            intCells_[Pstream::myProcNo()].append(cCell);
            ibPartialVolume_[Pstream::myProcNo()] += 1;
            body[cCell] = 1.0;
        }

        // ---- cross-rank exchange  ----
        // Note (MI): two lists: receiver knows if sender was surf or 
        //           shell-int
        exchangeFaceLabels(facesFromSurf, facesRecvA);
        exchangeFaceLabels(facesFromShellInt, facesRecvB);

        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            if (proci == Pstream::myProcNo())
            {
                continue;
            }
            forAll(facesRecvA[proci], rC)                               // from surf senders: propagating frontier continuation
            {
                label rCell = facesRecvA[proci][rC];
                if (octreeField[rCell] == MARK_NONE)
                {
                    auxFromSurf.append(rCell);
                }
            }
            forAll(facesRecvB[proci], rC)                               // from shell-int senders: non-propagating frontier
            {
                label rCell = facesRecvB[proci][rC];
                if (octreeField[rCell] == MARK_NONE)
                {
                    auxFromShell.append(rCell);
                }
            }
        }

        // clean up and prep for next iter
        nextFromSurf = auxFromSurf;
        nextFromShell = auxFromShell;
        nextSize = nextFromSurf.size() + nextFromShell.size();
        reduce(nextSize, maxOp<label>());
    }

    // ================= phase 2: interior fill =================
    // pure connectivity from the filled-int layer; ZERO geometry tests.
    // spread through face-neighbors until no more unclassified cells 
    // are reachable.
    fillFrontier.clear();
    forAll(intCells_[Pstream::myProcNo()], cI)
    {
        label cCell = intCells_[Pstream::myProcNo()][cI];
        if (octreeField[cCell] == MARK_FILLED_INT)
        {
            fillFrontier.append(cCell);
        }
    }
    nextSize = fillFrontier.size();
    reduce(nextSize, maxOp<label>());
    iterCount = 0;
    while (nextSize > 0 && iterCount++ < iterMax)
    {
        DynamicLabelList auxFill;

        forAll(fillFrontier, cI)
        {
            label cCell = fillFrontier[cI];
            const labelList& cFaces = mesh_.cells()[cCell];
            forAll(cFaces, fI)
            {
                label faceI = cFaces[fI];
                if (mesh_.isInternalFace(faceI))
                {
                    label nCell(mesh_.owner()[faceI]);
                    if (nCell == cCell)
                    {
                        nCell = mesh_.neighbour()[faceI];
                    }
                    if (octreeField[nCell] == MARK_NONE)
                    {
                        octreeField[nCell] = MARK_FILLED_INT;
                        body[nCell] = 1.0;
                        intCells_[Pstream::myProcNo()].append(nCell);
                        ibPartialVolume_[Pstream::myProcNo()] += 1;
                        auxFill.append(nCell);
                    }
                }
                else
                {
                    label facePatchI(mesh_.boundaryMesh().whichPatch(faceI));
                    const polyPatch& cPatch = mesh_.boundaryMesh()[facePatchI];
                    if (cPatch.type() == "processor")
                    {
                        const processorPolyPatch& procPatch
                            = refCast<const processorPolyPatch>(cPatch);
                        label iProc = (Pstream::myProcNo() == procPatch.myProcNo())
                            ? procPatch.neighbProcNo() : procPatch.myProcNo();
                        label iFace = cPatch.whichFace(faceI);
                        facesFromSurf[iProc].append(iFace);             // plain fill spread
                    }
                }
            }
        }

        // cross-rank fill spread (sender class irrelevant, reuse one list)
        exchangeFaceLabels(facesFromSurf, facesRecvA);
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            if (proci == Pstream::myProcNo())
            {
                continue;
            }
            forAll(facesRecvA[proci], rC)
            {
                label rCell = facesRecvA[proci][rC];
                if (octreeField[rCell] == MARK_NONE)
                {
                    octreeField[rCell] = MARK_FILLED_INT;
                    body[rCell] = 1.0;
                    intCells_[Pstream::myProcNo()].append(rCell);
                    ibPartialVolume_[Pstream::myProcNo()] += 1;
                    auxFill.append(rCell);
                }
            }
        }

        fillFrontier = auxFill;
        nextSize = fillFrontier.size();
        reduce(nextSize, maxOp<label>());
    }

    // ================= phase 3: audit =================
    bool anomaly(false);

    // (a) filled-int cell adjacent to a classified-outside cell:
    //     impossible for a true int cell (density theorem) -> fill leaked
    {
        const DynamicLabelList& intLst = intCells_[Pstream::myProcNo()];
        forAll(intLst, cI)
        {
            label cCell = intLst[cI];
            if (octreeField[cCell] != MARK_FILLED_INT)
            {
                continue;
            }
            const labelList& cFaces = mesh_.cells()[cCell];
            forAll(cFaces, fI)
            {
                label faceI = cFaces[fI];
                if (mesh_.isInternalFace(faceI))
                {
                    label nCell(mesh_.owner()[faceI]);
                    if (nCell == cCell)
                    {
                        nCell = mesh_.neighbour()[faceI];
                    }
                    if (octreeField[nCell] == MARK_OUTSIDE)
                    {
                        anomaly = true;
                        break;
                    }
                }
            }
            if (anomaly)
            {
                break;
            }
        }
    }

    // (b) cell-count sanity vs the last successful creation
    //     (loose Co-driven growth bound; no baseline on the first step)
    {
        label nCellsLoc
            = surfCells_[Pstream::myProcNo()].size()
            + intCells_[Pstream::myProcNo()].size();
        label nCellsGlob(nCellsLoc);
        reduce(nCellsGlob, sumOp<label>());

        // (b0) a body must cover at least one cell; an empty result would
        // otherwise pass the audit and crash downstream on M0 = 0
        if (nCellsGlob == 0)
        {
            anomaly = true;
        }

        if (nCellsPrev_ > 0)
        {
            // a body can only grow by the boundary flux per step; allow a
            // loose factor of 4 above pure growth-proportional bound
            label growBound = nCellsPrev_ + max(4*nCellsPrev_/100, 10);
            if (nCellsGlob > growBound)
            {
                anomaly = true;
            }
        }

        // (c) surf list non-empty on at least one rank
        label surfLoc = surfCells_[Pstream::myProcNo()].size();
        reduce(surfLoc, sumOp<label>());
        if (surfLoc == 0 && nCellsGlob > 0)
        {
            anomaly = true;
        }

        nCellsPrev_ = nCellsGlob;
    }

    reduce(anomaly, orOp<bool>());
    if (anomaly)
    {
        // wipe written to start the legacy re-run clean
        DynamicLabelList writtenCells(surfCells_[Pstream::myProcNo()]);
        writtenCells.append(intCells_[Pstream::myProcNo()]);
        forAll(writtenCells, c)
        {
            body[writtenCells[c]] = 0;
        }
        surfCells_[Pstream::myProcNo()].clear();
        intCells_[Pstream::myProcNo()].clear();
        haloCells_[Pstream::myProcNo()].clear();
        resetSeeds();
        InfoH << iB_Info << "connectivity-based creation failed audit, "
            << "falling back to legacy" << endl;
        return false;
    }

    // ---- success: update caches for next step ----
    if (surfCells_[Pstream::myProcNo()].size() > 0)
    {
        surfSeed_ = surfCells_[Pstream::myProcNo()][0];
    }
    if (intCells_[Pstream::myProcNo()].size() > 0)
    {
        cellToStartInCreateIB_ = intCells_[Pstream::myProcNo()][0];
    }
    else if (surfCells_[Pstream::myProcNo()].size() > 0)
    {
        cellToStartInCreateIB_ = surfSeed_;
    }

    // ---- owner determination (as in nonConvexBody legacy) ----
    Pstream::gatherList(ibPartialVolume_, 0);
    Pstream::broadcast(ibPartialVolume_, 0);
    for (label i = 0; i < ibPartialVolume_.size(); i++)
    {
        if (ibPartialVolume_[i] == max(ibPartialVolume_))
        {
            owner_ = i;
            break;
        }
    }

    // ---- halo cells: unchanged post-pass ----
    findHaloCells(body);

    return true;
}
//---------------------------------------------------------------------------//
