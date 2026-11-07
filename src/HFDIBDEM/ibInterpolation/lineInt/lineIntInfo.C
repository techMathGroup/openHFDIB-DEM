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

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/
#include "lineIntInfo.H"

#include <fstream>
#include <iomanip>

#define ORDER 2

using namespace Foam;

namespace
{
void writeVector(std::ostream& os, const Foam::vector& value)
{
    os << value.x() << ' ' << value.y() << ' ' << value.z();
}
}

//---------------------------------------------------------------------------//
lineIntInfo::lineIntInfo
(
    const  fvMesh&   mesh,
    std::shared_ptr<geomModel>& gModel
)
:
interpolationInfo(mesh, gModel)
{}
lineIntInfo::~lineIntInfo()
{}
//---------------------------------------------------------------------------//
void lineIntInfo::setIntpInfo()
{
    const DynamicLabelList& cSurfCells = getSurfCells();

    resetIntpInfo(cSurfCells.size());
    setIbCellLabels(cSurfCells);
    List<point>& ibPoints = getIbPoints();
    List<vector>& ibNormals = getIbNormals();
    List<List<intPoint>>& intPoints = getIntPoints();
    const scalar& intDist = charCellSize_;                              //use IB characteristic cell size - const in method

    static const bool dumpMarch(getenv("HFDIB_MARCH_DUMP") != nullptr);
    static const bool strictMarch(getenv("HFDIB_LINEINT_STRICT") != nullptr);
    static std::ofstream marchFile;
    if (dumpMarch && !marchFile.is_open())
    {
        marchFile.open
        (
            "marchProbe." + std::to_string(Pstream::myProcNo()) + ".dat"
        );
        marchFile << std::setprecision(17);
    }

    // prepare lists of points to solve
    List<DynamicList<point>> ibPointsToSolve(Pstream::nProcs());   //Note (LK): include this as a part of the intPoint struct?
    List<DynamicList<vector>> ibNormalsToSolve(Pstream::nProcs()); //Note (LK): include this as a part of the intPoint struct?
    List<DynamicList<intPoint>> intPointsToSolve(Pstream::nProcs());
    List<DynamicList<label>> intEntryFacesToSolve(Pstream::nProcs());

    // create temporary unit surface normals
    forAll (cSurfCells, cellI)
    {
        // get surface cell label
        label scell = cSurfCells[cellI];

        geomModel_->getClosestPointAndNormal(                           //finds surface point and normal
            mesh_.C()[scell],
            intDist*2*vector::one,
            ibPoints[cellI],                              //stores surface point
            ibNormals[cellI]
        );

        intPoints[cellI].setSize(ORDER+1);
        intPoint cIntPoint
        (
            ibPoints[cellI],                                     //1st is surface point
            scell,                                                //mesh cell label for surface cell
            Pstream::myProcNo(),                                  //processor with current point
            Pstream::myProcNo(),                                  //processor with originating surface cell
            cellI                                                //local surface cell label (in cSurfCells list)
        );
        intPoints[cellI][0] = cIntPoint;                                //save as the first interpolation point

        // save for looping lists
        ibPointsToSolve[Pstream::myProcNo()].append(ibPoints[cellI]);
        ibNormalsToSolve[Pstream::myProcNo()].append(ibNormals[cellI]);
        intPointsToSolve[Pstream::myProcNo()].append(cIntPoint);
        intEntryFacesToSolve[Pstream::myProcNo()].append(-1);
    }

    // go by orders
    // intPointI - counter, after the loop, index of interpolation point
    // .. 0 -> surface point, then interpolation points up to ORDER
    // Note (MI): rethink the comments (to do MI + LK)
    for(label intPointI = 0; intPointI < ORDER; ++intPointI)
    {
        const label myProc = Pstream::myProcNo();

        // lists to send
        List<DynamicList<point>> ibPointsToSend(Pstream::nProcs());
        List<DynamicList<vector>> ibNormalsToSend(Pstream::nProcs());
        List<DynamicList<intPoint>> intPointsToSend(Pstream::nProcs());

        // lists to continue and invalid points to return
        List<DynamicList<point>> ibPointsToCont(Pstream::nProcs());
        List<DynamicList<vector>> ibNormalsToCont(Pstream::nProcs());
        List<DynamicList<intPoint>> intPointsToCont(Pstream::nProcs());
        List<DynamicList<intPoint>> intPointsInvalid(Pstream::nProcs());

        label nPending(intPointsToSolve[myProc].size());
        label nMarchRounds(0);
        label nUnresolved(0);
        const label maxMarchRounds(Pstream::nProcs() + 1);
        reduce(nPending, sumOp<label>());

        while(nPending > 0)
        {
            if (nMarchRounds == maxMarchRounds)
            {
                nUnresolved += intPointsToSolve[myProc].size();
                forAll(intPointsToSolve[myProc], iInfo)
                {
                    intPoint invalidPoint(intPointsToSolve[myProc][iInfo]);
                    invalidPoint.iProc_ = -1;
                    intPointsInvalid[myProc].append(invalidPoint);
                }
                ibPointsToSolve[myProc].clear();
                ibNormalsToSolve[myProc].clear();
                intPointsToSolve[myProc].clear();
                intEntryFacesToSolve[myProc].clear();
                break;
            }
            ++nMarchRounds;

            forAll(intPointsToSolve[myProc], iInfo)
            {
                intPoint cIntPoint = intPointsToSolve[myProc][iInfo];
                const label entryFace = intEntryFacesToSolve[myProc][iInfo];
                point cPoint(cIntPoint.iPoint_);

                // only a completed previous-order point gets a new target
                if (entryFace < 0)
                {
                    cPoint = advanceIntPoint
                    (
                        cIntPoint,
                        ibNormalsToSolve[myProc][iInfo]
                    );
                }

                bool traversalCapped(false);
                intPoint nIntPoint
                (
                    findIntPoint
                    (
                        cIntPoint,
                        cPoint,
                        entryFace,
                        traversalCapped
                    )
                );
                if (traversalCapped)
                {
                    ++nUnresolved;
                    nIntPoint.iProc_ = -1;
                }

                if (nIntPoint.iProc_ == -1)
                {
                    intPointsInvalid[myProc].append(nIntPoint);
                }
                else if (nIntPoint.iProc_ != myProc)
                {
                    ibPointsToSend[nIntPoint.iProc_].append
                    (
                        ibPointsToSolve[myProc][iInfo]
                    );
                    ibNormalsToSend[nIntPoint.iProc_].append
                    (
                        ibNormalsToSolve[myProc][iInfo]
                    );
                    intPointsToSend[nIntPoint.iProc_].append(nIntPoint);
                }
                else
                {
                    correctIntPoint
                    (
                        ibPointsToSolve[myProc][iInfo],
                        nIntPoint
                    );
                    const bool contained
                    (
                        pointInCell(nIntPoint.iPoint_, nIntPoint.iCell_)
                    );

                    if (dumpMarch && marchFile.is_open())
                    {
                        marchFile << "MARCH_OWNER " << mesh_.time().value()
                            << ' ' << intPointI + 1 << ' ' << nIntPoint.oProc_
                            << ' ' << nIntPoint.oLabel_ << ' ' << myProc
                            << ' ' << nIntPoint.iCell_ << ' ' << contained
                            << std::endl;
                    }

                    if (!contained)
                    {
                        ++nUnresolved;
                        nIntPoint.iProc_ = -1;
                        intPointsInvalid[myProc].append(nIntPoint);
                    }
                    else
                    {
                        ibPointsToCont[myProc].append
                        (
                            ibPointsToSolve[myProc][iInfo]
                        );
                        ibNormalsToCont[myProc].append
                        (
                            ibNormalsToSolve[myProc][iInfo]
                        );
                        intPointsToCont[myProc].append(nIntPoint);
                    }
                }
            }

            ibPointsToSolve[myProc].clear();
            ibNormalsToSolve[myProc].clear();
            intPointsToSolve[myProc].clear();
            intEntryFacesToSolve[myProc].clear();

            // every rank enters every exchange, including with no work
            List<DynamicList<point>> ibPointsRecv(Pstream::nProcs());
            List<DynamicList<vector>> ibNormalsRecv(Pstream::nProcs());
            List<DynamicList<intPoint>> intPointsRecv(Pstream::nProcs());
            sendAndRecvIntPoints
            (
                ibPointsToSend,
                ibNormalsToSend,
                intPointsToSend,
                ibPointsRecv,
                ibNormalsRecv,
                intPointsRecv
            );

            for (label proci = 0; proci < Pstream::nProcs(); proci++)
            {
                ibPointsToSend[proci].clear();
                ibNormalsToSend[proci].clear();
                intPointsToSend[proci].clear();
            }

            for (label proci = 0; proci < Pstream::nProcs(); proci++)
            {
                forAll(intPointsRecv[proci], iInfo)
                {
                    const label patchFace
                    (
                        intPointsRecv[proci][iInfo].iCell_
                    );
                    const label entryFace(getReceivingFace(patchFace, proci));
                    const label cellI
                    (
                        entryFace < 0 ? 0 : mesh_.faceOwner()[entryFace]
                    );

                    intPoint cIntPoint
                    (
                        intPointsRecv[proci][iInfo].iPoint_,
                        cellI,
                        myProc,
                        intPointsRecv[proci][iInfo].oProc_,
                        intPointsRecv[proci][iInfo].oLabel_
                    );

                    if (dumpMarch && marchFile.is_open())
                    {
                        marchFile << "MARCH_RECV " << mesh_.time().value()
                            << ' ' << intPointI + 1 << ' '
                            << cIntPoint.oProc_ << ' ' << cIntPoint.oLabel_
                            << std::endl;
                    }

                    // the final owner will correct this unchanged target
                    ibPointsToSolve[myProc].append
                    (
                        ibPointsRecv[proci][iInfo]
                    );
                    ibNormalsToSolve[myProc].append
                    (
                        ibNormalsRecv[proci][iInfo]
                    );
                    intPointsToSolve[myProc].append(cIntPoint);
                    intEntryFacesToSolve[myProc].append(entryFace);
                }
            }

            nPending = intPointsToSolve[myProc].size();
            reduce(nPending, sumOp<label>());
        }

        reduce(nUnresolved, sumOp<label>());
        if (nUnresolved > 0)
        {
            if (strictMarch)
            {
                FatalErrorInFunction
                    << nUnresolved << " unresolved marched interpolation "
                    << "points remain" << exit(FatalError);
            }

            static bool warnedUnresolved(false);
            if (Pstream::master() && !warnedUnresolved)
            {
                WarningInFunction
                    << "Unresolved marched interpolation points are "
                    << "invalidated and excluded from later orders. "
                    << "Further warnings of this type are suppressed."
                    << endl;
                warnedUnresolved = true;
            }
        }

        // return solved interpolation points to processor of origin
        List<DynamicList<intPoint>> intPointsSolved(Pstream::nProcs());
        returnSolvedIntPoints
        (
            intPointsToCont,
            intPointsInvalid,
            intPointsSolved
        );

#ifdef FULLDEBUG
        label recordErrors(0);
        HashSet<label> returnedLabels;
#endif

        // save solved points
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            forAll(intPointsSolved[proci], iInfo)
            {
                // get the original label
                label oLabel = intPointsSolved[proci][iInfo].oLabel_;

#ifdef FULLDEBUG
                if (!returnedLabels.insert(oLabel))
                {
                    ++recordErrors;
                }
#endif

                // save int point
                intPoints[oLabel][intPointI+1] = intPointsSolved[proci][iInfo];
                if (dumpMarch && marchFile.is_open())
                {
                    const intPoint& saved = intPoints[oLabel][intPointI+1];
                    marchFile << "MARCH_FINAL " << mesh_.time().value()
                        << ' ' << intPointI + 1 << ' ';
                    writeVector(marchFile, mesh_.C()[cSurfCells[oLabel]]);
                    marchFile << ' ' << saved.iProc_ << ' ';
                    writeVector(marchFile, saved.iPoint_);
                    marchFile << std::endl;
                }
            }
        }

#ifdef FULLDEBUG
        if (returnedLabels.size() != cSurfCells.size())
        {
            ++recordErrors;
        }
        reduce(recordErrors, sumOp<label>());
        if (recordErrors > 0)
        {
            FatalErrorInFunction
                << "Duplicate or missing marched interpolation records"
                << exit(FatalError);
        }
#endif

        // the final owner constructs the next order
        forAll(intPointsToCont[myProc], iInfo)
        {
            ibPointsToSolve[myProc].append(ibPointsToCont[myProc][iInfo]);
            ibNormalsToSolve[myProc].append
            (
                ibNormalsToCont[myProc][iInfo]
            );
            intPointsToSolve[myProc].append(intPointsToCont[myProc][iInfo]);
            intEntryFacesToSolve[myProc].append(-1);
        }
    }
}
//---------------------------------------------------------------------------//
point lineIntInfo::advanceIntPoint
(
    const intPoint& fromP,
    const vector& normal
) const
{
    point cPoint(fromP.iPoint_);
    do
    {
        cPoint += normal*charCellSize_;
    } while(pointInCell(cPoint, fromP.iCell_));

    return cPoint;
}
//---------------------------------------------------------------------------//
label lineIntInfo::getReceivingFace
(
    const label faceI,
    const label srcProc
) const
{
    label entryFace(-1);
    forAll(mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& cPatch = mesh_.boundaryMesh()[patchi];
        if (cPatch.type() == "processor")
        {
            const processorPolyPatch& procPatch
                = refCast<const processorPolyPatch>(cPatch);

            label sProc = (Pstream::myProcNo() == procPatch.myProcNo())
                ? procPatch.neighbProcNo() : procPatch.myProcNo();

            if (sProc == srcProc)
            {
                entryFace = cPatch.start() + faceI;
            }
        }
    }

    return entryFace;
}
//---------------------------------------------------------------------------//
void lineIntInfo::correctIntPoint
(
    point ibPoint,
    intPoint& cPoint
)
{
    if(cPoint.iProc_ != Pstream::myProcNo())
    {
        return;
    }

    vector closestPoint = getClosestPoint(ibPoint, cPoint);

    if(pointInCell(closestPoint, cPoint.iCell_))
    {
        cPoint.iPoint_ = closestPoint;
    }
    else
    {
        const labelList& cellFaces(mesh_.cells()[cPoint.iCell_]);
        point clippedPoint(cPoint.iPoint_);
        scalar nearestDistSqr(GREAT);
        bool clipped(false);

        forAll (cellFaces, fi)
        {
            const face faceI = mesh_.faces()[cellFaces[fi]];
            vector dir = closestPoint - cPoint.iPoint_;
            dir /= mag(dir);
            // Note (LK): faceI.ray is defined in file src/OpenFOAM/meshes/meshShapes/face/faceIntersection.C
            // its readme says that the direction vector does not have to be normalized
            // but we can do it, it should not matter and it will be consistent
            // Note (MI): added dir normalization for consistency

            pointHit pHit = faceI.ray(
                cPoint.iPoint_,
                dir,
                mesh_.points()
            );

            if(pHit.hit())
            {
                vector newP = 0.95*(pHit.hitPoint() - cPoint.iPoint_);
                newP += cPoint.iPoint_;
                const scalar distSqr(magSqr(newP - cPoint.iPoint_));

                // face order can differ between decompositions
                if
                (
                    pointInCell(newP, cPoint.iCell_)
                 && distSqr < nearestDistSqr
                )
                {
                    clippedPoint = newP;
                    nearestDistSqr = distSqr;
                    clipped = true;
                }
            }
        }
        if (clipped)
        {
            cPoint.iPoint_ = clippedPoint;
        }
    }
}
//---------------------------------------------------------------------------//
vector lineIntInfo::getClosestPoint
(
    const vector& ibPoint,
    const intPoint& cPoint
) const
{
    vector dir = cPoint.iPoint_ - ibPoint;
    dir /= mag(dir);


    vector dirToC = mesh_.C()[cPoint.iCell_] - ibPoint;

    return ibPoint + dir*(dirToC & dir);
}
//---------------------------------------------------------------------------//
intPoint lineIntInfo::findIntPoint
(
    const intPoint& fromP,
    const point& endP,
    const label prevFace,
    bool& traversalCapped
) const
{
    traversalCapped = false;
    if(prevFace < 0 && fromP.iPoint_ == endP)
    {
        return intPoint();
    }

    intPoint retP
    (
        endP,
        fromP.iCell_,
        fromP.iProc_,
        fromP.oProc_,
        fromP.oLabel_
    );

    if(fromP.iProc_ == Pstream::myProcNo())
    {
        label faceInDir = prevFace;
        label nTraversed(0);
        const label maxTraversed(max(label(1), mesh_.nCells()) + 1);
        while(!pointInCell(retP.iPoint_, retP.iCell_))
        {
            // bound pathological face cycles
            if (nTraversed == maxTraversed)
            {
                traversalCapped = true;
                retP.iProc_ = -1;
                return retP;
            }
            ++nTraversed;

            faceInDir = getFaceInDir(retP, faceInDir);
            if (!mesh_.isInternalFace(faceInDir))
            {
                label facePatchId(mesh_.boundaryMesh().whichPatch(faceInDir));
                const polyPatch& cPatch = mesh_.boundaryMesh()[facePatchId];

                if (cPatch.type() == "processor")
                {
                    const processorPolyPatch& procPatch
                        = refCast<const processorPolyPatch>(cPatch);

                    label sProc = (Pstream::myProcNo() == procPatch.myProcNo())
                        ? procPatch.neighbProcNo() : procPatch.myProcNo();

                    retP.iCell_ = cPatch.whichFace(faceInDir);
                    retP.iProc_ = sProc;

                    return retP;
                }
                else
                {
                    retP.iProc_ = -1;
                    return retP;
                }
            }

            label owner(mesh_.owner()[faceInDir]);
            label neighbour(mesh_.neighbour()[faceInDir]);
            retP.iCell_ = (retP.iCell_ == neighbour) ? owner : neighbour;
        }

        return retP;
    }

    return retP;
}
//---------------------------------------------------------------------------//
label lineIntInfo::getFaceInDir
(
    const intPoint& retPoint,
    const label prevFace
) const
{
    label faceToReturn = -1;
    vector dir = retPoint.iPoint_ - mesh_.C()[retPoint.iCell_];
    dir /= mag(dir);

    const labelList& cellFaces(mesh_.cells()[retPoint.iCell_]);
    scalar dotProd(-GREAT);

    forAll (cellFaces, faceI)
    {
        label fI = cellFaces[faceI];
        if(fI != prevFace)
        {
            vector outNorm(mesh_.faceAreas()[fI]);

            outNorm = (mesh_.faceOwner()[fI] == retPoint.iCell_)
                ? outNorm : -outNorm;
            outNorm /= mag(outNorm);

            scalar auxDotProd(outNorm & dir);

            if (auxDotProd > dotProd)
            {
                dotProd = auxDotProd;
                faceToReturn = fI;
            }
        }
    }

    return faceToReturn;
}
//---------------------------------------------------------------------------//
bool lineIntInfo::pointInCell
(
    const point& pToCheck,
    const label cToCheck
) const
{
    const labelList& cellFaces(mesh_.cells()[cToCheck]);
    forAll(cellFaces, faceI)
    {
        label fI = cellFaces[faceI];
        vector outNorm(mesh_.faceAreas()[fI]);
        vector faceC(mesh_.faceCentres()[fI]);

        outNorm = (mesh_.faceOwner()[fI] == cToCheck) ? outNorm : -outNorm;   
        if (((pToCheck - faceC) & outNorm) > 0)
        {
            return false;
        }        
    }
    return true;
}

//---------------------------------------------------------------------------//
void lineIntInfo::sendAndRecvIntPoints
(
    List<DynamicList<point>>& ibPointsToSend,
    List<DynamicList<vector>>& ibNormalsToSend,
    List<DynamicList<intPoint>>& intPointsToSend,
    List<DynamicList<point>>& ibPointsRecv,
    List<DynamicList<vector>>& ibNormalsRecv,
    List<DynamicList<intPoint>>& intPointsRecv
)
{
    PstreamBuffers pBufsIbPoints(Pstream::commsTypes::nonBlocking);
    PstreamBuffers pBufsIbNormals(Pstream::commsTypes::nonBlocking);
    PstreamBuffers pBufsIntPoints(Pstream::commsTypes::nonBlocking);
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UOPstream sendIbPoints(proci, pBufsIbPoints);
            UOPstream sendIbNormals(proci, pBufsIbNormals);
            UOPstream sendIntPoints(proci, pBufsIntPoints);
            sendIbPoints << ibPointsToSend[proci];
            sendIbNormals << ibNormalsToSend[proci];
            sendIntPoints << intPointsToSend[proci];
        }
    }

    pBufsIbPoints.finishedSends();
    pBufsIbNormals.finishedSends();
    pBufsIntPoints.finishedSends();

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recvIbPoints(proci, pBufsIbPoints);
            UIPstream recvIbNormals(proci, pBufsIbNormals);
            UIPstream recvIntPoints(proci, pBufsIntPoints);
            DynamicList<point> recIbPoints (recvIbPoints);
            DynamicList<point> recIbNormals (recvIbNormals);
            DynamicList<intPoint> recIntPoints (recvIntPoints);
            ibPointsRecv[proci] = recIbPoints;
            ibNormalsRecv[proci] = recIbNormals;
            intPointsRecv[proci] = recIntPoints;
        }
    }

    pBufsIbPoints.clear();
    pBufsIbNormals.clear();
    pBufsIntPoints.clear();
}
//---------------------------------------------------------------------------//
void lineIntInfo::returnSolvedIntPoints
(
    List<DynamicList<intPoint>>& intPointsToCont,
    List<DynamicList<intPoint>>& intPointsInvalid,
    List<DynamicList<intPoint>>& intPointsSolved
)
{
    // prepare list to return
    List<DynamicList<intPoint>> intPointsToRetr(Pstream::nProcs());
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        forAll(intPointsToCont[proci], iInfo)
        {
            // get the label of the processor of origin
            label oProc = intPointsToCont[proci][iInfo].oProc_;

            // save to solved 
            if (Pstream::myProcNo() == oProc)
            {
                intPointsSolved[oProc].append(intPointsToCont[proci][iInfo]);
            }

            // add to return
            else
            {
                intPointsToRetr[oProc].append(intPointsToCont[proci][iInfo]);
            }
        }

        forAll(intPointsInvalid[proci], iInfo)
        {
            // get the label of the processor of origin
            label oProc = intPointsInvalid[proci][iInfo].oProc_;

            // save to solved
            if (Pstream::myProcNo() == oProc)
            {
                intPointsSolved[oProc].append(intPointsInvalid[proci][iInfo]);
            }

            // add to return
            else
            {
                intPointsToRetr[oProc].append(intPointsInvalid[proci][iInfo]);
            }
        }
    }

    // send
    PstreamBuffers pBufsIntPoints(Pstream::commsTypes::nonBlocking);
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UOPstream sendIntPoints(proci, pBufsIntPoints);
            sendIntPoints << intPointsToRetr[proci];
        }
    }

    pBufsIntPoints.finishedSends();

    // recieve
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recvIntPoints(proci, pBufsIntPoints);
            DynamicList<intPoint> recIntPoints (recvIntPoints);
            intPointsSolved[proci] = recIntPoints;
        }
    }

    pBufsIntPoints.clear();
}
//---------------------------------------------------------------------------//
