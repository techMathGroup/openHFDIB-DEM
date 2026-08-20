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
#include "lineIntInfo.H"

#define ORDER 2

using namespace Foam;

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
    // Info << "!! --> entry check !!" << endl;
    const DynamicLabelList& cSurfCells = getSurfCells();

    resetIntpInfo(cSurfCells.size());
    List<point>& ibPoints = getIbPoints();
    List<vector>& ibNormals = getIbNormals();
    List<List<intPoint>>& intPoints = getIntPoints();
    const scalar& intDist = charCellSize_;                              //use IB characteristic cell size - const in method

    // prepare lists of points to solve
    List<DynamicList<point>> ibPointsToSolve(Pstream::nProcs());        //Note (LK): include this as a part of the intPoint struct?
    List<DynamicList<vector>> ibNormalsToSolve(Pstream::nProcs());      //Note (LK): include this as a part of the intPoint struct?
    List<DynamicList<intPoint>> intPointsToSolve(Pstream::nProcs());

    // create temporary unit surface normals
    forAll (cSurfCells, cellI)
    {
        // get surface cell label
        label scell = cSurfCells[cellI];

        geomModel_->getClosestPointAndNormal(                           //finds surface point and normal
            mesh_.C()[scell],
            intDist*2*vector::one,
            ibPoints[cellI],                                            //stores surface point
            ibNormals[cellI]
        );

        intPoints[cellI].setSize(ORDER+1);
        intPoint cIntPoint
        (
            ibPoints[cellI],                                            //1st is surface point
            scell,                                                      //mesh cell label for surface cell
            Pstream::myProcNo(),                                        //processor with current point
            Pstream::myProcNo(),                                        //processor with originating surface cell
            cellI                                                       //local surface cell label (in cSurfCells list)
        );
        intPoints[cellI][0] = cIntPoint;                                //save as the first interpolation point

        // save for looping lists
        ibPointsToSolve[Pstream::myProcNo()].append(ibPoints[cellI]);
        ibNormalsToSolve[Pstream::myProcNo()].append(ibNormals[cellI]);
        intPointsToSolve[Pstream::myProcNo()].append(cIntPoint);
    }

    // go by orders
    // Info << "!! --> finding int points check !!" << endl;
    // intPointI - counter, after the loop, index of interpolation point
    // .. 0 -> surface point, then interpolation points up to ORDER
    // Note (MI): rethink the comments (to do MI + LK)
    for(label intPointI = 0; intPointI < ORDER; ++intPointI)
    {
        // lists to send
        List<DynamicList<point>> ibPointsToSend(Pstream::nProcs());
        List<DynamicList<vector>> ibNormalsToSend(Pstream::nProcs());
        List<DynamicList<intPoint>> intPointsToSend(Pstream::nProcs());

        // lists to continue
        List<DynamicList<point>> ibPointsToCont(Pstream::nProcs());
        List<DynamicList<vector>> ibNormalsToCont(Pstream::nProcs());
        List<DynamicList<intPoint>> intPointsToCont(Pstream::nProcs());
        

        // loop over processors
        // Note (MI): rethink indexing variables - iInfo vs proci vs cellI
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            // Info << "!!   |-> loop over intPointsToSolve !!" << endl;
            // loop over interpolation points
            forAll(intPointsToSolve[proci], iInfo)
            {
                // latest interpolation point
                intPoint cIntPoint = intPointsToSolve[proci][iInfo];
                point cPoint = cIntPoint.iPoint_;

                do {
                    cPoint += ibNormalsToSolve[proci][iInfo]*intDist;
                } while(pointInCell(cPoint, cIntPoint.iCell_));

                // new interpolation points
                intPoint nIntPoint = findIntPoint(cIntPoint, cPoint);
                correctIntPoint(ibPointsToSolve[proci][iInfo], nIntPoint);

                // check for cells at domain boundary
                if (nIntPoint.iProc_ == -1)
                {
                    // pass to next interpolation point, keeping default iProc_ at -1
                    // iProc_ check by getIntOrder function
                    continue;
                }

                // check if to send or keep
                else if (Pstream::myProcNo() != nIntPoint.iProc_)
                {
                    ibPointsToSend[nIntPoint.iProc_].append(ibPointsToSolve[proci][iInfo]);
                    ibNormalsToSend[nIntPoint.iProc_].append(ibNormalsToSolve[proci][iInfo]);
                    intPointsToSend[nIntPoint.iProc_].append(nIntPoint);
                }

                else
                {
                    ibPointsToCont[proci].append(ibPointsToSolve[proci][iInfo]);
                    ibNormalsToCont[proci].append(ibNormalsToSolve[proci][iInfo]);
                    intPointsToCont[proci].append(nIntPoint);
                }
            }
        }

        // sync with others
        List<DynamicList<point>> ibPointsRecv(Pstream::nProcs());
        List<DynamicList<vector>> ibNormalsRecv(Pstream::nProcs());
        List<DynamicList<intPoint>> intPointsRecv(Pstream::nProcs());
        sendAndRecvIntPoints(
                ibPointsToSend,
                ibNormalsToSend,
                intPointsToSend,
                ibPointsRecv,
                ibNormalsRecv,
                intPointsRecv
            );

        // clear lists
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            ibPointsToSend[proci].clear();
            ibNormalsToSend[proci].clear();
            intPointsToSend[proci].clear();
        }

        // finished solving of received points
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            forAll(intPointsRecv[proci], iInfo)
            {
                // get face label
                label faceI = intPointsRecv[proci][iInfo].iCell_;
                label cellI(0);

                // get cell label
                forAll(mesh_.boundaryMesh(), patchi)
                {
                    const polyPatch& cPatch = mesh_.boundaryMesh()[patchi];
                    if (cPatch.type() == "processor")
                    {
                        const processorPolyPatch& procPatch
                            = refCast<const processorPolyPatch>(cPatch);

                        label sProc = (Pstream::myProcNo() == procPatch.myProcNo())
                            ? procPatch.neighbProcNo() : procPatch.myProcNo();

                        if (sProc == proci)
                        {
                            cellI = mesh_.faceOwner()[cPatch.start() + faceI];
                        }
                    }
                }

                // save
                intPointsRecv[proci][iInfo].iCell_ = cellI;

                // find interpolation point
                intPoint cIntPoint
                (
                    intPointsRecv[proci][iInfo].iPoint_,
                    cellI,
                    Pstream::myProcNo(),
                    intPointsRecv[proci][iInfo].oProc_,
                    intPointsRecv[proci][iInfo].oLabel_
                );

                vector dir = cIntPoint.iPoint_ - ibPointsRecv[proci][iInfo];
                dir /= mag(dir);
                correctIntPoint(ibPointsRecv[proci][iInfo], cIntPoint);

                intPointsRecv[proci][iInfo] = cIntPoint;
            }
        }

        // return solved interpolation points to processor of origin
        List<DynamicList<intPoint>> intPointsSolved(Pstream::nProcs());
        returnSolvedIntPoints(intPointsToCont, intPointsRecv, intPointsSolved);

        // save solved points
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            forAll(intPointsSolved[proci], iInfo)
            {
                // get the original label
                label oLabel = intPointsSolved[proci][iInfo].oLabel_;

                // save int point
                intPoints[oLabel][intPointI+1] = intPointsSolved[proci][iInfo];
            }
        }

        // clear lists and prepare for next order solution
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            ibPointsToSolve[proci].clear();
            ibNormalsToSolve[proci].clear();
            intPointsToSolve[proci].clear();

            forAll(intPointsToCont[proci], iInfo)
            {
                ibPointsToSolve[proci].append(ibPointsToCont[proci][iInfo]);
                ibNormalsToSolve[proci].append(ibNormalsToCont[proci][iInfo]);
                intPointsToSolve[proci].append(intPointsToCont[proci][iInfo]);
            }

            // clear lists
            ibPointsToCont[proci].clear();
            ibNormalsToCont[proci].clear();
            intPointsToCont[proci].clear();
        }

        // add recieved points
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            forAll(intPointsRecv[proci], iInfo)
            {
                ibPointsToSolve[proci].append(ibPointsRecv[proci][iInfo]);
                ibNormalsToSolve[proci].append(ibNormalsRecv[proci][iInfo]);
                intPointsToSolve[proci].append(intPointsRecv[proci][iInfo]);
            }

            // clear lists
            ibPointsRecv[proci].clear();
            ibNormalsRecv[proci].clear();
            intPointsRecv[proci].clear();
        }
    }
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

        forAll (cellFaces, fi)
        {
            const face faceI = mesh_.faces()[cellFaces[fi]];
            vector dir = closestPoint - cPoint.iPoint_;

            pointHit pHit = faceI.ray(
                cPoint.iPoint_,
                dir,
                mesh_.points()
            );

            if(pHit.hit())
            {
                vector newP = 0.95*(pHit.hitPoint() - cPoint.iPoint_);
                newP += cPoint.iPoint_;

                if(pointInCell(newP, cPoint.iCell_))
                {
                    cPoint.iPoint_ = newP;
                    break;
                }
            }
        }
    }
}
//---------------------------------------------------------------------------//
vector lineIntInfo::getClosestPoint
(
    vector ibPoint,
    intPoint& cPoint
)
{
    vector dir = cPoint.iPoint_ - ibPoint;
    dir /= mag(dir);

    vector dirToC = mesh_.C()[cPoint.iCell_] - ibPoint;

    return ibPoint + dir*(dirToC&dir);
}
//---------------------------------------------------------------------------//
intPoint lineIntInfo::findIntPoint
(
    intPoint& fromP,
    point& endP
)
{
    if(fromP.iPoint_ == endP)
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
        label faceInDir = -1;
        while(!pointInCell(retP.iPoint_, retP.iCell_))
        {
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
)
{
    label faceToReturn = -1;
    vector dir = retPoint.iPoint_ - mesh_.C()[retPoint.iCell_];

    const labelList& cellFaces(mesh_.cells()[retPoint.iCell_]);
    scalar dotProd(-GREAT);

    forAll (cellFaces, faceI)
    {
        label fI = cellFaces[faceI];
        if(fI != prevFace)
        {
            vector outNorm = (mesh_.faceOwner()[fI] == retPoint.iCell_)
                ? mesh_.Sf()[fI] : (-1*mesh_.Sf()[fI]);

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
    point pToCheck,
    label cToCheck
)
{
    const labelList& cellFaces(mesh_.cells()[cToCheck]);
    forAll(cellFaces, faceI)
    {
        label fI = cellFaces[faceI];
        vector outNorm = mesh_.Sf()[fI];
        outNorm = (mesh_.faceOwner()[fI] == cToCheck) ? outNorm : (-1*outNorm);

        if (((pToCheck - mesh_.Cf()[fI]) & outNorm) > 0)
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
    List<DynamicList<intPoint>>& intPointsRecv,
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

        forAll(intPointsRecv[proci], iInfo)
        {
            // get the label of the processor of origin
            label oProc = intPointsRecv[proci][iInfo].oProc_;

            // save to solved 
            if (Pstream::myProcNo() == oProc)
            {
                intPointsSolved[oProc].append(intPointsRecv[proci][iInfo]);
            }

            // add to return
            else
            {
                intPointsToRetr[oProc].append(intPointsRecv[proci][iInfo]);
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
