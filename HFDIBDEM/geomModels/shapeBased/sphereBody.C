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
#include "sphereBody.H"

using namespace Foam;

//---------------------------------------------------------------------------//
bool sphereBody::canAddBody
(
    const volScalarField& body
)
{
    boundBox ibBound(getBounds());
    bool ibInsideMesh(false);
    pointField ibBoundPoints(ibBound.points());

    forAll(ibBoundPoints,point)
    {
        bool dirOk(true);
        forAll(geometricD,dir)
        {
            if(geometricD[dir] == 1)
            {
                if(!(curMeshBounds_.min()[dir] < ibBoundPoints[point][dir] && curMeshBounds_.max()[dir] > ibBoundPoints[point][dir]))
                {
                    dirOk = false;
                }
            }
        }

        if(dirOk)
        {
            ibInsideMesh = true;
            break;
        }
    }

    if(!ibInsideMesh)
        return true;

    Field<label> octreeField(mesh_.nCells(),0);

    const pointField& pp = mesh_.points();

    labelList nextToCheck(1,0);
    label iterCount(0);label iterMax(mesh_.nCells());
    labelList vertexLabels;
    boolList vertexesInside;
    pointField pointPos;
    bool insideIB(false);
    DynamicLabelList auxToCheck;
    while (nextToCheck.size() > 0 and iterCount < iterMax)
    {
        iterCount++;
        auxToCheck.clear();

        forAll (nextToCheck,cellToCheck)
        {
            if (octreeField[nextToCheck[cellToCheck]] == 0)
            {
                octreeField[nextToCheck[cellToCheck]] = 1;

                vertexLabels = mesh_.cellPoints()[nextToCheck[cellToCheck]];
                pointPos = filterField(pp,vertexLabels);
                vertexesInside = pointInside(pointPos);
                bool cellInside(false);
                forAll (vertexesInside, verIn)
                {
                    if (vertexesInside[verIn]==true)
                    {
                        cellInside = true;
                        insideIB = true;
                        if(body[nextToCheck[cellToCheck]] > SMALL)
                        {
                            return false;
                        }

                        if(case3D)
                        {
                            const labelList& cFaces = mesh_.cells()[nextToCheck[cellToCheck]];

                            forAll (cFaces,faceI)
                            {
                                if (!mesh_.isInternalFace(cFaces[faceI]))
                                {
                                    // get reference to the patch which is in contact with IB. There is contact only if the patch is marked as a wall
                                    label facePatchId(-1);
                                    facePatchId = mesh_.boundaryMesh().whichPatch(cFaces[faceI]);
                                    const polyPatch& cPatch = mesh_.boundaryMesh()[facePatchId];
                                    if (cPatch.type()=="wall" || cPatch.type()=="patch")
                                    {
                                        pointField points = mesh_.faces()[cFaces[faceI]].points(pp);
                                        boolList faceVertexesInside = pointInside(pointPos);
                                        forAll (faceVertexesInside, verIn)
                                        {
                                            if (faceVertexesInside[verIn]==true)
                                            {
                                                return false;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                if (!insideIB || cellInside)
                {
                    auxToCheck.append(mesh_.cellCells()[nextToCheck[cellToCheck]]);
                }
            }
        }
        nextToCheck = auxToCheck;
    }

    return true;
}
//---------------------------------------------------------------------------//
// create immersed body for convex body
void sphereBody::createImmersedBody
(
    volScalarField& body,
    Field<label>& octreeField,
    List<DynamicLabelList>& surfCells,
    List<DynamicLabelList>& intCells,
    List<labelList>& cellPoints
)
{
    clockTime stopWatch;
    stopWatch.timeIncrement();
    // clear old list contents
    intCells[Pstream::myProcNo()].clear();
    surfCells[Pstream::myProcNo()].clear();
    // find the processor with most of this IB inside
    ibPartialVolume_[Pstream::myProcNo()] = 0;

    if(!isBBoxInMesh())
    {
        return;
    }

    label cellInIB = getCellInBody(octreeField);
    if(cellInIB == -1)
    {
        return;
    }

    // get the list of cell centroids
    const pointField& cp = mesh_.C();
    // octreeField *= 0;

    autoPtr<DynamicLabelList> nextToCheck(
        new DynamicLabelList(1,cellInIB));
    autoPtr<DynamicLabelList> auxToCheck(
        new DynamicLabelList);

    label tableSize = 128;
    if(cachedNeighbours_.valid())
    {
        tableSize = cachedNeighbours_().toc().size()*1.5;
    }
    else
    {
        cachedNeighbours_ = new HashTable<const labelList&, label, Hash<label>>;
    }

    HashTable<bool, label, Hash<label>> cellInside(tableSize);
    // autoPtr<HashTable<const labelList&, label, Hash<label>>> newNeigh( 
    //     new HashTable<const labelList&, label, Hash<label>>(tableSize));

    label iterCount(0);label iterMax(mesh_.nCells());
    while (nextToCheck().size() > 0 and iterCount++ < iterMax)
    {
        auxToCheck().clear();
        forAll (nextToCheck(),cellToCheck)
        {
            label cCell = nextToCheck()[cellToCheck];
            if (!cellInside.found(cCell))
            {
                iterCount++;

                if(pointInside(cp[cCell]))
                {
                    cellInside.set(cCell, true);

                    if(cachedNeighbours_.valid() && cachedNeighbours_().found(cCell))
                    {
                        // const labelList& neigh = cachedNeighbours_()[cCell];
                        // newNeigh().insert(cCell, neigh);
                        auxToCheck().append(cachedNeighbours_()[cCell]);
                    }
                    else
                    {
                        const labelList& neigh = mesh_.cellCells(cCell);
                        // newNeigh().insert(cCell, neigh);
                        cachedNeighbours_().insert(cCell, neigh);
                        auxToCheck().append(neigh);
                    }
                }
                else
                {
                    cellInside.set(cCell, false);
                }
            }
        }
        const autoPtr<DynamicLabelList> helpPtr(nextToCheck.ptr());
        nextToCheck.set(auxToCheck.ptr());
        auxToCheck = helpPtr;
    }

    // cachedNeighbours_ = newNeigh;
    DynamicLabelList keyToErase;
    for(auto it = cachedNeighbours_().begin(); it != cachedNeighbours_().end(); ++it)
    {
        if(!cellInside.found(it.key()))
        {
            keyToErase.append(it.key());
        }
    }
    cachedNeighbours_().erase(keyToErase);

    geomModel::t1 += stopWatch.timeIncrement();

    DynamicLabelList potentSurfCells = 
        getPotentSurfCells(
            body,
            intCells,
            cellInside,
            cellPoints
        );

    geomModel::t2 += stopWatch.timeIncrement();
    correctSurfCells
    (
        body,
        intCells,
        surfCells,
        potentSurfCells,
        cellInside,
        cellPoints
    );

    if(intCells[Pstream::myProcNo()].size() > 0)
    {
        cellToStartInCreateIB_ = min(intCells[Pstream::myProcNo()]);
    }
    geomModel::t3 += stopWatch.timeIncrement();
}
//---------------------------------------------------------------------------//
// Find first cell with center inside the body
label sphereBody::getCellInBody
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
        const autoPtr<DynamicLabelList> helpPtr(nextToCheck.ptr());
        nextToCheck.set(auxToCheck.ptr());
        auxToCheck = helpPtr;
    }
    return -1;
}
//---------------------------------------------------------------------------//
void sphereBody::synchronPos()
{
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    if (owner_ == Pstream::myProcNo())
    {
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            UOPstream send(proci, pBufs);
            send << position_;
        }
    }

    pBufs.finishedSends();
    // move body to points calculated by owner_
    UIPstream recv(owner_, pBufs);
    vector pos (recv);

    // move mesh
    position_ = pos;
}
//---------------------------------------------------------------------------//
boolList sphereBody::pointInside(pointField pointI)
{
    boolList inside(pointI.size());

    forAll(pointI,point)
    {
        inside[point] = mag(position_-pointI[point]) < radius_;
    }

    return inside;
}
//---------------------------------------------------------------------------//
bool sphereBody::pointInside(point pointI)
{
    return mag(position_-pointI) < radius_;
}
//---------------------------------------------------------------------------//
pointField sphereBody::sampleSurfacePoints()
{
    pointField returnField(6);

    vector a(1,0,0);
    vector b(0,1,0);
    vector c(0,0,1);

    List<vector> listV(3);
    listV[0] = a;
    listV[1] = b;
    listV[2] = c;

    forAll (listV,v)
    {
        returnField[v] = position_ + listV[v] * radius_;
    }
    forAll (listV,v)
    {
        returnField[3+v] = position_ - listV[v] * radius_;
    }

    return returnField;
}
