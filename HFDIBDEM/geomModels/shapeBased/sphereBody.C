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
    label nGeomDir(0);

    forAll(ibBoundPoints,point)
    {
        bool dirOk(true);
        forAll(geometricD_,dir)
        {
            if(geometricD_[dir] == 1)
            {
                if(!(curMeshBounds_.min()[dir] < ibBoundPoints[point][dir] && curMeshBounds_.max()[dir] > ibBoundPoints[point][dir]))
                {
                    dirOk = false;
                }
                nGeomDir += 1;
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

                        if(nGeomDir == 3)
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
    List<pointField>& cellPoints
)
{
    boundBox ibBound(getBounds());
    bool ibInsideMesh(false);

    pointField ibBoundPoints(ibBound.points());

    forAll(ibBoundPoints,point)
    {
        bool dirOk(true);
        forAll(geometricD_,dir)
        {
            if(geometricD_[dir] == 1)
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

    // clear old list contents
    intCells[Pstream::myProcNo()].clear();
    surfCells[Pstream::myProcNo()].clear();
    // find the processor with most of this IB inside
    ibPartialVolume_[Pstream::myProcNo()] = 0;
    octreeField *= 0;

    if(ibInsideMesh)
    {
        // get the list of cell centroids
        const pointField& cp = mesh_.C();

        bool insideIB(false);
        bool insideIbBound(false);

        if(cellToStartInCreateIB_ >= octreeField.size())
            cellToStartInCreateIB_ = 0;

        labelList nextToCheck(1,cellToStartInCreateIB_);
        label iterCount(0);label iterMax(mesh_.nCells());
        labelList vertexLabels;
        boolList vertexesInside;
        pointField pointPos;
        bool centerInside;
        DynamicLabelList auxToCheck;
        HashTable<pointField,label,Hash<label>> lastIbPoints(0);
        while (nextToCheck.size() > 0 and iterCount < iterMax)
        {
            iterCount++;
            auxToCheck.clear();

            forAll (nextToCheck,cellToCheck)
            {
                if (octreeField[nextToCheck[cellToCheck]] == 0)
                {
                    octreeField[nextToCheck[cellToCheck]] = 1;

                    if(lastIbPoints_.found(nextToCheck[cellToCheck]))
                    {
                        pointPos = lastIbPoints_[nextToCheck[cellToCheck]];
                    }
                    else
                    {
                        pointPos = cellPoints[nextToCheck[cellToCheck]];
                    }
                    bool cellInsideBB(false);
                    forAll(pointPos,pos)
                    {
                        if(ibBound.contains(pointPos[pos]))
                        {
                            cellInsideBB = true;
                            break;
                        }
                    }

                    if(cellInsideBB)
                    {
                        insideIbBound = true;
                        cellToStartInCreateIB_ = nextToCheck[cellToCheck];
                        lastIbPoints.insert(nextToCheck[cellToCheck],pointPos);
                        vertexesInside = pointInside(pointPos);
                        pointField centerPoint(1,cp[nextToCheck[cellToCheck]]);
                        centerInside = (pointInside(centerPoint))[0];
                        if(std::any_of(vertexesInside.begin(),vertexesInside.end(),[](bool b){return b;}) || centerInside || !insideIB)
                        {
                            auxToCheck.append(
                                createImmersedBodyByOctTree(
                                    nextToCheck[cellToCheck],
                                    insideIB,
                                    centerInside,
                                    vertexesInside,
                                    body,
                                    surfCells,
                                    intCells
                                )
                            );
                        }
                    }
                    else if(!insideIB && !insideIbBound)
                    {
                        auxToCheck.append(mesh_.cellCells()[nextToCheck[cellToCheck]]);
                    }
                }
            }
            nextToCheck = auxToCheck;
        }
        if(intCells[Pstream::myProcNo()].size() > 0)
            cellToStartInCreateIB_ = min(intCells[Pstream::myProcNo()]);
        lastIbPoints_ = lastIbPoints;
    }
}
//---------------------------------------------------------------------------//
// create immersed body info
labelList sphereBody::createImmersedBodyByOctTree
(
    label cellToCheck,
    bool& insideIB,
    bool& centerInside,
    boolList& vertexesInside,
    volScalarField& body,
    List<DynamicLabelList>& surfCells,
    List<DynamicLabelList>& intCells
)
{
    labelList retList;
    scalar rVInSize(0.5/vertexesInside.size());
    // Note: weight of a single vertex in the cell

    scalar cBody(0);
    forAll (vertexesInside, verIn)
    {
        if (vertexesInside[verIn]==true)
        {
            cBody  += rVInSize;                                         //fraction of cell covered
        }
    }

    vector sDSpan(4.0*(mesh_.bounds().max()-mesh_.bounds().min()));
    if (centerInside)
    {
        cBody+=0.5;
    }
    bool cellInside(false);
    if (cBody > thrSurf_)
    {
        if (cBody > (1.0-thrSurf_))
        {
            intCells[Pstream::myProcNo()].append(cellToCheck);
        }
        else if (cBody  <= (1.0-thrSurf_))
        {
            surfCells[Pstream::myProcNo()].append(cellToCheck);
            if (sdBasedLambda_)
            {
                vector nearestDir(mesh_.C()[cellToCheck]-position_);
                scalar signedDist(0.0);
                if (mag(nearestDir) > SMALL )
                {
                    nearestDir /= mag(nearestDir);
                    vector nearest(position_ + nearestDir*radius_);
                    signedDist = mag(nearest-mesh_.C()[cellToCheck]);
                }
                else
                {
                    Info << "Missed point in signedDist computation !!" << endl;
                }
                if (centerInside)
                {
                    cBody = 0.5*(Foam::tanh(intSpan_*signedDist/Foam::pow(mesh_.V()[cellToCheck],0.333))+1.0);
                }
                else
                {
                    cBody = 0.5*(-1.0*Foam::tanh(intSpan_*signedDist/Foam::pow(mesh_.V()[cellToCheck],0.333))+1.0);
                }
            }
        }
        ibPartialVolume_[Pstream::myProcNo()] += 1;
        cellInside = true;
        insideIB = true;
    }
    body[cellToCheck]+= cBody;
    
    // clip the body field values
    body[cellToCheck] = min(max(0.0,body[cellToCheck]),1.0);
    if (!insideIB || cellInside)
    {
        retList = mesh_.cellCells()[cellToCheck];
    }

    return retList;
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
