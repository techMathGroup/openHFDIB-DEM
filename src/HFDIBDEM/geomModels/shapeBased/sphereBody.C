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
// create immersed body for convex body
void sphereBody::createImmersedBody
(
    volScalarField& body,
    Field<label>& octreeField,
    List<labelList>& cellPoints
)
{
    // clear old list contents
    intCells_[Pstream::myProcNo()].clear();
    surfCells_[Pstream::myProcNo()].clear();
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
        // cachedNeighbours_ = new HashTable<const labelList&, label, Hash<label>>;
        cachedNeighbours_.reset(new HashTable<labelList, label, Hash<label>>);
    }

    HashTable<bool, label, Hash<label>> cellInside(tableSize);

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
                        auxToCheck().append(cachedNeighbours_()[cCell]);
                    }
                    else
                    {
                        const labelList& neigh = mesh_.cellCells(cCell);
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
        autoPtr<DynamicLabelList> helpPtr(nextToCheck.ptr()); // removed const Type pointer
        nextToCheck.reset(auxToCheck.ptr()); // issue set -> reset compiler warning
        auxToCheck = std::move(helpPtr); // added std::move
    }

    DynamicLabelList keyToErase;
    for(auto it = cachedNeighbours_().begin(); it != cachedNeighbours_().end(); ++it)
    {
        if(!cellInside.found(it.key()))
        {
            keyToErase.append(it.key());
        }
    }
    cachedNeighbours_().erase(keyToErase);

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
        autoPtr<DynamicLabelList> helpPtr(nextToCheck.ptr());
        nextToCheck.reset(auxToCheck.ptr());
        auxToCheck = std::move(helpPtr);
    }
    return -1;
}
//---------------------------------------------------------------------------//
void sphereBody::synchronPos(label owner)
{
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    owner = (owner == -1) ? owner_ : owner;

    if (owner == Pstream::myProcNo())
    {
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            UOPstream send(proci, pBufs);
            send << position_;
        }
    }

    pBufs.finishedSends();
    // move body to points calculated by owner_
    UIPstream recv(owner, pBufs);
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
