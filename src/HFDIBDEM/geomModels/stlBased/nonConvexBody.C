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

    bool isInsideBB(false);
    labelList nextToCheck(1,0);
    label iterCount(0);label iterMax(mesh_.nCells());
    while ((nextToCheck.size() > 0 or not isInsideBB) && iterCount < iterMax)
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
        }
        nextToCheck = auxToCheck;
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
    forAll (bBoxCells[Pstream::myProcNo()],bCellI)                       //go only through bBox
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
    Pstream::scatter(ibPartialVolume_, 0);

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
