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
#include "leastSquaresIntInfo.H"
#include "scalarMatrices.H"

using namespace Foam;

//---------------------------------------------------------------------------//
leastSquaresIntInfo::leastSquaresIntInfo
(
    const  fvMesh&   mesh,
    std::shared_ptr<geomModel>& gModel,
    scalar distFactor,
    scalar radiusFactor,
    scalar angleFactor,
    scalar maxCCRows
)
:
interpolationInfo
(
    mesh,
    gModel
),
distFactor_(distFactor),
radiusFactor_(radiusFactor),
angleFactor_(angleFactor),
maxCCRows_(maxCCRows)
{}
leastSquaresIntInfo::~leastSquaresIntInfo()
{}
//---------------------------------------------------------------------------//
void leastSquaresIntInfo::setIntpInfo()
{
    const DynamicLabelList& cSurfCells = getSurfCells();

    resetIntpInfo(cSurfCells.size());
    List<point>& ibPoints = getIbPoints();
    List<vector>& ibNormals = getIbNormals();
    labelListList& cellCells = getCellCells();

    const vectorField& C = mesh_.cellCentres();
    forAll(cellCells, cellI)
    {
        labelList currCells;
        scalar centerMeanDist;
        findCellCells
        (
            cSurfCells[cellI],
            currCells,
            centerMeanDist
        );
        centerMeanDist *= radiusFactor_;

        vector span(centerMeanDist, centerMeanDist, centerMeanDist);

        geomModel_->getClosestPointAndNormal(
            C[cSurfCells[cellI]],
            span*2,
            ibPoints[cellI],
            ibNormals[cellI]
        );

        scalar angleLimit =
            Foam::cos(angleFactor_*Foam::constant::mathematical::pi/180);

        cellCells[cellI] = labelList(currCells.size(), -1);
        label nUsedCells = 0;
        forAll(currCells, cCellI)
        {
            label currCell = currCells[cCellI];
            scalar r = mag(C[currCell] - C[cSurfCells[cellI]]);

            if(r <= centerMeanDist)
            {
                vector dir = C[currCell] - ibPoints[cellI];
                if(mag(dir) > 0)
                    dir /= mag(dir);

                if(mag(ibNormals[cellI] & dir) >= angleLimit)
                {
                    cellCells[cellI][nUsedCells++] = currCell;
                }
            }
        }

        cellCells[cellI].setSize(nUsedCells);
    }

    getIbCellsFaces
    (
        cSurfCells
    );

    getInvDirichletMatrix
    (
        cSurfCells
    );
}
//---------------------------------------------------------------------------//
void leastSquaresIntInfo::getInvDirichletMatrix
(
    const DynamicLabelList& ibCells
)
{
    const List<point>& ibPoints = getIbPoints();
    //const List<vector>& ibNormals = getIbNormals();
    const labelListList& cellCells = getCellCells();
    resetInvDirMat(ibCells.size());
    PtrList<scalarRectangularMatrix>& cInvDirMats = getInvDirMats();
    const vectorField& C = mesh_.cellCentres();

    label nCoeffs = 5;
    if(case3D)
    {
        nCoeffs += 4;
    }

    forAll(ibCells, cellI)
    {
        const labelList& curCells = cellCells[cellI];

        vectorField allPoints(curCells.size(), vector::zero);

        if(allPoints.size() < nCoeffs)
        {
            Info << "Not eneough interp Points" << endl;
        }
        else
        {
            forAll(curCells, cI)
            {
                allPoints[cI] = C[curCells[cI]];
            }

            //Weights
            point origin = C[ibCells[cellI]];
            scalarField curDist = mag(allPoints - origin);
            scalarField W = 0.5*
                (
                    1+cos(Foam::constant::mathematical::pi*curDist/
                    (1.1*max(curDist)))
                );

            cInvDirMats.set
            (
                cellI,
                new scalarRectangularMatrix
                (
                    nCoeffs,
                    allPoints.size(),
                    0.0
                )
            );
            scalarRectangularMatrix& curMatrix = cInvDirMats[cellI];
            scalarRectangularMatrix Mi
            (
                allPoints.size(),
                nCoeffs,
                0.0
            );

            origin = ibPoints[cellI];

            for(label i = 0; i < allPoints.size(); i++)
            {
                scalar X = allPoints[i].x() - origin.x();
                scalar Y = allPoints[i].y() - origin.y();
                scalar Z = allPoints[i].z() - origin.z();
                if(case3D)
                {
                    label coeff = 0;
                    Mi[i][coeff++] = X;
                    Mi[i][coeff++] = Y;
                    Mi[i][coeff++] = X*Y;
                    Mi[i][coeff++] = sqr(X);
                    Mi[i][coeff++] = sqr(Y);
                    Mi[i][coeff++] = Z;
                    Mi[i][coeff++] = X*Z;
                    Mi[i][coeff++] = Y*Z;
                    Mi[i][coeff++] = sqr(Z);
                }
                else
                {
                    labelList A (2,0.0);
                    scalarList dists(3,0.0);
                    dists[0] = X;
                    dists[1] = Y;
                    dists[2] = Z;
                    label validDim = 0;
                    for(label dim = 0; dim < 3; dim++)
                    {
                        if(dim != emptyDim)
                        {
                            A[validDim++] = dim;
                        }
                    }
                    label coeff = 0;
                    Mi[i][coeff++] = dists[A[0]];
                    Mi[i][coeff++] = dists[A[1]];
                    Mi[i][coeff++] = dists[A[0]]*dists[A[1]];
                    Mi[i][coeff++] = sqr(dists[A[0]]);
                    Mi[i][coeff++] = sqr(dists[A[1]]);
                }
            }

            for(label i = 0; i < Mi.m(); i++)
            {
                for(label j = 0; j < Mi.n(); j++)
                {
                    Mi[i][j] *= W[i];
                }
            }

            scalarSquareMatrix lsM(nCoeffs, 0.0);
            for(label i = 0; i < lsM.m(); i++)
            {
                for(label j = 0; j < lsM.n(); j++)
                {
                    for(label k = 0; k < Mi.m(); k++)
                    {
                        lsM[i][j] += Mi[k][i]*Mi[k][j];
                    }
                }
            }

            scalarSquareMatrix invLsM = LUinvert(lsM);

            for(label i = 0; i < lsM.m(); i++)
            {
                for(label j = 0; j < Mi.m(); j++)
                {
                    for(label k = 0; k < lsM.m(); k++)
                    {
                        curMatrix[i][j] += invLsM[i][k]*Mi[j][k]*W[j];
                    }
                }
            }
        }
    }

}
//---------------------------------------------------------------------------//
void leastSquaresIntInfo::getIbCellsFaces
(
    const DynamicLabelList& ibCells
)
{
    labelList ibCellIndic(mesh_.nCells(), -1);

    forAll(ibCells, cellI)
    {
        ibCellIndic[ibCells[cellI]] = cellI;
    }

    DynamicLabelList& ibFaces = getIbFaces();
    DynamicList<bool>& ibFaceFlips = getIbFacesFlips();

    const labelUList& owner = mesh_.owner();                            //OF.com: unallocLabelList& -> labelUList&
    const labelUList& neighbour = mesh_.neighbour();

    forAll(neighbour, faceI)
    {
        if
        (
            ibCellIndic[owner[faceI]] == -1
        &&  ibCellIndic[neighbour[faceI]] > -1
        )
        {
            ibFaces.append(faceI);
            ibFaceFlips.append(false);
        }
        else if
        (
            ibCellIndic[owner[faceI]] > -1
        &&  ibCellIndic[neighbour[faceI]] == -1
        )
        {
            ibFaces.append(faceI);
            ibFaceFlips.append(true);
        }
    }
}
//---------------------------------------------------------------------------//
void leastSquaresIntInfo::findCellCells
(
    const label cellId,
    labelList& cellCells,
    scalar& centerMeanDist
)
{
    labelHashSet cellSet;
    cellSet.insert(cellId);
    labelList currCells;
    labelList auxCells = cellSet.toc();

    for(label nRow = 0; nRow < maxCCRows_; nRow++)
    {
        currCells = auxCells;
        auxCells.clear();
        forAll(currCells, cellI)
        {
            label curCell = currCells[cellI];
            const labelList& curCellPoints = mesh_.cellPoints()[curCell];
            forAll(curCellPoints, pointI)
            {
                label curPoint = curCellPoints[pointI];
                const labelList& curPointCells = mesh_.pointCells()[curPoint];
                forAll(curPointCells, nCellI)
                {
                    label nCellId = curPointCells[nCellI];
                    if(!cellSet.found(nCellId))
                    {
                        cellSet.insert(nCellId);
                        auxCells.append(nCellId);
                    }
                }
            }
        }

        if (nRow == 0 && auxCells.size())
        {
            const vectorField& C = mesh_.cellCentres();
            forAll(auxCells, cellI)
            {
                centerMeanDist += mag(C[cellId] - C[auxCells[cellI]]);
            }
            centerMeanDist /= auxCells.size();
        }
    }

    cellSet.erase(cellId);
    cellCells = cellSet.toc();
}
//---------------------------------------------------------------------------//
scalarSquareMatrix leastSquaresIntInfo::LUinvert
(
    scalarSquareMatrix& matrix
)
{
    scalarSquareMatrix luInvert(matrix.m());
    scalarField column(matrix.m());

    labelList pivotIndices(matrix.m());
    LUDecompose(matrix, pivotIndices);

    for(label j = 0; j < matrix.m(); j++)
    {
        for(label i = 0; i < matrix.m(); i++)
        {
            column[i] = 0.0;
        }

        column[j] = 1.0;

        LUBacksubstitute(matrix, pivotIndices, column);

        for(label i = 0; i < matrix.m(); i++)
        {
            luInvert[i][j] = column[i];
        }
    }

    return luInvert;
}
//---------------------------------------------------------------------------//
