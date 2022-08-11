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
#include "geomModel.H"

using namespace Foam;

//---------------------------------------------------------------------------//
geomModel::geomModel
(
    const  fvMesh&   mesh,
    const contactType cType,
    scalar  thrSurf
)
:
contactType_(cType),
mesh_(mesh),
owner_(0),
cellToStartInCreateIB_(0),
thrSurf_(thrSurf),
intSpan_(2.0),
sdBasedLambda_(false),
curMeshBounds_(mesh_.points(),false),
M_(0.0),
M0_(0.0),
CoM_(vector::zero),
I_(symmTensor::zero),
dC_(0.0),
rhoS_("rho",dimensionSet(1,-3,0,0,0,0,0),1.0)
{
    surfCells_.setSize(Pstream::nProcs());
    intCells_.setSize(Pstream::nProcs());
    ibPartialVolume_.setSize(Pstream::nProcs());
}
geomModel::~geomModel()
{
}
//---------------------------------------------------------------------------//
void geomModel::calculateGeometricalProperties
(
    volScalarField& body
)
{
    //vector CoMOld  = CoM_;
    M_      = scalar(0);
    //CoM_    = vector::zero;
    setCoM();
    I_      = symmTensor::zero;
    //vector tmpCom(vector::zero);

    addToMAndI(body,surfCells_[Pstream::myProcNo()]);
    addToMAndI(body,intCells_[Pstream::myProcNo()]);

    // collect from processors
    reduce(M_, sumOp<scalar>());
    //reduce(tmpCom,  sumOp<vector>());
    reduce(I_,  sumOp<symmTensor>());

    /*if(M_ > 0)
    {
        CoM_ = tmpCom / M_;
    }*/
}
//---------------------------------------------------------------------------//
void geomModel::addToMAndI
(
    volScalarField& body,
    DynamicLabelList& labelCellLst
)
{
    forAll (labelCellLst,cell)
    {
        label cellI  = labelCellLst[cell];

        scalar Mi    = body[cellI]*rhoS_.value()*mesh_.V()[cellI];
        // add to M_

        M_      += Mi;
        scalar xLoc = mesh_.C()[cellI].x() - CoM_.x();
        scalar yLoc = mesh_.C()[cellI].y() - CoM_.y();
        scalar zLoc = mesh_.C()[cellI].z() - CoM_.z();

        // add to I_
        I_.xx() += Mi*(yLoc*yLoc + zLoc*zLoc);
        I_.yy() += Mi*(xLoc*xLoc + zLoc*zLoc);
        I_.zz() += Mi*(xLoc*xLoc + yLoc*yLoc);

        I_.xy() -= Mi*(xLoc*yLoc);
        I_.xz() -= Mi*(xLoc*zLoc);
        I_.yz() -= Mi*(yLoc*zLoc);
    }
}
//---------------------------------------------------------------------------//
void geomModel::computeBodyCharPars()
{
    forAll (surfCells_[Pstream::myProcNo()],sCellI)
    {
        label cellI = surfCells_[Pstream::myProcNo()][sCellI];
        dC_ = max(dC_,mag(CoM_-mesh_.C()[cellI]));
    }
    M0_ = M_;
    reduce(dC_, maxOp<scalar>());
}
//---------------------------------------------------------------------------//
void geomModel::resetBody(volScalarField& body)
{
    forAll (intCells_[Pstream::myProcNo()],cellI)
    {
        body[intCells_[Pstream::myProcNo()][cellI]] = 0;
    }
    forAll (surfCells_[Pstream::myProcNo()],cellI)
    {
        body[surfCells_[Pstream::myProcNo()][cellI]] = 0;
    }

    surfCells_[Pstream::myProcNo()].clear();
    intCells_[Pstream::myProcNo()].clear();
    }
//---------------------------------------------------------------------------//
bool geomModel::isBBoxInMesh()
{
    boundBox ibBound(getBounds());

    forAll(geometricD,dir)
    {
        if(geometricD[dir] == 1)
        {
            if(!(curMeshBounds_.max()[dir] >= ibBound.min()[dir]
                && curMeshBounds_.min()[dir] <= ibBound.max()[dir]))
            {
                return false;
            }
        }
    }
    return true;
}
//---------------------------------------------------------------------------//
DynamicList<label> geomModel::getPotentSurfCells
(
    volScalarField& body,
    HashTable<bool, label, Hash<label>>& cellInside,
    List<labelList>& cellPoints
)
{
    const labelList foundCells = cellInside.toc();
    DynamicLabelList potentSurfCells(foundCells.size()/2);
    forAll(foundCells, cellI)
    {
        label cCell = foundCells[cellI];

        if(cellInside[cCell])
        {
            const labelList& neigh = cachedNeighbours_()[cCell];

            bool anyOutside = false;
            forAll(neigh, neighI)
            {
                if(!cellInside[neigh[neighI]])
                {
                    anyOutside = true;
                    break;
                }
            }

            if(anyOutside)
            {
                potentSurfCells.append(cCell);
            }
            else
            {
                intCells_[Pstream::myProcNo()].append(cCell);
                ibPartialVolume_[Pstream::myProcNo()] += 1;
                body[cCell] = 1.0;
            }
        }
        else
        {
            potentSurfCells.append(cCell);
        }
    }
    return potentSurfCells;
}
//---------------------------------------------------------------------------//
void geomModel::correctSurfCells
(
    volScalarField& body,
    DynamicLabelList& potentSurfCells,
    HashTable<bool, label, Hash<label>>& cellInside,
    List<labelList>& cellPoints
)
{
    HashTable<bool, label, Hash<label>> verticesStatus(potentSurfCells.size()*6);
    forAll(potentSurfCells, cellLabel)
    {    
        label cCell = potentSurfCells[cellLabel];
        scalar cBody(0);

        if(cellInside[cCell])
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

        if (cBody > thrSurf_)
        {
            if (cBody > (1.0-thrSurf_))
            {
                intCells_[Pstream::myProcNo()].append(cCell);
            }
            else
            {
                surfCells_[Pstream::myProcNo()].append(cCell);
            }
            ibPartialVolume_[Pstream::myProcNo()] += 1;

            body[cCell] += cBody;

            // clip the body field values
            body[cCell] = min(max(0.0,body[cCell]),1.0);
        }
    }
}
//---------------------------------------------------------------------------//
List<boundBox*> geomModel::getBBoxes()
{   
    boundBox cBBox = getBounds();
    bBox_.min() = cBBox.min();
    bBox_.max() = cBBox.max();

    List<boundBox*> retList(1);
    retList[0] = &bBox_;
    return retList;
}
//---------------------------------------------------------------------------//
