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
#include "geomModel.H"

using namespace Foam;

//---------------------------------------------------------------------------//
geomModel::geomModel
(
    const  dynamicFvMesh&   mesh,
    contactType cType,
    scalar  thrSurf,
    Vector<label> geometricD
)
:
mesh_(mesh),
contactType_(cType),
owner_(0),
cellToStartInCreateIB_(0),
thrSurf_(thrSurf),
geometricD_(geometricD),
intSpan_(2.0),
sdBasedLambda_(false),
curMeshBounds_(mesh_.points(),false),
lastIbPoints_(0),
M_(0.0),
M0_(0.0),
CoM_(vector::zero),
I_(symmTensor::zero),
dC_(0.0),
rhoS_("rho",dimensionSet(1,-3,0,0,0,0,0),1.0)
{
    ibPartialVolume_.setSize(Pstream::nProcs());
}
geomModel::~geomModel()
{
}
//---------------------------------------------------------------------------//
void geomModel::calculateGeometricalProperties( volScalarField& body,List<DynamicLabelList>& surfCells, List<DynamicLabelList>& intCells )
{
    vector CoMOld  = CoM_;
    M_      = scalar(0);
    CoM_    = vector::zero;
    I_      = symmTensor::zero;
    vector tmpCom(vector::zero);

    addToMAndI(body,surfCells[Pstream::myProcNo()],tmpCom, CoMOld);
    addToMAndI(body,intCells[Pstream::myProcNo()],tmpCom, CoMOld);

    // collect from processors
    reduce(M_, sumOp<scalar>());
    reduce(tmpCom,  sumOp<vector>());
    reduce(I_,  sumOp<symmTensor>());

    CoM_ = tmpCom / (M_+SMALL);
}
//---------------------------------------------------------------------------//
void geomModel::addToMAndI
(
    volScalarField& body,
    DynamicLabelList& labelCellLst,
    vector& tmpCom,
    vector CoMOld
)
{
    forAll (labelCellLst,cell)
    {
        label cellI  = labelCellLst[cell];
        
        scalar Mi    = body[cellI]*rhoS_.value()*mesh_.V()[cellI];
        // add to M_
        
        M_      += Mi;
        tmpCom  += Mi*mesh_.C()[cellI];
        scalar xLoc = mesh_.C()[cellI].x() - CoMOld.x();
        scalar yLoc = mesh_.C()[cellI].y() - CoMOld.y();
        scalar zLoc = mesh_.C()[cellI].z() - CoMOld.z();
        
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
void geomModel::computeBodyCharPars(List<DynamicLabelList>& surfCells)
{
    forAll (surfCells[Pstream::myProcNo()],sCellI)
    {
        label cellI = surfCells[Pstream::myProcNo()][sCellI];
        dC_ = max(dC_,mag(CoM_-mesh_.C()[cellI]));
    }
    M0_ = M_;
    reduce(dC_, maxOp<scalar>());
}
//---------------------------------------------------------------------------//
