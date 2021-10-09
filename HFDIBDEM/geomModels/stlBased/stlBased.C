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
#include "stlBased.H"
#include "meshSearch.H"

using namespace Foam;

//---------------------------------------------------------------------------//
stlBased::stlBased
(
    const  dynamicFvMesh&   mesh,
    contactType cType,
    word      stlPath,
    scalar  thrSurf,
    Vector<label> geometricD
)
:
geomModel(mesh,cType,thrSurf,geometricD),
bodySurfMesh_
(
    IOobject
    (
        stlPath,
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
stlPath_(stlPath)
{
    historyPoints_ = bodySurfMesh_.points();
    triSurf_.set(new triSurface(bodySurfMesh_));
    triSurfSearch_.set(new triSurfaceSearch(triSurf_()));
}
//---------------------------------------------------------------------------//
vector stlBased::addModelReturnRandomPosition
(
    const bool allActiveCellsInMesh,
    const boundBox  cellZoneBounds,
    Random&          randGen
)
{
    vector ranVec(vector::zero);

    label nGeometricD(0);
    forAll (geometricD_, direction)
    {
        if (geometricD_[direction] == 1)
        {
            nGeometricD++;
        }
    }

    meshSearch searchEng(mesh_);
    pointField bSMeshPts = bodySurfMesh_.points();

    // get its center of mass
    vector CoM(vector::zero);
    forAll(bSMeshPts,point)
    {
        CoM += bSMeshPts[point];
    }
    CoM/= bSMeshPts.size();

    const vector validDirs = (geometricD_ + Vector<label>::one)/2;
    vector dirCorr(cmptMultiply((vector::one - validDirs),CoM));
    dirCorr += cmptMultiply((vector::one - validDirs),0.5*(mesh_.bounds().max() + mesh_.bounds().min()));

    boundBox bodySurfBounds(bSMeshPts);
    // compute the max scales to stay in active bounding box
    vector maxScales(cellZoneBounds.max() - bodySurfBounds.max());
    maxScales -= cellZoneBounds.min() - bodySurfBounds.min();
    maxScales *= 0.5*0.9;//0.Y is there just to be sure

    Info << "-- addModelMessage-- " << "acceptable movements: " << maxScales << endl;

    scalar ranNum = 0;
    for (int i=0;i<3;i++)
    {
        ranNum = 2.0*maxScales[i]*randGen.scalar01() - 1.0*maxScales[i];
        ranVec[i] = ranNum;
    }

    ranVec = cmptMultiply(validDirs,ranVec);                            //translate only with respect to valid directions
    ranVec += dirCorr;

    return ranVec;
}
//---------------------------------------------------------------------------//
void stlBased::bodyMovePoints
(
    vector translVec
)
{
    pointField bodyPoints(bodySurfMesh_.points());
    bodyPoints += translVec;

    bodySurfMesh_.movePoints(bodyPoints);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}
//---------------------------------------------------------------------------//
void stlBased::bodyScalePoints
(
    scalar scaleFac
)
{
    pointField bodyPoints(bodySurfMesh_.points());

    // get its center of mass
    vector CoM(vector::zero);
    forAll(bodyPoints,point)
    {
        CoM += bodyPoints[point];
    }
    CoM/= bodyPoints.size();

    bodyPoints -= CoM;
    bodySurfMesh_.movePoints(bodyPoints);
    bodySurfMesh_.scalePoints(scaleFac);
    bodyPoints = bodySurfMesh_.points();
    bodyPoints += CoM;
    bodySurfMesh_.movePoints(bodyPoints);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}
//---------------------------------------------------------------------------//
void stlBased::bodyRotatePoints
(
    scalar rotAngle,
    vector axisOfRot
)
{
    pointField bodyPoints(bodySurfMesh_.points());
    // get its center of mass
    vector CoM(vector::zero);
    forAll(bodyPoints,point)
    {
        CoM += bodyPoints[point];
    }
    CoM/= bodyPoints.size();

    tensor rotMatrix(Foam::cos(rotAngle)*tensor::I);

    rotMatrix += Foam::sin(rotAngle)*tensor(
            0.0,      -axisOfRot.z(),  axisOfRot.y(),
            axisOfRot.z(), 0.0,       -axisOfRot.x(),
        -axisOfRot.y(), axisOfRot.x(),  0.0
    );

    rotMatrix += (1.0-Foam::cos(rotAngle))*(axisOfRot * axisOfRot);

    bodyPoints -= CoM;
    bodyPoints = rotMatrix & bodyPoints;
    bodyPoints += CoM;
    bodySurfMesh_.movePoints(bodyPoints);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}
//---------------------------------------------------------------------------//
void stlBased::synchronPos()
{
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    if (owner_ == Pstream::myProcNo())
    {
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            UOPstream send(proci, pBufs);
            send << bodySurfMesh_.points();
        }
    }

    pBufs.finishedSends();
    // move body to points calculated by owner_
    UIPstream recv(owner_, pBufs);
    pointField bodyPoints (recv);

    // move mesh
    bodySurfMesh_.movePoints(bodyPoints);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}
//---------------------------------------------------------------------------//
