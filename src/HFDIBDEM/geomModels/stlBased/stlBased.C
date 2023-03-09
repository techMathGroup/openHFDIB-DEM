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
#include "stlBased.H"

using namespace Foam;

//---------------------------------------------------------------------------//
stlBased::stlBased
(
    const  fvMesh&   mesh,
    const contactType cType,
    word      stlPath,
    scalar  thrSurf
)
:
geomModel(mesh,cType,thrSurf),
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

    //meshSearch searchEng(mesh_);
    pointField bSMeshPts = bodySurfMesh_.points();

    // get its center of mass
    vector CoM(vector::zero);
    forAll(bSMeshPts,point)
    {
        CoM += bSMeshPts[point];
    }
    CoM/= bSMeshPts.size();

    const vector validDirs = (geometricD + vector::one)/2;
    vector dirCorr(cmptMultiply((vector::one - validDirs),CoM));
    dirCorr += cmptMultiply((vector::one - validDirs),0.5*(mesh_.bounds().max() + mesh_.bounds().min()));

    boundBox bodySurfBounds(bSMeshPts);
    // compute the max scales to stay in active bounding box
    vector maxScales(cellZoneBounds.max() - bodySurfBounds.max());
    maxScales -= cellZoneBounds.min() - bodySurfBounds.min();
    maxScales *= 0.5*0.9;//0.Y is there just to be sure

    InfoH << addModel_Info << "-- addModelMessage-- "
        << "acceptable movements: " << maxScales << endl;

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
volumeType stlBased::getVolumeType(subVolume& sv, bool cIb)
{
    const indexedOctree<treeDataTriSurface>& tree = triSurfSearch_->tree();
    const treeDataTriSurface& shapes = tree.shapes();

    autoPtr<labelList>& shapesIn = sv.getVolumeInfo(cIb).shapesIn_;

    if (shapesIn.empty())
    {
        std::shared_ptr<subVolume> parentSV = sv.parentSV();
        if (parentSV)
        {
            const labelList& parentShapesIn = *(parentSV->getVolumeInfo(cIb).shapesIn_);
            labelHashSet shapesInSV;

            forAll(parentShapesIn, i)
            {
                label shapeI = parentShapesIn[i];
                if (shapes.overlaps(shapeI, sv))
                {
                    shapesInSV.insert(shapeI);
                }
            }

            shapesIn.reset(new labelList(shapesInSV.toc()));
        }
        else
        {
            shapesIn.reset(new labelList(tree.findBox(sv)));
        }
    }

    if (shapesIn->size() > 0)
    {
        return volumeType::mixed;
    }

    return pointInside(sv.midpoint())
            ? volumeType::inside : volumeType::outside;
}
//---------------------------------------------------------------------------//
bool stlBased::limitFinalSubVolume
(
    const subVolume& sv,
    bool cIb,
    boundBox& limBBox
)
{
    limBBox = boundBox(sv.min(), sv.max());
    const autoPtr<labelList>& shapesIn = sv.getVolumeInfo(cIb).shapesIn_;
    if (shapesIn.empty())
    {
        return false;
    }

    vector normal = vector::zero;
    DynamicPointList intersectionPoints;
    forAll(*shapesIn, i)
    {
        getIntersectionPoints((*shapesIn)[i], sv, intersectionPoints);
        normal += (*triSurf_)[(*shapesIn)[i]].area(triSurf_->points());
    }

    if (intersectionPoints.size() == 0)
    {
        return false;
    }

    point closestPoint = vector::zero;
    forAll(intersectionPoints, i)
    {
        closestPoint += intersectionPoints[i];
    }
    closestPoint /= intersectionPoints.size();

    // Round normal to closest axis
    scalar magX = mag(normal.x());
    scalar magY = mag(normal.y());
    scalar magZ = mag(normal.z());

    if (magX > magY && magX > magZ)
    {
        if (sign(normal.x()) < 0)
        {
            limBBox.min().x() = closestPoint.x();
            return true;
        }
        else
        {
            limBBox.max().x() = closestPoint.x();
            return true;
        }
    }
    else if (magY > magX && magY > magZ)
    {
        if (sign(normal.y()) < 0)
        {
            limBBox.min().y() = closestPoint.y();
            return true;
        }
        else
        {
            limBBox.max().y() = closestPoint.y();
            return true;
        }
    }
    else if (magZ > magX && magZ > magY)
    {
        if (sign(normal.z()) < 0)
        {
            limBBox.min().z() = closestPoint.z();
            return true;
        }
        else
        {
            limBBox.max().z() = closestPoint.z();
            return true;
        }
    }

    return false;
    // Info << "sv: " << sv << " shapesIn size: " << shapesIn->size() << endl;
    // Info << "intersectionPoints: " << intersectionPoints << endl;

    // scalar nearestDistSqr = GREAT;
    // label minIndex = -1;
    // closestPoint = point::max;

    // const indexedOctree<treeDataTriSurface>& tree = triSurfSearch_->tree();
    // treeDataTriSurface::findNearestOp nearestOp(tree);
    // nearestOp(
    //     *shapesIn,
    //     sv.midpoint(),
    //     nearestDistSqr,
    //     minIndex,
    //     closestPoint
    // );

    // if (minIndex == -1)
    // {
    //     Info << "minIndex is -1" << endl;
    //     return false;
    // }

    // normal = (*triSurf_)[minIndex].area(triSurf_->points());

    // if (intersectionPoints.size() == 0)
    // {
    //     Info << "intersectionPoints is empty" << endl;
    //     return false;
    // }

    // forAll(intersectionPoints, i)
    // {
    //     closestPoint += intersectionPoints[i];
    // }
    // closestPoint /= intersectionPoints.size();

    // Info << "cIb: " << cIb << " sv: " << sv << " normal: " << normal << " nearest: " << closestPoint << " intersectionPoints: " << intersectionPoints << endl;

    // return true;
}
//---------------------------------------------------------------------------//
void stlBased::getIntersectionPoints
(
    const label index,
    const treeBoundBox& cubeBb,
    DynamicPointList& intersectionPoints
)
{
    const pointField& points = triSurf_->points();
    const typename triSurface::FaceType& f = (*triSurf_)[index];

    for (auto ind : f)
    {
        if (cubeBb.contains(points[ind]))
        {
            intersectionPoints.append(points[ind]);
        }
    }

    const point fc = f.centre(points);

    if (f.size() == 3)
    {
        return intersectBb
        (
            points[f[0]],
            points[f[1]],
            points[f[2]],
            cubeBb,
            intersectionPoints
        );
    }
    else
    {
        forAll(f, fp)
        {
            intersectBb
            (
                points[f[fp]],
                points[f[f.fcIndex(fp)]],
                fc,
                cubeBb,
                intersectionPoints
            );
        }
    }

    return;
}
//---------------------------------------------------------------------------//
void stlBased::intersectBb
(
    const point& p0,
    const point& p1,
    const point& p2,
    const treeBoundBox& cubeBb,
    DynamicPointList& intersectionPoints
)
{
    const vector p10 = p1 - p0;
    const vector p20 = p2 - p0;

    // cubeBb points; counted as if cell with vertex0 at cubeBb.min().
    const point& min = cubeBb.min();
    const point& max = cubeBb.max();

    const point& cube0 = min;
    const point  cube1(min.x(), min.y(), max.z());
    const point  cube2(max.x(), min.y(), max.z());
    const point  cube3(max.x(), min.y(), min.z());

    const point  cube4(min.x(), max.y(), min.z());
    const point  cube5(min.x(), max.y(), max.z());
    const point  cube7(max.x(), max.y(), min.z());

    //
    // Intersect all 12 edges of cube with triangle
    //

    point pInter;
    pointField origin(4);
    // edges in x direction
    origin[0] = cube0;
    origin[1] = cube1;
    origin[2] = cube5;
    origin[3] = cube4;

    scalar maxSx = max.x() - min.x();

    if (triangleFuncs::intersectAxesBundle(p0, p10, p20, 0, origin, maxSx, pInter))
    {
        intersectionPoints.append(pInter);
    }

    // edges in y direction
    origin[0] = cube0;
    origin[1] = cube1;
    origin[2] = cube2;
    origin[3] = cube3;

    scalar maxSy = max.y() - min.y();

    if (triangleFuncs::intersectAxesBundle(p0, p10, p20, 1, origin, maxSy, pInter))
    {
        intersectionPoints.append(pInter);
    }

    // edges in z direction
    origin[0] = cube0;
    origin[1] = cube3;
    origin[2] = cube7;
    origin[3] = cube4;

    scalar maxSz = max.z() - min.z();

    if (triangleFuncs::intersectAxesBundle(p0, p10, p20, 2, origin, maxSz, pInter))
    {
        intersectionPoints.append(pInter);
    }


    // Intersect triangle edges with bounding box
    if (cubeBb.intersects(p0, p1, pInter))
    {
        intersectionPoints.append(pInter);
    }
    if (cubeBb.intersects(p1, p2, pInter))
    {
        intersectionPoints.append(pInter);
    }
    if (cubeBb.intersects(p2, p0, pInter))
    {
        intersectionPoints.append(pInter);
    }
}
//---------------------------------------------------------------------------//
