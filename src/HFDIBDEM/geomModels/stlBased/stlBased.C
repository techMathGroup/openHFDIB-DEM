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

    openHFDIB-DEM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (Version 3) as published
    by the Free Software Foundation.

    openHFDIB-DEM is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with openHFDIB-DEM. If not, see <http://www.gnu.org/licenses/>.

InNamespace
    Foam

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/
#include "stlBased.H"

#include <fstream>                                                      //std::ofstream cover visualization
#include <iomanip>                                                      //std::setprecision cover visualization

using namespace Foam;

//---------------------------------------------------------------------------//
stlBased::stlBased
(
    const  fvMesh&   mesh,
    const contactType cType,
    word      stlPath
)
:
geomModel(mesh,cType),
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
stlPath_(stlPath),
coverActive_(false)
{
    historyPoints_ = bodySurfMesh_.points();
    triSurf_.reset(new triSurface(bodySurfMesh_));                      //OF.com: set -> reset
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
    computeVolumeCoM(triSurf_());
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
    const pointField& bSMeshPts = bodySurfMesh_.points();

    // get its center of mass
    vector CoM(getCoM());

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
        ranNum = 2.0*maxScales[i]*randGen.sample01<scalar>() - 1.0*maxScales[i];
        ranVec[i] = ranNum;
    }

    ranVec = cmptMultiply(validDirs,ranVec);                            // translate only with respect to valid directions
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
    bodyFieldValid_ = false;                                            // points moved: cell lists and the body field are stale

    CoM_ += translVec;
}
//---------------------------------------------------------------------------//
void stlBased::bodyScalePoints
(
    scalar scaleFac
)
{
    pointField bodyPoints(bodySurfMesh_.points());

    // scale about the center of mass so that CoM_ stays valid
    vector CoM(getCoM());

    bodyPoints -= CoM;
    bodySurfMesh_.movePoints(bodyPoints);
    bodySurfMesh_.scalePoints(scaleFac);
    bodyPoints = bodySurfMesh_.points();
    bodyPoints += CoM;
    bodySurfMesh_.movePoints(bodyPoints);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
    bodyFieldValid_ = false;                                            // points moved/scaled/rotated/synced
}
//---------------------------------------------------------------------------//
void stlBased::bodyRotatePoints
(
    scalar rotAngle,
    vector axisOfRot
)
{
    pointField bodyPoints(bodySurfMesh_.points());
    // rotate about the tracked center of mass
    vector CoM(getCoM());

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
    bodyFieldValid_ = false;                                            // points moved/scaled/rotated/synced
}
//---------------------------------------------------------------------------//
void stlBased::synchronPos(label owner)
{
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    owner = (owner == -1) ? owner_ : owner;

    if (owner == Pstream::myProcNo())
    {
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            UOPstream send(proci, pBufs);
            send << bodySurfMesh_.points();
        }
    }

    pBufs.finishedSends();
    // move body to points calculated by owner_
    UIPstream recv(owner, pBufs);
    pointField bodyPoints (recv);

    // move mesh
    bodySurfMesh_.movePoints(bodyPoints);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
    bodyFieldValid_ = false;                                            // points moved/scaled/rotated/synced

    // re-track the centroid
    computeVolumeCoM(triSurf_());                                       // points were replaced wholesale
}
//---------------------------------------------------------------------------//
void stlBased::getClosestPointAndNormal
(
    const point& startPoint,
    const vector& span,
    point& closestPoint,
    vector& normal
)
{
    // get nearest point on surface from contact center
    pointIndexHit ibPointIndexHit = triSurfSearch_().nearest(startPoint, span);
    List<pointIndexHit> ibPointIndexHitList(1,ibPointIndexHit);
    vectorField normalVectorField;

    // get contact normal direction
    const triSurfaceMesh& ibTempMesh( bodySurfMesh_);
    ibTempMesh.getNormal(ibPointIndexHitList,normalVectorField);

    if(ibPointIndexHit.hit())
    {
        normal = normalVectorField[0];
        closestPoint = ibPointIndexHit.hitPoint();
    }
    else
    {
        InfoH << basic_Info << "triSurfSearch: Missing the closest point on stl surface" << endl;
        normal = startPoint - getCoM();
        closestPoint = getCoM();
    }
}
//---------------------------------------------------------------------------//
volumeType stlBased::getVolumeType(subVolume& sv, bool cIb)
{
    auto& info = sv.getVolumeInfo(cIb);

    if (!info.shapesIn_.valid())
    {
        const indexedOctree<treeDataTriSurface>& tree = triSurfSearch_->tree();

        std::shared_ptr<subVolume> parentSV = sv.parentSV();
        if (parentSV)
        {
            const auto& infoParent = parentSV->getVolumeInfo(cIb);
            if (infoParent.shapesIn_.valid())
            {
                const labelList& parentShapesIn = infoParent.shapesIn_();
                if (parentShapesIn.size() == 0)
                {
                    info.shapesIn_.reset(new labelList(0));
                }
                else
                {
                    const treeDataTriSurface& shapes = tree.shapes();
                    labelList shapesIn(parentShapesIn.size());
                    label nShapesIn = 0;

                    forAll(parentShapesIn, i)
                    {
                        const label shapeI = parentShapesIn[i];
                        if (shapes.overlaps(shapeI, sv))
                        {
                            shapesIn[nShapesIn++] = shapeI;
                        }
                    }

                    shapesIn.setSize(nShapesIn);
                    info.shapesIn_.reset(new labelList(shapesIn));
                }
            }
            else
            {
                info.shapesIn_.reset(new labelList(tree.findBox(sv)));
            }
        }
        else
        {
            info.shapesIn_.reset(new labelList(tree.findBox(sv)));
        }
    }

    const labelList& shapesIn = info.shapesIn_();

    if (shapesIn.size() > 0)                                            //OF.com: mixed, inside, outside -> MIXED, INSIDE, OUTSIDE
    {
        return volumeType::MIXED;
    }

    return pointInside(sv.midpoint())
            ? volumeType::INSIDE : volumeType::OUTSIDE;
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
    if (shapesIn->empty())
    {
        return false;
    }

    vector normal = vector::zero;
    DynamicPointList intersectionPoints;
    forAll(*shapesIn, i)
    {
        getIntersectionPoints((*shapesIn)[i], sv, intersectionPoints);
        normal += (*triSurf_)[(*shapesIn)[i]].areaNormal(triSurf_->points());
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
bool stlBased::getLeafSubVolumePlane
(
    subVolume& sv,
    bool cIb,
    point& p,
    vector& n
)
{
    // Note (MI): shapesIn_ must be filled - getVolumeType was called on
    //             this leaf before
    const auto& info = sv.getVolumeInfo(cIb);

    if (!info.shapesIn_.valid() || info.shapesIn_->size() == 0)
    {
        return false;
    }

    return getShapesSurfacePlane(info.shapesIn_(), sv, p, n);
}
//---------------------------------------------------------------------------//
bool stlBased::getBoxSurfacePlane
(
    const boundBox& leaf,
    point& p,
    vector& n
)
{
    // wall-path variant: no octree parent chain, so search the
    // facets overlapping the leaf box directly
    const indexedOctree<treeDataTriSurface>& tree = triSurfSearch_->tree();

    const labelList shapesIn(tree.findBox(treeBoundBox(leaf)));

    if (shapesIn.size() == 0)                                           //internal finds nothing -> returns full leaf volume
    {
        return false;
    }

    return getShapesSurfacePlane(shapesIn, leaf, p, n);
}
//---------------------------------------------------------------------------//
bool stlBased::getShapesSurfacePlane
(
    const labelList& shapesIn,
    const boundBox& leaf,
    point& p,
    vector& n
)
{
    // dominant-plane rule: the facet with the largest in-box
    // triangle-box overlap area approximates the surface crossing the
    // leaf
    const pointField& surfPts(triSurf_->points());

    label bestFacet(-1);
    scalar bestArea(0);

    forAll(shapesIn, i)
    {
        const label facetI(shapesIn[i]);
        const triPointRef tri
        (
            surfPts[(*triSurf_)[facetI][0]],
            surfPts[(*triSurf_)[facetI][1]],
            surfPts[(*triSurf_)[facetI][2]]
        );

        // overlap weight: triangle-bbox vs leaf-box overlap volume;
        const boundBox triBBox
        (
            min(min(tri.a(), tri.b()), tri.c()),
            max(max(tri.a(), tri.b()), tri.c())
        );
        const boundBox overlap
        (
            max(leaf.min(), triBBox.min()),
            min(leaf.max(), triBBox.max())
        );

        if (overlap.valid() && overlap.volume() > bestArea)
        {
            bestArea = overlap.volume();                                //largest dominates
            bestFacet = facetI;
        }
        else if
        (
            overlap.valid()
         && bestFacet > -1
         && overlap.volume() == bestArea
         && facetI < bestFacet
        )
        {
            // bit-identical overlap ties are the norm on flat
            // faces of axis-aligned or lattice-rotated bodies
            // (every diagonal-split face pair shares one bbox);
            // the smallest facet index keeps the pick independent
            // of the candidate order the caller happens to use
            bestFacet = facetI;
        }
    }

    if (bestFacet == -1)
    {
        return false;
    }

    const vector areaN
    (
        (*triSurf_)[bestFacet].areaNormal(surfPts)
    );
    const scalar aNmag(mag(areaN));

    if (aNmag < VSMALL)
    {
        return false;
    }

    // outward facet normal; inward orientation is resolved with a
    // single pointInside probe offset into the leaf from the plane,
    // so the stl file orientation is not trusted
    // Note (MI): the cost of this "not trusted" should be checked
    //            -> if this becomes a bottleneck, just trust the user
    //           (stl orientation)
    vector nOut(areaN/aNmag);

    p = (*triSurf_)[bestFacet].centre(surfPts);
    const vector offset(nOut*(leaf.avgDim()/8.0 + VSMALL));

    // probe the two sides of the plane: whichever lands inside the
    // body tells the inward direction
    const bool posInside(pointInside(p + offset));
    const bool negInside(pointInside(p - offset));

    if (posInside == negInside)
    {
        // both or neither inside => no reliable plane
        return false;
    }

    // n must point into the body: the planePolyClip half-spaces keep
    // the (x - p) & n >= 0 side, which is the body interior here
    n = posInside ? nOut : -nOut;

    return true;
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
void stlBased::setBodyPosition(pointField pos)
{
    // the DEM broadcast calls this for every body each subcycle with
    // the points gathered from rank 0; for static bodies
    // skip the rebuild
    if (pos == bodySurfMesh_.points())
    {
        return;
    }

    bodySurfMesh_.movePoints(pos);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
    bodyFieldValid_ = false;                                            // points moved: cell lists and the body field are stale

    // re-track the centroid
    computeVolumeCoM(triSurf_());                                       // points were replaced wholesale
}
//---------------------------------------------------------------------------//
boundBox stlBased::triSetBBox(const labelList& tris) const
{
    const pointField& points = triSurf_->points();

    point minP = point::max;
    point maxP = point::min;

    forAll(tris, i)
    {
        const typename triSurface::FaceType& f = (*triSurf_)[tris[i]];

        for (auto ind : f)
        {
            minP = min(minP, points[ind]);
            maxP = max(maxP, points[ind]);
        }
    }

    return boundBox(minP, maxP);
}
//---------------------------------------------------------------------------//
scalar stlBased::coverBoxVol(const boundBox& bBox) const
{
    // direction-masked volume: geometricD is +1 in active, -1 in
    // empty directions, so the empty component becomes 1
    vector span(bBox.span());
    const vector validDirs((geometricD + vector::one)/2);

    for (label dir = 0; dir < 3; ++dir)
    {
        if (validDirs[dir] < 0.5)
        {
            span[dir] = 1.0;
        }
    }

    return cmptProduct(span);
}
//---------------------------------------------------------------------------//
void stlBased::coverSplit
(
    const labelList& tris,
    scalar stallCoeff,
    label& nLeaves,
    label budget,
    List<labelList>& leaves
) const
{
    const boundBox parentBBox(triSetBBox(tris));

    const vector span(parentBBox.max() - parentBBox.min());
    const vector validDirs((geometricD + vector::one)/2);
    const vector spanDirs(cmptMultiply(span, validDirs));

    label splitDir(-1);
    scalar widest(0);
    for (label dir = 0; dir < 3; ++dir)
    {
        if (spanDirs[dir] > widest)
        {
            widest = spanDirs[dir];
            splitDir = dir;
        }
    }

    // cannot split further in a non-degenerate direction
    if (splitDir == -1)
    {
        leaves[nLeaves++] = tris;
        return;
    }

    // a set of one triangle cannot be halved into two non-empty
    // children - it is a leaf (an empty child would produce an
    // inverted boundBox and overflow in the volume computation)
    if (tris.size() < 2)
    {
        leaves[nLeaves++] = tris;
        return;
    }

    // budget exhausted: this subtree must emit exactly one box
    // (the budget is decremented by the emitted leaves, so the
    // write below is always in bounds)
    if (budget < 2)
    {
        leaves[nLeaves++] = tris;
        return;
    }

    // split at the median triangle centroid along the widest axis;
    // ties are broken by the triangle INDEX (sortedOrder is stable
    // w.r.t. it), never by float comparison of equal centroids, so
    // every rank builds the identical partition
    const pointField& points = triSurf_->points();
    scalarField centroidCoord(tris.size());

    forAll(tris, i)
    {
        centroidCoord[i] =
            (*triSurf_)[tris[i]].centre(points)[splitDir];
    }

    labelList order(sortedOrder(centroidCoord));

    labelList leftSet(order.size()/2);
    labelList rightSet(order.size() - order.size()/2);

    forAll(order, i)
    {
        if (i < order.size()/2)
        {
            leftSet[i] = tris[order[i]];
        }
        else
        {
            rightSet[i - order.size()/2] = tris[order[i]];
        }
    }

    // stall criterion: if the children barely reduce the bound, the
    // geometry fills this box - stop splitting here
    const scalar parentVol(coverBoxVol(parentBBox));
    const scalar childrenVol
    (
        coverBoxVol(triSetBBox(leftSet)) + coverBoxVol(triSetBBox(rightSet))
    );

    if (childrenVol > stallCoeff*parentVol)
    {
        leaves[nLeaves++] = tris;
        return;
    }

    // hand half the budget to the left subtree, the remainder to
    // the right: both always get at least one slot, and the total
    // emitted can never exceed the parent budget (integer split,
    // deterministic across ranks)
    const label leftBudget((budget + 1)/2);
    const label nLeavesBefore(nLeaves);

    coverSplit(leftSet, stallCoeff, nLeaves, leftBudget, leaves);
    coverSplit
    (
        rightSet,
        stallCoeff,
        nLeaves,
        budget - (nLeaves - nLeavesBefore),
        leaves
    );
}
//---------------------------------------------------------------------------//
scalar stlBased::boxUnionVol(const List<boundBox>& boxes) const
{
    // sweep along x over all box x-boundaries; inside each slab
    // take the exact union area of the yz rectangles
    if (boxes.size() == 0)
    {
        return 0.0;
    }

    scalarField xs(2*boxes.size());
    forAll(boxes, bI)
    {
        xs[2*bI] = boxes[bI].min().x();
        xs[2*bI + 1] = boxes[bI].max().x();
    }
    sort(xs);

    scalar totalVol(0);

    for (label i = 0; i < xs.size() - 1; ++i)
    {
        if (xs[i + 1] - xs[i] <= VSMALL)
        {
            continue;
        }

        // yz rectangles of every box spanning this slab
        List<boundBox> rects;
        forAll(boxes, bI)
        {
            if
            (
                boxes[bI].min().x() <= xs[i]
             && boxes[bI].max().x() >= xs[i + 1]
            )
            {
                rects.append
                (
                    boundBox
                    (
                        point(boxes[bI].min().y(), boxes[bI].min().z(), 0),
                        point(boxes[bI].max().y(), boxes[bI].max().z(), 0)
                    )
                );
            }
        }

        if (rects.size() == 0)
        {
            continue;
        }

        // exact union area of 2-d rectangles: sweep along y,
        // merge z intervals per y slab
        scalarField ys(2*rects.size());
        forAll(rects, rI)
        {
            ys[2*rI] = rects[rI].min().x();
            ys[2*rI + 1] = rects[rI].max().x();
        }
        sort(ys);

        scalar area(0);

        for (label j = 0; j < ys.size() - 1; ++j)
        {
            if (ys[j + 1] - ys[j] <= VSMALL)
            {
                continue;
            }

            // z intervals of every rectangle spanning this
            // y slab, merged in sorted order
            List<Tuple2<scalar,scalar>> zInts;
            forAll(rects, rI)
            {
                if
                (
                    rects[rI].min().x() <= ys[j]
                 && rects[rI].max().x() >= ys[j + 1]
                )
                {
                    zInts.append
                    (
                        Tuple2<scalar,scalar>
                        (
                            rects[rI].min().y(),
                            rects[rI].max().y()
                        )
                    );
                }
            }

            if (zInts.size() == 0)
            {
                continue;
            }

            sort(zInts);
            scalar zLo(zInts[0].first());
            scalar zHi(zInts[0].second());
            scalar merged(0);

            for (label k = 1; k < zInts.size(); ++k)
            {
                if (zInts[k].first() <= zHi)
                {
                    zHi = max(zHi, zInts[k].second());
                }
                else
                {
                    merged += zHi - zLo;
                    zLo = zInts[k].first();
                    zHi = zInts[k].second();
                }
            }
            merged += zHi - zLo;

            area += (ys[j + 1] - ys[j])*merged;
        }

        totalVol += (xs[i + 1] - xs[i])*area;
    }

    return totalVol;
}
//---------------------------------------------------------------------------//
void stlBased::spatialSplit
(
    const labelList& tris,
    const boundBox& clampBox,
    label& nLeaves,
    label budget,
    label depth,
    bool prune,
    List<labelList>& leaves,
    List<boundBox>& leafBoxes
) const
{
    // effective node box: triangle-set bbox clamped by the
    // half-space chain of all ancestors (cumulative clamping)
    const boundBox triBox(triSetBBox(tris));
    const boundBox nodeBox
    (
        max(triBox.min(), clampBox.min()),
        min(triBox.max(), clampBox.max())
    );

    const vector span(nodeBox.span());
    const vector validDirs((geometricD + vector::one)/2);

    if (budget < 2)
    {
        leaves[nLeaves] = tris;
        leafBoxes[nLeaves] = nodeBox;
        ++nLeaves;
        return;
    }

    // valid split axes (non-empty, non-degenerate span)
    labelList axes;
    for (label dir = 0; dir < 3; ++dir)
    {
        if (validDirs[dir] > 0.5 && span[dir] > 0)
        {
            axes.append(dir);
        }
    }

    if (axes.size() == 0)
    {
        leaves[nLeaves] = tris;
        leafBoxes[nLeaves] = nodeBox;
        ++nLeaves;
        return;
    }

    // cyclic rule: rotate through the valid axes by depth so the
    // thin axis gets its turn regardless of local aspect ratio
    const label splitDir(axes[depth % axes.size()]);
    const scalar mid(0.5*(nodeBox.min()[splitDir]
        + nodeBox.max()[splitDir]));

    // classify: a triangle touches a half when its bbox overlaps
    // the half-interval along the split axis; straddlers are
    // duplicated into both halves
    const pointField& points = triSurf_->points();

    labelList leftTris(tris.size());
    labelList rightTris(tris.size());
    label nLeft(0);
    label nRight(0);

    forAll(tris, i)
    {
        const typename triSurface::FaceType& f = (*triSurf_)[tris[i]];
        scalar tMin(GREAT);
        scalar tMax(-GREAT);

        for (auto ind : f)
        {
            tMin = min(tMin, points[ind][splitDir]);
            tMax = max(tMax, points[ind][splitDir]);
        }

        if (tMax <= mid)
        {
            leftTris[nLeft++] = tris[i];
        }
        else if (tMin >= mid)
        {
            rightTris[nRight++] = tris[i];
        }
        else
        {
            leftTris[nLeft++] = tris[i];
            rightTris[nRight++] = tris[i];
        }
    }

    leftTris.setSize(nLeft);
    rightTris.setSize(nRight);

    if (nLeft == 0 || nRight == 0)
    {
        // criterion (a): empty half - the body fills only one
        // side; splitting cannot tighten anything here
        leaves[nLeaves] = tris;
        leafBoxes[nLeaves] = nodeBox;
        ++nLeaves;
        return;
    }

    // child clamp boxes: the node box cut at the split plane
    boundBox leftClamp(nodeBox);
    boundBox rightClamp(nodeBox);
    leftClamp.max()[splitDir] = mid;
    rightClamp.min()[splitDir] = mid;

    const label nLeavesBefore(nLeaves);
    spatialSplit
    (
        leftTris,
        leftClamp,
        nLeaves,
        (budget + 1)/2,
        depth + 1,
        prune,
        leaves,
        leafBoxes
    );
    spatialSplit
    (
        rightTris,
        rightClamp,
        nLeaves,
        budget - (nLeaves - nLeavesBefore),
        depth + 1,
        prune,
        leaves,
        leafBoxes
    );

    if (prune)
    {
        // bottom-up prune: if the subtree's leaf-box union does
        // not improve on this node's own box, collapse it - the
        // boxes only tile the node box without tightening
        List<boundBox> subBoxes
        (
            leafBoxes.slice
            (
                nLeavesBefore,
                nLeaves - nLeavesBefore
            )
        );

        if
        (
            boxUnionVol(subBoxes)
         >= coverBoxVol(nodeBox)*(1.0 - SMALL)
        )
        {
            // roll the subtree back to a single leaf
            leaves[nLeavesBefore] = tris;
            leafBoxes[nLeavesBefore] = nodeBox;
            nLeaves = nLeavesBefore + 1;
        }
    }
}
//---------------------------------------------------------------------------//
bool stlBased::computeStaticCover
(
    label maxBoxes,
    scalar stallCoeff,
    bool writeCover,
    word bodyName,
    word coverAlgorithm,
    word axisRule,
    bool prune
)
{
    if (coverActive_)
    {
        return true;
    }

    const label nTris(triSurf_->size());

    if (nTris < 1)
    {
        return false;
    }

    labelList allTris(nTris);
    forAll(allTris, i)
    {
        allTris[i] = i;
    }

    List<labelList> leaves(maxBoxes);
    label nLeaves(0);

    if (coverAlgorithm == "spatial")
    {
        // spatial-halving tree with clamped leaf boxes; the
        // leaves/leafBoxes arrays are exactly maxBoxes long (the
        // budget bounds the emitted leaves), and prune can only
        // shrink nLeaves further
        List<boundBox> leafBoxes(maxBoxes);
        spatialSplit
        (
            allTris,
            getBounds(),
            nLeaves,
            maxBoxes,
            0,
            prune,
            leaves,
            leafBoxes
        );

        coverPartition_.setSize(nLeaves);
        coverBoxes_.setSize(nLeaves);

        for (label leafI = 0; leafI < nLeaves; ++leafI)
        {
            coverPartition_[leafI] = leaves[leafI];
            coverBoxes_[leafI] = std::make_shared<boundBox>
            (
                leafBoxes[leafI]
            );
        }
    }
    else
    {
        coverSplit(allTris, stallCoeff, nLeaves, maxBoxes, leaves);

        coverPartition_.setSize(nLeaves);
        coverBoxes_.setSize(nLeaves);

        for (label leafI = 0; leafI < nLeaves; ++leafI)
        {
            coverPartition_[leafI] = leaves[leafI];
            coverBoxes_[leafI] = std::make_shared<boundBox>
            (
                triSetBBox(leaves[leafI])
            );
        }
    }

    // build-time self-check: every surface point must lie inside
    // the union of the cover boxes. kdTree: a point of a triangle
    // lies in its own leaf box. spatial: with duplication a
    // straddling triangle's point can fall in the sibling's
    // clamped box, so the check tests the union membership (the
    // actual broad-phase invariant)
    if (coverAlgorithm == "spatial")
    {
        const pointField& points = triSurf_->points();

        forAll(allTris, tI)
        {
            const typename triSurface::FaceType& f =
                (*triSurf_)[tI];

            for (auto ind : f)
            {
                bool inside(false);
                for (label leafI = 0; leafI < nLeaves; ++leafI)
                {
                    if (coverBoxes_[leafI]->contains(points[ind]))
                    {
                        inside = true;
                        break;
                    }
                }

                if (!inside)
                {
                    FatalErrorInFunction
                        << "cover self-check failed: point "
                        << points[ind] << " of triangle " << tI
                        << " outside the cover box union"
                        << abort(FatalError);
                }
            }
        }
    }
    else
    {
        for (label leafI = 0; leafI < nLeaves; ++leafI)
        {
            const boundBox& bBox(*coverBoxes_[leafI]);
            const pointField& points = triSurf_->points();

            forAll(coverPartition_[leafI], i)
            {
                const typename triSurface::FaceType& f =
                    (*triSurf_)[coverPartition_[leafI][i]];

                for (auto ind : f)
                {
                    if (!bBox.contains(points[ind]))
                    {
                        FatalErrorInFunction
                            << "cover self-check failed: point "
                            << points[ind] << " of triangle "
                            << coverPartition_[leafI][i]
                            << " outside its cover box " << bBox
                            << abort(FatalError);
                    }
                }
            }
        }
    }

    coverActive_ = true;

    // cover statistics: total box volume vs the full AABB
    scalar coverVol(0);
    for (label leafI = 0; leafI < nLeaves; ++leafI)
    {
        coverVol += coverBoxVol(*coverBoxes_[leafI]);
    }

    InfoH << basic_Info << "cover built for " << stlPath_
        << ": " << nLeaves << " boxes, total volume " << coverVol
        << " vs AABB volume " << coverBoxVol(getBounds())
        << " (ratio " << coverVol/max(coverBoxVol(getBounds()), SMALL)
        << ")" << endl;

    if (writeCover)
    {
        writeCoverVtk
        (
            bodyName,
            coverBoxes_,
            coverPartition_,
            triSurf_(),
            getBounds()
        );
    }

    return true;
}
//---------------------------------------------------------------------------//
void stlBased::refreshCoverBounds()
{
    const pointField& points = triSurf_->points();

    forAll(coverPartition_, leafI)
    {
        point minP = point::max;
        point maxP = point::min;

        forAll(coverPartition_[leafI], i)
        {
            const typename triSurface::FaceType& f =
                (*triSurf_)[coverPartition_[leafI][i]];

            for (auto ind : f)
            {
                minP = min(minP, points[ind]);
                maxP = max(maxP, points[ind]);
            }
        }

        // write through the aliased boxes in place - never
        // reallocate (in-place mutation contract)
        coverBoxes_[leafI]->min() = minP;
        coverBoxes_[leafI]->max() = maxP;
    }
}
//---------------------------------------------------------------------------//
void stlBased::writeCoverVtk
(
    const word& bodyName,
    const List<std::shared_ptr<boundBox>>& boxes,
    const labelListList& partition,
    const triSurface& surf,
    const boundBox& fullAABB
) const
{
    // master rank only: the deterministic build guarantees every
    // rank holds the identical partition, no gather is needed
    if (!Pstream::master())
    {
        return;
    }

    fileName coverDir
    (
        mesh_.time().rootPath() + "/"
        + mesh_.time().globalCaseName()
        + "/bodiesInfo/" + mesh_.time().timeName()
        + "/coverFiles"
    );

    if (!isDir(coverDir))
    {
        mkDir(coverDir);
    }

    const label nBoxes(boxes.size());
    const scalar fullVol(max(coverBoxVol(fullAABB), SMALL));

    // ---- one hexahedron per cover box ----------------
    {
        std::ofstream boxFile
        (
            (coverDir + "/" + bodyName + "_boxes.vtk").c_str(),
            std::ofstream::trunc
        );

        if (!boxFile.is_open())
        {
            FatalErrorInFunction
                << "cannot open cover file "
                << coverDir + "/" + bodyName + "_boxes.vtk"
                << " for writing"
                << abort(FatalError);
        }

        boxFile << std::setprecision(16);

        boxFile << "# vtk DataFile Version 3.0\n"
            << "cover boxes for " << bodyName << "\n"
            << "ASCII\n"
            << "DATASET UNSTRUCTURED_GRID\n"
            << "POINTS " << 8*nBoxes << " double\n";

        // 8 points per cell, contiguous layout, VTK hexahedron
        // order (counterclockwise quads)
        for (label boxI = 0; boxI < nBoxes; ++boxI)
        {
            const point& minP(boxes[boxI]->min());
            const point& maxP(boxes[boxI]->max());

            boxFile << minP.x() << " " << minP.y() << " " << minP.z() << "\n";
            boxFile << maxP.x() << " " << minP.y() << " " << minP.z() << "\n";
            boxFile << maxP.x() << " " << maxP.y() << " " << minP.z() << "\n";
            boxFile << minP.x() << " " << maxP.y() << " " << minP.z() << "\n";
            boxFile << minP.x() << " " << minP.y() << " " << maxP.z() << "\n";
            boxFile << maxP.x() << " " << minP.y() << " " << maxP.z() << "\n";
            boxFile << maxP.x() << " " << maxP.y() << " " << maxP.z() << "\n";
            boxFile << minP.x() << " " << maxP.y() << " " << maxP.z() << "\n";
        }

        boxFile << "CELLS " << nBoxes << " " << 9*nBoxes << "\n";

        for (label boxI = 0; boxI < nBoxes; ++boxI)
        {
            boxFile << " 8";

            for (label ptI = 0; ptI < 8; ++ptI)
            {
                boxFile << " " << 8*boxI + ptI;
            }

            boxFile << "\n";
        }

        boxFile << "CELL_TYPES " << nBoxes << "\n";

        for (label boxI = 0; boxI < nBoxes; ++boxI)
        {
            boxFile << "12\n";
        }

        boxFile << "CELL_DATA " << nBoxes << "\n";

        boxFile << "SCALARS boxIndex int 1\n"
            << "LOOKUP_TABLE default\n";

        for (label boxI = 0; boxI < nBoxes; ++boxI)
        {
            boxFile << boxI << "\n";
        }

        boxFile << "SCALARS nTriangles int 1\n"
            << "LOOKUP_TABLE default\n";

        for (label boxI = 0; boxI < nBoxes; ++boxI)
        {
            boxFile << partition[boxI].size() << "\n";
        }

        boxFile << "SCALARS boxVolume double 1\n"
            << "LOOKUP_TABLE default\n";

        for (label boxI = 0; boxI < nBoxes; ++boxI)
        {
            boxFile << coverBoxVol(*boxes[boxI]) << "\n";
        }

        boxFile << "SCALARS volumeRatio double 1\n"
            << "LOOKUP_TABLE default\n";

        for (label boxI = 0; boxI < nBoxes; ++boxI)
        {
            boxFile << coverBoxVol(*boxes[boxI])/fullVol << "\n";
        }
    }

    // ---- surface triangles colored by partition ------
    {
        // reverse map triangle -> box; the build-time self-check
        // guarantees every triangle feeds exactly one box, a
        // surviving -1 is a partition bug
        const label nTris(surf.size());
        labelList triBox(nTris, -1);

        forAll(partition, boxI)
        {
            forAll(partition[boxI], i)
            {
                triBox[partition[boxI][i]] = boxI;
            }
        }

        forAll(triBox, triI)
        {
            if (triBox[triI] < 0)
            {
                WarningInFunction
                    << "triangle " << triI << " of " << bodyName
                    << " lies in no cover box"
                    << endl;
            }
        }

        const pointField& points(surf.points());

        std::ofstream partFile
        (
            (coverDir + "/" + bodyName + "_partition.vtk").c_str(),
            std::ofstream::trunc
        );

        if (!partFile.is_open())
        {
            FatalErrorInFunction
                << "cannot open cover file "
                << coverDir + "/" + bodyName + "_partition.vtk"
                << " for writing"
                << abort(FatalError);
        }

        partFile << std::setprecision(16);

        partFile << "# vtk DataFile Version 3.0\n"
            << "cover partition for " << bodyName << "\n"
            << "ASCII\n"
            << "DATASET UNSTRUCTURED_GRID\n"
            << "POINTS " << points.size() << " double\n";

        forAll(points, ptI)
        {
            partFile << points[ptI].x() << " "
                << points[ptI].y() << " "
                << points[ptI].z() << "\n";
        }

        partFile << "CELLS " << nTris << " " << 4*nTris << "\n";

        for (label triI = 0; triI < nTris; ++triI)
        {
            const typename triSurface::FaceType& f = surf[triI];

            partFile << " 3";

            for (auto ind : f)
            {
                partFile << " " << ind;
            }

            partFile << "\n";
        }

        partFile << "CELL_TYPES " << nTris << "\n";

        for (label triI = 0; triI < nTris; ++triI)
        {
            partFile << "5\n";
        }

        partFile << "CELL_DATA " << nTris << "\n"
            << "SCALARS boxIndex int 1\n"
            << "LOOKUP_TABLE default\n";

        forAll(triBox, triI)
        {
            partFile << triBox[triI] << "\n";
        }
    }
}
//---------------------------------------------------------------------------//
List<std::shared_ptr<boundBox>> stlBased::getBBoxes()
{
    if (!coverActive_)
    {
        return geomModel::getBBoxes();
    }

    // aliasing-contract watch: the same boundBox objects must be
    // returned every call - verletPoints hold references into them
    if (coverBoxWatch_.size() != coverBoxes_.size())
    {
        coverBoxWatch_.setSize(coverBoxes_.size());
        forAll(coverBoxes_, boxI)
        {
            coverBoxWatch_[boxI] = coverBoxes_[boxI].get();
        }
    }
    else
    {
        forAll(coverBoxes_, boxI)
        {
            if (coverBoxWatch_[boxI] != coverBoxes_[boxI].get())
            {
                FatalErrorInFunction
                    << "cover boxes of " << stlPath_
                    << " were reallocated after verlet registration"
                    << " (in-place mutation contract violated - live"
                    << " verletPoints alias these boxes)"
                    << abort(FatalError);
            }
        }
    }

    refreshCoverBounds();

    return coverBoxes_;
}
//---------------------------------------------------------------------------//
