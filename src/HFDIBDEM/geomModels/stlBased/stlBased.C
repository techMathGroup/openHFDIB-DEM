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
    // dominant-plane rule: the facet with the largest in-box
    // triangle-box overlap area approximates the surface crossing the
    // subVolume (aka leaf);
    // Note (MI): shapesIn_ must be filled - getVolumeType was called on
    //             this leaf before
    const auto& info = sv.getVolumeInfo(cIb);

    if (!info.shapesIn_.valid() || info.shapesIn_->size() == 0)
    {
        return false;
    }

    const labelList& shapesIn(info.shapesIn_());
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
            max(sv.min(), triBBox.min()),
            min(sv.max(), triBBox.max())
        );

        if (overlap.valid() && overlap.volume() > bestArea)
        {
            bestArea = overlap.volume();                                //largest dominates
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
    const vector offset(nOut*(sv.mag()/8.0 + VSMALL));

    // probe the two sides of the plane: whichever lands inside the
    // body tells the inward direction
    const bool posInside(pointInside(p + offset));
    const bool negInside(pointInside(p - offset));

    if (posInside == negInside)
    {
        // both or neither inside => no reliable plane
        return false;
    }

    n = posInside ? -nOut : nOut;

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
void stlBased::coverSplit
(
    const labelList& tris,
    scalar stallCoeff,
    label& nLeaves,
    label maxBoxes,
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
    // geometry fills this box - stop splitting here. The box cap is
    // checked against the leaf budget including both prospective
    // children so nLeaves never exceeds maxBoxes
    const scalar parentVol(parentBBox.mag());
    const scalar childrenVol
    (
        triSetBBox(leftSet).mag() + triSetBBox(rightSet).mag()
    );

    if
    (
        childrenVol > stallCoeff*parentVol
        || nLeaves + 2 > maxBoxes
    )
    {
        leaves[nLeaves++] = tris;
        return;
    }

    coverSplit(leftSet, stallCoeff, nLeaves, maxBoxes, leaves);
    coverSplit(rightSet, stallCoeff, nLeaves, maxBoxes, leaves);
}
//---------------------------------------------------------------------------//
bool stlBased::computeStaticCover(label maxBoxes, scalar stallCoeff)
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

    // build-time self-check: every surface point of every triangle
    // must lie inside its own leaf box (hence inside the union)
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

    coverActive_ = true;

    // cover statistics: total box volume vs the full AABB
    scalar coverVol(0);
    for (label leafI = 0; leafI < nLeaves; ++leafI)
    {
        coverVol += coverBoxes_[leafI]->mag();
    }

    InfoH << basic_Info << "cover built for " << stlPath_
        << ": " << nLeaves << " boxes, total volume " << coverVol
        << " vs AABB volume " << getBounds().mag()
        << " (ratio " << coverVol/max(getBounds().mag(), SMALL)
        << ")" << endl;

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
