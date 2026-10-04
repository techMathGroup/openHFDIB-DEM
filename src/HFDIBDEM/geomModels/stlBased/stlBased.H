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

Description
    class for addition and removal of immersedBodies

SourceFiles
    stlBased.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef stlBased_H
#define stlBased_H

#include "geomModel.H"
#include "momentOfInertia.H"                                            //reqiured for computeVolumeCoM

namespace Foam
{
    class stlBased
    : public geomModel
    {

        private:

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        protected:
        triSurfaceMesh      bodySurfMesh_;
        word                stlPath_;
        autoPtr<triSurface> triSurf_;
        autoPtr<triSurfaceSearch> triSurfSearch_;

        pointField historyPoints_;

        // tight multi-box cover (broad-phase bound): k-d median
        // partition of the surface triangles, built once per body
        // Note (MI): the shared boundBox objects in coverBoxes_ are
        // handed to verlet registration and are aliased by live
        // verletPoints - they must NEVER be reset, reallocated or
        // resized afterwards, only mutated in place through
        // min()/max() (in-place mutation contract, verletPoint.H)
        List<std::shared_ptr<boundBox>> coverBoxes_;
        labelListList coverPartition_;
        bool coverActive_;

        // watch list of raw box pointers from the previous
        // getBBoxes() call - catches aliasing-contract violation
        List<const boundBox*> coverBoxWatch_;

        // AABB of a subset of triangles (indices into triSurf_)
        boundBox triSetBBox(const labelList& tris) const;

        // box volume with empty (2D) directions masked to 1;
        // boundBox::mag() is the diagonal, NOT the volume
        scalar coverBoxVol(const boundBox& bBox) const;

        // recursive median split of a triangle set along its widest
        // axis, stopped by the stall criterion or the box cap; the
        // budget bounds the leaves this subtree may emit, so the
        // cap is honored by construction
        void coverSplit
        (
            const labelList& tris,
            scalar stallCoeff,
            label& nLeaves,
            label budget,
            List<labelList>& leaves
        ) const;

        // re-derive cover box bounds from the partition and the
        // current (synchronized) surface points - the single
        // refresh mechanism for every motion class
        void refreshCoverBounds();

        // ascii legacy-vtk cover writer (debug machinery)
        // writes:
        // <bodyName>_boxes.vtk (hexahedron per box) and
        // <bodyName>_partition.vtk (triangles colored by box);
        // Note (MI): after the createCover project is finished and
        //            debugged, this function can be removed
        void writeCoverVtk
        (
            const word& bodyName,
            const List<std::shared_ptr<boundBox>>& boxes,
            const labelListList& partition,
            const triSurface& surf,
            const boundBox& fullAABB
        ) const;

        void getIntersectionPoints
        (
            const label index,
            const treeBoundBox& cubeBb,
            DynamicPointList& intersectionPoints
        );

        void intersectBb
        (
            const point& p0,
            const point& p1,
            const point& p2,
            const treeBoundBox& cubeBb,
            DynamicPointList& intersectionPoints
        );

        public:

        // Constructors ------------------------------------------------
        stlBased
        (
            const  fvMesh&   mesh,
            const contactType cType,
            word      stlPath
        );

        // Destructors -------------------------------------------------
        virtual ~stlBased(){}

        virtual std::shared_ptr<geomModel> getCopy() = 0;

        vector addModelReturnRandomPosition(
            const bool allActiveCellsInMesh,
            const boundBox  cellZoneBounds,
            Random&         randGen
        ) final override;

        void bodyMovePoints(vector translVec) final override;

        void bodyScalePoints(scalar scaleFac) final override;

        void bodyRotatePoints
        (
            scalar rotAngle,
            vector axisOfRot
        ) final override;

        void resetBody() final override
        {
            triSurfaceMesh tmpSurf(IOobject
                (
                    stlPath_,
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            bodySurfMesh_.movePoints(tmpSurf.points());
            triSurf_.reset(new triSurface(bodySurfMesh_));
            computeVolumeCoM(triSurf_());
        }

        vector getCoM() final override
        {
            return CoM_;
        }

        boundBox getBounds() final override
        {
            boundBox bounds(bodySurfMesh_.points(), false);
            return bounds;
        }

        // multi-box broad-phase bound; refreshes box bounds from
        // the partition + current points on every call (the
        // verletList::update side-effect refresh path)
        List<std::shared_ptr<boundBox>> getBBoxes() final override;

        // build the cover partition + boxes once; returns true on
        // success, false leaves the legacy single-AABB behavior;
        // writeCover (default off) writes the vtk visualization,
        // bodyName is used for the output file names
        bool computeStaticCover
        (
            label maxBoxes,
            scalar stallCoeff,
            bool writeCover = false,
            word bodyName = word::null
        ) final override;

        void synchronPos(label owner = -1) final override;

        void writeBodyInfo(label bodyId, IOdictionary& dict, word path) final override
        {
            triSurface triToRet(bodySurfMesh_);
            triToRet.write(path+"/stlFiles/"+name(bodyId) + ".stl",".stl");
        }

        boolList pointInside(pointField pointI) final override
        {
            boolList returnList = triSurfSearch_().calcInside( pointI);
            return returnList;
        }

        bool pointInside(point pointI) final override
        {
            pointField points(1,pointI);
            boolList returnList = triSurfSearch_().calcInside( points);
            return returnList[0];
        }

        volumeType getVolumeType(subVolume& sv, bool cIb)  final override;

        void getClosestPointAndNormal
        (
            const point& startPoint,
            const vector& span,
            point& closestPoint,
            vector& normal
        ) final override;

        bool limitFinalSubVolume
        (
            const subVolume& sv,
            bool cIb,
            boundBox& limBBox
        ) final override;

        bool getLeafSubVolumePlane
        (
            subVolume& sv,
            bool cIb,
            point& p,
            vector& n
        ) final override;

        bool getBoxSurfacePlane
        (
            const boundBox& leaf,
            point& p,
            vector& n
        ) final override;

        private:

        // shared dominant-facet tangent-plane helper: picks the facet
        // with the largest triangle-bbox/leaf-box overlap from
        // shapesIn and resolves the inward orientation with
        // pointInside probes either side of the plane
        bool getShapesSurfacePlane
        (
            const labelList& shapesIn,
            const boundBox& leaf,
            point& p,
            vector& n
        );

        public:

        void initVars(pointField surfPoints, scalar intSpan, bool sdBasedLambda)
        {
            bodySurfMesh_.movePoints(surfPoints);
            triSurf_.reset(new triSurface(bodySurfMesh_));
            triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
            historyPoints_ = surfPoints;
            intSpan_ = intSpan;
            sdBasedLambda_ = sdBasedLambda;
            computeVolumeCoM(triSurf_());
        }

        // volume centroid of a closed surface (divergence theorem)
        // Note (MI): based on standard OpenFOAM function
        void computeVolumeCoM(const triSurface& surf)
        {
            scalar m(0);
            tensor J(tensor::zero);
            momentOfInertia::massPropertiesSolid
            (
                surf,
                1.0,
                m,
                CoM_,
                J
            );

            if (m < SMALL)                                              // fallback for singular cases
            {
                WarningInFunction
                    << "stl surface " << stlPath_
                    << " yields non-positive volume (" << m
                    << "): open or inside-out. Falling back to the"
                    << " vertex average - rotation bookkeeping may"
                    << " be inaccurate" << endl;
                const pointField& bSMeshPts = bodySurfMesh_.points();
                CoM_    = vector::zero;
                forAll(bSMeshPts,point)
                {
                    CoM_ += bSMeshPts[point];
                }
                CoM_/= bSMeshPts.size();
                // Note (MI): these kinds of loops can be replaced with
                // CoM_ = sum(bSMeshPts)/bSMeshPts.size();
                // -- and this happens on multiple occasions in the codebase
            }
        }

        void setCoM() final override
        {
            computeVolumeCoM(triSurf_());
        }

        pointField getBodyPoints() final override
        {
            return bodySurfMesh_.points();
        }

        void setBodyPosition(pointField pos) final override;
    };
}
#endif
