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
        }

        vector getCoM() final override
        {
            pointField bSMeshPts = bodySurfMesh_.points();
            // get its center of mass
            vector CoM(vector::zero);
            forAll(bSMeshPts,point)
            {
                CoM += bSMeshPts[point];
            }
            CoM/= bSMeshPts.size();
            return CoM;
        }

        boundBox getBounds() final override
        {
            boundBox bounds(bodySurfMesh_.points(), false);
            return bounds;
        }

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

        void initVars(pointField surfPoints, scalar intSpan, bool sdBasedLambda)
        {
            bodySurfMesh_.movePoints(surfPoints);
            triSurf_.reset(new triSurface(bodySurfMesh_));
            triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
            historyPoints_ = surfPoints;
            intSpan_ = intSpan;
            sdBasedLambda_ = sdBasedLambda;
        }

        void setCoM() final override
        {
            pointField bSMeshPts = bodySurfMesh_.points();
            CoM_    = vector::zero;
            forAll(bSMeshPts,point)
            {
                CoM_ += bSMeshPts[point];
            }
            CoM_/= bSMeshPts.size();
        }

        pointField getBodyPoints() final override
        {
            return bodySurfMesh_.points();
        }

        void setBodyPosition(pointField pos) final override;
    };
}
#endif
