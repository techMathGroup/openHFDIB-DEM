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

Description
    class for addition and removal of immersedBodies

SourceFiles
    stlBased.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
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
            word      stlPath,
            scalar  thrSurf
        );

        // Destructors -------------------------------------------------
        virtual ~stlBased(){}

        virtual void createImmersedBody
        (
            volScalarField& body,
            Field<label>& octreeField,
            List<labelList>& cellPoints
        )
        {
            InfoH << iB_Info << "-- geomModelMessage-- "
                << "notIplemented"<<endl;
        }

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
        ) final override
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
                InfoH << basic_Info << "Missing the closest point!" << endl;
                normal = startPoint - getCoM();
                closestPoint = getCoM();
            }
        }

        bool limitFinalSubVolume
        (
            const subVolume& sv,
            bool cIb,
            boundBox& limBBox
        ) final override;

        void initVars(pointField surfPoints, scalar thrSurf, scalar intSpan, bool sdBasedLambda)
        {
            bodySurfMesh_.movePoints(surfPoints);
            triSurf_.reset(new triSurface(bodySurfMesh_));
            triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
            historyPoints_ = surfPoints;
            thrSurf_ = thrSurf;
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
