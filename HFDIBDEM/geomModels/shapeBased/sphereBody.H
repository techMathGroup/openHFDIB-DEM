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
    addModels.C

Contributors
    Martin Isoz (2019-*), Martin Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef sphereBody_H
#define sphereBody_H

#include "shapeBased.H"
#include "OFstream.H"

namespace Foam
{
    class sphereBody
    : public shapeBased
    {

        private:

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        vector position_;
        scalar radius_;

        vector startPosition_;
        scalar startRadius_;

        vector historyPosition_;
        scalar historyRadius_;

        protected:


        public:

        // Constructors ------------------------------------------------
        sphereBody
        (
            const  dynamicFvMesh&   mesh,
            vector startPos,
            scalar radius,
            scalar  thrSurf = 1,
            Vector<label> geometricD = {1,1,1}
        )
        :
        shapeBased
        (
            mesh,
            sphere,
            thrSurf,
            geometricD
        ),
        position_(startPos),
        radius_(radius)
        {
            startPosition_ = position_;
            startRadius_ = radius_;

            historyPosition_ = position_;
            historyRadius_ = radius_;
        }

        // Destructors -------------------------------------------------
        ~sphereBody(){}

        // Public subclasses -------------------------------------------

        // Public variables --------------------------------------------

        // Public methods ----------------------------------------------

        bool canAddBody
        (
            const volScalarField& body
        ) final override;

        void createImmersedBody
        (
            volScalarField& body,
            Field<label>& octreeField,
            List<DynamicLabelList>& surfCells,
            List<DynamicLabelList>& intCells,
            List<pointField>& cellPoints
        ) final override;

        labelList createImmersedBodyByOctTree
        (
            label cellToCheck,
            bool& insideIB,
            bool& centerInside,
            boolList& vertexesInside,
            volScalarField& body,
            List<DynamicLabelList>& surfCells,
            List<DynamicLabelList>& intCells
        );

        geomModel* getGeomModel() final override
        {
            sphereBody* sphereBodyModel = new sphereBody(mesh_,position_,radius_,thrSurf_,geometricD_);
            sphereBodyModel->initVars(intSpan_,sdBasedLambda_);
            geomModel* bodyGeomModel = sphereBodyModel;

            return bodyGeomModel;
        }

        void bodyMovePoints(vector translVec) final override
        {position_ += translVec;}

        void bodyScalePoints(scalar scaleFac) final override
        {radius_ *= scaleFac;}

        void bodyRotatePoints
        (
            scalar rotAngle,
            vector axisOfRot
        ) final override
        {}

        void resetBody() final override
        {
            position_ = startPosition_;
            radius_ = startRadius_;
        }

        vector getCoM() final override
        {return position_;}

        boundBox getBounds() final override
        {
            vector one(1,1,1);
            one *= radius_;
            boundBox bounds(position_-one,position_+one);
            return bounds;
        }

        void recordHistory() final override
        {
            historyPosition_ = position_;
            historyRadius_ = radius_;
        }

        void returnHistory() final override
        {
            position_ = historyPosition_;
            radius_ = historyRadius_;
        }

        void synchronPos() final override;

        void writeBodyInfo(label bodyId, IOdictionary& dict, word path) final override
        {
            dictionary sphereDict("sphere");
            sphereDict.add("position", position_);
            sphereDict.add("radius", radius_);

            dict.add("sphere",sphereDict);
        }

        boolList pointInside(pointField pointI) final override;


        bool pointInside(point pointI) final override;

        vector getClosestNormal(point startPoint, vector span) final override
        {
            if(startPoint == position_)
                return vector::zero;
            else
            {
                vector normal(startPoint-position_);
                normal /= mag(normal);
                return normal;
            }
        }

        void initVars(scalar intSpan, bool sdBasedLambda)
        {
            intSpan_ = intSpan;
            sdBasedLambda_ = sdBasedLambda;
        }

        pointField sampleSurfacePoints() final override;

        scalar getDC() final override
        {
            return 2*radius_;
        }
    };
}

#endif
