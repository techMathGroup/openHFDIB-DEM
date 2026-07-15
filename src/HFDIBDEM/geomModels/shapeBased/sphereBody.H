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
    sphereBody.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef sphereBody_H
#define sphereBody_H

#include "shapeBased.H"
#include "OFstream.H"

#include <assert.h>

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
        scalar sphereM_;

        vector startPosition_;
        scalar startRadius_;

        vector historyPosition_;
        scalar historyRadius_;

        public:

        // Constructors ------------------------------------------------
        sphereBody
        (
            const  fvMesh&   mesh,
            vector startPos,
            scalar radius,
            scalar  thrSurf = 1
        )
        :
        shapeBased
        (
            mesh,
            sphere,
            thrSurf
        ),
        position_(startPos),
        radius_(radius),
        sphereM_((4*Foam::constant::mathematical::pi*pow(radius,3))/3)
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

        void createImmersedBody
        (
            volScalarField& body,
            Field<label>& octreeField,
            List<labelList>& cellPoints
        ) final override;

        label getCellInBody
        (
            Field<label>& octreeField
        );

        std::shared_ptr<geomModel> getCopy() final override
        {
            std::shared_ptr<sphereBody> sphereBodyModel =
                std::make_shared<sphereBody>(mesh_,position_,radius_,thrSurf_);

            sphereBodyModel->initVars(intSpan_,sdBasedLambda_);
            sphereBodyModel->setRhoS(rhoS_);
            sphereBodyModel->setM0(M0_);
            sphereBodyModel->getSurfaceCellList() = surfCells_;
            sphereBodyModel->getInternalCellList() = intCells_;

            return sphereBodyModel;
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

        void synchronPos(label owner = -1) final override;

        void writeBodyInfo(label bodyId, IOdictionary& dict, word path) final override
        {
            dictionary sphereDict("sphere");
            sphereDict.add("position", position_);
            sphereDict.add("radius", radius_);

            dict.add("sphere",sphereDict);
        }

        boolList pointInside(pointField pointI) final override;


        bool pointInside(point pointI) final override;

        void getClosestPointAndNormal
        (
            const point& startPoint,
            const vector& span,
            point& closestPoint,
            vector& normal
        ) final override
        {
            if(startPoint == position_)
            {
                closestPoint = startPoint;
                normal = vector::zero;
            }
            else
            {
                normal = startPoint-position_;
                normal /= mag(normal);
                closestPoint = position_ + normal * radius_;
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

        scalar& getM0() final override
        {
            if(case3D)
            {
                return sphereM_;
            }
            else
            {
                return M0_;
            }
        }

        void setRhoS(dimensionedScalar rhoS)
        {
            rhoS_ = rhoS;
            sphereM_ = ((4*Foam::constant::mathematical::pi*pow(radius_,3))/3)*rhoS_.value();
        }

        void setCoM() final override
        {
            CoM_ = position_;
        }

        void setBodyPosition(pointField pos) final override
        {
            assert(pos.size() == 1 && "Only one point is allowed for sphere body");
            position_ = pos[0];
        }

        pointField getBodyPoints() final override
        {
            pointField points(1,position_);
            return points;
        }
    };
}

#endif
