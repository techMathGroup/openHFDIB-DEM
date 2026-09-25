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
    shapeBased.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef shapeBased_H
#define shapeBased_H

#include "geomModel.H"
#include "meshSearch.H"

namespace Foam
{
    class shapeBased
    : public geomModel
    {

        private:

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        protected:

        public:

        // Constructors ------------------------------------------------
        shapeBased
        (
            const  fvMesh&   mesh,
            const contactType cType,
            scalar  thrSurf
        );

        // Destructors -------------------------------------------------
        virtual ~shapeBased(){}

        virtual void createImmersedBody
        (
            volScalarField& body,
            Field<label>& octreeField,
            List<labelList>& cellPoints
        )
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual std::shared_ptr<geomModel> getCopy() = 0;

        vector addModelReturnRandomPosition(
            const bool allActiveCellsInMesh,
            const boundBox  cellZoneBounds,
            Random&         randGen
        ) final override;

        virtual void bodyMovePoints(vector translVec)
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual void bodyScalePoints(scalar scaleFac)
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual void bodyRotatePoints
        (
            scalar rotAngle,
            vector axisOfRot
        )
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual void resetBody()
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual vector getCoM()
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
            return vector::zero;
        }

        virtual boundBox getBounds()
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
            boundBox bound;
            return bound;
        }

        virtual void synchronPos(label owner)
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual void writeBodyInfo(label bodyId, IOdictionary& dict, word path)
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual boolList pointInside(pointField pointI)
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
            return {false};
        }


        virtual bool pointInside(point pointI)
        {
            InfoH << iB_Info << "-- geomModelMessage-- "
                << "notIplemented"<<endl;
            return false;
        }

        virtual void getClosestPointAndNormal
        (
            const point& startPoint,
            const vector& span,
            point& closestPoint,
            vector& normal
        )
        {
            InfoH << iB_Info << "-- geomModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual pointField sampleSurfacePoints()
        {
            InfoH << iB_Info << "-- shapeModelMessage-- "
                << "notIplemented"<<endl;
            pointField returnField;
            return returnField;
        }

        virtual void setCoM()
        {
            InfoH << iB_Info << "-- geomModelMessage-- "
                << "notIplemented"<<endl;
        }
    };
}
#endif
