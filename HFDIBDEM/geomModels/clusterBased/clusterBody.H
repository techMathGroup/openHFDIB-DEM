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
    class for clustering of immersedBodies

SourceFiles
    clusterBody.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef clusterBody_H
#define clusterBody_H

#include "geomModel.H"
#include "OFstream.H"

namespace Foam
{
    class clusterBody
    : public geomModel
    {

        private:

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        PtrList<geomModel> ibGeomModelList;

        protected:


        public:

        // Constructors ------------------------------------------------
        clusterBody
        (
            const  fvMesh&   mesh,
            scalar  thrSurf = 1
        )
        :
        geomModel
        (
            mesh,
            arbShape,
            thrSurf
        )
        {}

        // Destructors -------------------------------------------------
        ~clusterBody(){}

        // Public subclasses -------------------------------------------

        // Public variables --------------------------------------------

        // Public methods ----------------------------------------------

        void createImmersedBody
        (
            volScalarField& body,
            Field<label>& octreeField,
            List<pointField>& cellPoints
        ) final override;

        geomModel* getGeomModel() final override
        {
            return this;
        }

        void getReferencedLists(
            List<DynamicLabelList>& intLists,
            List<DynamicLabelList>& surfLists,
            DynamicVectorList& referenceCoM
        ) final override;

        void calculateGeometricalProperties(
            volScalarField& body
        ) final override;

        void bodyMovePoints(vector translVec) final override;

        void bodyScalePoints(scalar scaleFac) final override;

        void bodyRotatePoints
        (
            scalar rotAngle,
            vector axisOfRot
        ) final override;


        boundBox getBounds() final override;

        void synchronPos() final override;

        boolList pointInside(pointField pointI) final override;


        bool pointInside(point pointI) final override;

        scalar getDC() final override;

        void getClosestPointAndNormal
        (
            const point& startPoint,
            const vector& span,
            point& closestPoint,
            vector& normal
        ) final override;

        scalar& getM0() final override;

        bool isCluster() const final override
        {
            return true;
        }

        void addBodyToCluster(autoPtr<geomModel> bodyToAdd)
        {
            ibGeomModelList.append(bodyToAdd);
        }
    };
}

#endif
