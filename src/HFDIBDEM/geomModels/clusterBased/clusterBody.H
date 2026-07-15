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

        protected:

        std::vector<std::shared_ptr<geomModel>> ibGeomModelList;

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
            cluster,
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
            List<labelList>& cellPoints
        ) final override;

        virtual std::shared_ptr<geomModel> getCopy() = 0;

        // Connect int/surf cells with correct center of mass
        // for fluid-solid interaction
        virtual void getReferencedLists(
            List<DynamicLabelList>& intLists,
            List<DynamicLabelList>& surfLists,
            DynamicVectorList& referenceCoM
        ) = 0;

        void updateSurfList();

        virtual List<DynamicLabelList>& getSurfaceCellList() final override
        {
            updateSurfList();
            return surfCells_;
        };

        void updateIntList();

        virtual List<DynamicLabelList>& getInternalCellList() final override
        {
            updateIntList();
            return intCells_;
        };

        void calculateGeometricalProperties(
            volScalarField& body
        ) final override;

        void calculateGeometricalPropertiesParallel(
            volScalarField& body
        ) final override;

        void setMassAndInertia();

        virtual void bodyMovePoints(vector translVec) = 0;

        virtual void bodyScalePoints(scalar scaleFac) = 0;

        virtual void bodyRotatePoints
        (
            scalar rotAngle,
            vector axisOfRot
        ) = 0;

        virtual vector getCoM() = 0;

        boundBox getBounds() final override;

        void synchronPos(label owner = -1) final override;

        boolList pointInside(pointField pointI) final override;

        bool pointInside(point pointI) final override;

        // Get characteristic particle diameter
        virtual scalar getDC() = 0;

        void getClosestPointAndNormal
        (
            const point& startPoint,
            const vector& span,
            point& closestPoint,
            vector& normal
        ) final override;

        // Get processor owner of the particle
        virtual label getOwner() = 0;

        // Get initial mass of the cluster
        virtual scalar& getM0() = 0;

        // Get vector from correct center of mass to toPoint
        virtual vector getLVec(const point& toPoint) = 0;

        void resetBody(volScalarField& body) final override;

        bool isCluster() const final override
        {
            return true;
        }

        std::vector<std::shared_ptr<geomModel>>& getClusterBodies()
        {
            return ibGeomModelList;
        }

        List<std::shared_ptr<boundBox>> getBBoxes() final override;

        virtual void writeBodyInfo(label bodyId, IOdictionary& dict, word path)
        {
            InfoH << iB_Info << "-- geomModelMessage-- "
                << "custer body writeBodyInfo notIplemented"<<endl;
        }
    };
}

#endif
