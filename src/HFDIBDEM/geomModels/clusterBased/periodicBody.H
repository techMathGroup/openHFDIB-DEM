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
    periodicBody.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef periodicBody_H
#define periodicBody_H

#include "clusterBody.H"
#include "OFstream.H"

namespace Foam
{
    class periodicBody
    : public clusterBody
    {

        private:

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        protected:


        public:

        // Constructors ------------------------------------------------
        periodicBody
        (
            const  fvMesh&   mesh,
            scalar  thrSurf = 1
        )
        :
        clusterBody
        (
            mesh,
            thrSurf
        )
        {}

        // Destructors -------------------------------------------------
        ~periodicBody(){}

        // Public subclasses -------------------------------------------

        // Public variables --------------------------------------------

        // Public methods ----------------------------------------------

        std::shared_ptr<geomModel> getCopy() final override
        {
            std::shared_ptr<periodicBody> periodicGeomModel = std::make_shared<periodicBody>
            (
                mesh_,
                thrSurf_
            );

            periodicGeomModel->ibGeomModelList = ibGeomModelList;
            return periodicGeomModel;
        }

        void getReferencedLists(
            List<DynamicLabelList>& intLists,
            List<DynamicLabelList>& surfLists,
            DynamicVectorList& referenceCoM
        ) final override;

        void bodyMovePoints(vector translVec) final override;

        void bodyScalePoints(scalar scaleFac) final override;

        void bodyRotatePoints
        (
            scalar rotAngle,
            vector axisOfRot
        ) final override;

        vector getCoM() final override;

        scalar getDC() final override;

        // Set processor owner of the particle
        void setOwner() final override;

        // Get processor owner of the particle
        label getOwner() final override;

        // Get initial mass of the cluster
        scalar& getM0() final override;

        vector getLVec(const point& toPoint) final override;

        void addBodyToCluster(std::shared_ptr<geomModel> bodyToAdd)
        {
            ibGeomModelList.push_back(std::move(bodyToAdd));
        }

        bool shouldBeUnclustered();

        std::shared_ptr<geomModel> getRemGeomModel();

        void writeBodyInfo(label bodyId, IOdictionary& dict, word path) final override
        {
            ibGeomModelList[0]->writeBodyInfo(bodyId, dict, path);
        }
    };
}

#endif
