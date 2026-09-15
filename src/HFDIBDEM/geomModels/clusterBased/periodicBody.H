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
    class for clustering of immersedBodies

SourceFiles
    periodicBody.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
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
