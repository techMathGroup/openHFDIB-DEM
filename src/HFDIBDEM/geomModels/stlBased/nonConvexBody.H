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
    nonConvexBody.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef nonConvexBody_H
#define nonConvexBody_H

#include "stlBased.H"

namespace Foam
{
    class nonConvexBody
    : public stlBased
    {

        private:

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        protected:


        public:

        // Constructors ------------------------------------------------
        nonConvexBody
        (
            const  fvMesh&   mesh,
            word      stlPath,
            scalar  thrSurf = 1
        )
        :
        stlBased
        (
            mesh,
            nonConvex,
            stlPath,
            thrSurf
        ){}

        // Destructors -------------------------------------------------
        ~nonConvexBody(){}

        // Public subclasses -------------------------------------------

        // Public variables --------------------------------------------

        // Public methods ----------------------------------------------

        void createImmersedBody
        (
            volScalarField& body,
            Field<label>& octreeField,
            List<labelList>& cellPoints
        ) final override;

        labelList getBBoxCellsByOctTree
        (
            label cellToCheck,
            bool& insideBB,
            vector& bBoxMin,
            vector& bBoxMax,
            List<DynamicLabelList>& bBoxCells,
            Field<label>& octreeField
        );

        std::shared_ptr<geomModel> getCopy() final override
        {
            std::shared_ptr<nonConvexBody> nonConvexBodyModel =
                std::make_shared<nonConvexBody>(mesh_,stlPath_,thrSurf_);

            nonConvexBodyModel->initVars(bodySurfMesh_.points(), thrSurf_, intSpan_, sdBasedLambda_);
            nonConvexBodyModel->setRhoS(rhoS_);
            nonConvexBodyModel->setM0(M0_);
            nonConvexBodyModel->getSurfaceCellList() = surfCells_;
            nonConvexBodyModel->getHaloCellList() = haloCells_;
            nonConvexBodyModel->getInternalCellList() = intCells_;

            return nonConvexBodyModel;
        }
    };
}

#endif
