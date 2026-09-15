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
    convexBody.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef convexBody_H
#define convexBody_H

#include "stlBased.H"

namespace Foam
{
    class convexBody
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
        convexBody
        (
            const  fvMesh&   mesh,
            word      stlPath,
            scalar  thrSurf = 1
        )
        :
        stlBased
        (
            mesh,
            convex,
            stlPath,
            thrSurf
        ){}

        // Destructors -------------------------------------------------
        ~convexBody(){}

        // Public subclasses -------------------------------------------

        // Public variables --------------------------------------------

        // Public methods ----------------------------------------------
        //~ void init();

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
            std::shared_ptr<convexBody> convexBodyModel =
                std::make_shared<convexBody>(mesh_,stlPath_,thrSurf_);

            convexBodyModel->initVars(bodySurfMesh_.points(), thrSurf_, intSpan_, sdBasedLambda_);
            convexBodyModel->setRhoS(rhoS_);
            convexBodyModel->setM0(M0_);
            convexBodyModel->getSurfaceCellList() = surfCells_;
            convexBodyModel->getHaloCellList() = haloCells_;
            convexBodyModel->getInternalCellList() = intCells_;

            return convexBodyModel;
        }
    };
}

#endif
