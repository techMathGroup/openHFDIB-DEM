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
        //~ void init();

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
            std::shared_ptr<nonConvexBody> convexBodyModel =
                std::make_shared<nonConvexBody>(mesh_,stlPath_,thrSurf_);

            convexBodyModel->initVars(bodySurfMesh_.points(), thrSurf_, intSpan_, sdBasedLambda_);
            convexBodyModel->setRhoS(rhoS_);
            convexBodyModel->setM0(M0_);
            convexBodyModel->getSurfaceCellList() = surfCells_;
            convexBodyModel->getInternalCellList() = intCells_;

            return convexBodyModel;
        }
    };
}

#endif
