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
            const  dynamicFvMesh&   mesh,
            word      stlPath,
            scalar  thrSurf = 1,
            Vector<label> geometricD = {1,1,1}
        )
        :
        stlBased
        (
            mesh,
            arbShape,
            stlPath,
            thrSurf,
            geometricD
        ){}
        
        // Destructors -------------------------------------------------
        ~nonConvexBody(){}
        
        // Public subclasses -------------------------------------------
        
        // Public variables --------------------------------------------
        
        // Public methods ---------------------------------------------- 
        //~ void init();
        
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
        
        labelList getBBoxCellsByOctTree
        (
            label cellToCheck,
            bool& insideBB,
            vector& bBoxMin,
            vector& bBoxMax,
            List<DynamicLabelList>& bBoxCells,
            Field<label>& octreeField
        ); 
        
        geomModel* getGeomModel() final override
        {
            nonConvexBody* convexBodyModel = new nonConvexBody(mesh_,stlPath_,thrSurf_,geometricD_);
            stlBased* stlBasedModel = convexBodyModel;
            stlBasedModel->initVars(bodySurfMesh_.points(), thrSurf_, geometricD_, intSpan_, sdBasedLambda_);
            geomModel* bodyGeomModel = stlBasedModel;
            return bodyGeomModel;
        }
    };
}

#endif
