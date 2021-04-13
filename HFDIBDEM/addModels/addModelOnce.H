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

#ifndef addModelOnce_H
#define addModelOnce_H

#include "addModel.H"

namespace Foam
{    
    class addModelOnce
    : public addModel
    {
    
        private:
        const dictionary&       addModelDict_;
        word                    addMode_;
        const word              stlName_;
        bool                    bodyAdded_;
        
        const  dynamicFvMesh&   mesh_;
        
        protected:
        
        
        public:
        
        // Constructors ------------------------------------------------
        addModelOnce
        (
            const dictionary& addModelDict,
            const word        stlName,
            const Foam::dynamicFvMesh& mesh
        );
        
        // Destructors -------------------------------------------------
        ~addModelOnce();
        
        // Public subclasses -------------------------------------------
        
        // Public variables --------------------------------------------
        
        // Public methods ---------------------------------------------- 
        //~ void init();
        
        bool shouldAddBody(const volScalarField& body){return not bodyAdded_;}
        
        triSurface addBody(const volScalarField& body);
        
        bool canAddBody
        (
            const volScalarField& body,
            const triSurfaceMesh& bodySurfMesh
        );
        
        bool getBodyAdded(){return bodyAdded_;}
        
                       
    };
}

#endif
