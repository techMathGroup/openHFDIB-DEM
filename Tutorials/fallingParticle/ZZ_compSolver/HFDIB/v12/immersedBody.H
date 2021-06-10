/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________
                       | | | ||  ___|  _  \_   _| ___ \     H ybrid
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /     F ictitious
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \     D omain
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /     I mmersed
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/      B oundary
      | |
      |_|
-------------------------------------------------------------------------------
License

    openHFDIB is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

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
    class for immersed bodies representation.

SourceFiles
    immersedBodies.C

Contributors
    Federico Municchi (2016)
\*---------------------------------------------------------------------------*/

#ifndef immersedBody_H
#define immersedBody_H

#include "dictionary.H"
#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedTensor.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "fvMatricesFwd.H"
//~ #include "fvMesh.H"
#include "dynamicFvMesh.H"
#include "triSurfaceMesh.H"

namespace Foam
{
    enum
    {
                STATICBODY, //Body does not move
         TRANSROTATINGBODY, //Trans-rotational imposed motion
             FLUIDCOUPLING  //Fluid driven motion
    };

    class immersedBody
    {
    
        private:
        
        // Private data types
        
        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;
        
        //first time update?
        bool    isFirstUpdate_;
        
        //is the particle in contact with other particles
        bool    isInPrtContact_;
        
        //is the particle in contact with walls
        bool    isInWallContact_;
        
        //Dictionary
        dictionary&                                  immersedDict_;
        
        //Mesh
        const  dynamicFvMesh&                               mesh_;
        
        //Dictionary
        dictionary&                        transportProperties_;
        
        
        //Body mass
        scalar M_;
        
        //Location of center of mass
        vector CoM_;
        
        //Rotation axis
        vector Axis_;
        
        //Angular velocity
        scalar omega_;
        
        //Body velocity
        vector Vel_;
        
        //Moment of intertia tensor
        symmTensor  I_;
        
        //Forces acting on the particle
        vector  F_;
        
        //Total torque acting on the particle
        vector  T_;
        
        //spring stiffness for normal contact force
        scalar  kN_;
        
        //dumping coefficient for normal contact force
        scalar gammaN_;
        
        //threshold to distinguish surface cells [thrSurf_,1-thrSurf_]
        scalar  thrSurf_;
        
        //integer to identify current body
        label   bodyId_;
        
        //Operation to perform
        int     bodyOperation_;
        
        //base geometry
        triSurfaceMesh      bodySurfMesh_;
        pointField          bodyPoints_;
        
        //surface cells
        DynamicLabelList    surfCells_;
        
        //internal cells
        DynamicLabelList    intCells_;
        
        //cells in contact with other particles
        DynamicLabelList    prtContactCells_;
        
        //cells in contact with walls
        DynamicLabelList    wallContactCells_;
        
        //interpolation points, point[0] is a surface point
        DynamicList<DynamicPointList>  interpolationPoints_;
        
        //interpolation cells
        DynamicList<DynamicLabelList>   interpolationCells_;
        
        
        //Create immersed body info
        void  createImmersedBody( volScalarField& body );
        
        //Update immersed body info
        void  updateImmersedBody( volScalarField& body,
                                 volVectorField& f
                               );
        
        //Create interpolation points
        void  calculateInterpolationPoints(volScalarField& body);
        
        //Evaluate center of mass and inertia tensor
        vector addToMAndI
        (
            volScalarField& body,
            DynamicLabelList& labelCellLst,
            vector tmpCom
        );
        void calculateGeometricalProperties(volScalarField& body);
        
        //Rotate immersed body   
        void moveImmersedBody();
        
        //Reset body field
        void resetBody(volScalarField& body);
        
        //Refine body field using MC-like integration
        //~ void refineBody(volScalarField& body,
                       //~ triSurfaceSearch* ibTriSurfSearch,
                       //~ const pointField* pp
                   //~ );
        
        //Update Fluid-Solid fluid
        void updateCoupling
        ( 
            volScalarField& body,
            volVectorField& f
        );
        
        //update body movement based on F_ and T_
        void updateMovement();
        
        
        public:
        
        immersedBody
        (
            word fileName,
            const Foam::dynamicFvMesh& mesh,
            dictionary& HFDIBDict,
            dictionary& transportProperties,
            label bodyId
        );
        
        ~immersedBody();
        
        //Getters
        
        const bool& checkPrtContact() const
        {
            return isInPrtContact_;
        }
        
        const bool& checkWallContact() const
        {
            return isInWallContact_;
        }
        
        const label& getBodyId() const
        {
            return bodyId_;
        }
        
        const scalar& getM()    const
        {
            return M_;
        }
        
        const vector& getVel()  const
        {
            return Vel_;
        }
        
        const vector& getCoM()  const
        {
            return CoM_;
        }
        
        const scalar& getOmega()    const
        {
            return omega_;
        }
        
        const vector& getAxis()     const
        {
            return Axis_;
        }
        
        const vector& getF()    const
        {
            return F_;
        }
        
        const vector& getT()    const
        {
            return T_;
        }
        
        const symmTensor& getI()    const
        {
            return I_;
        }
        
        const scalar& getKN()   const
        {
            return kN_;
        }
        
        const scalar& getGammaN()   const
        {
            return gammaN_;
        }
        
        const DynamicLabelList& getSurfaceCellList()        const 
        {
            return surfCells_;
        };
        
        const DynamicLabelList& getPrtContactCellList()     const 
        {
            return prtContactCells_;
        };
        
        const DynamicLabelList& getWallContactCellList()    const 
        {
            return wallContactCells_;
        };
        
        const List<DynamicPointList>& getInterpolationPoints() const
        {
            return interpolationPoints_;
        };
        
        const List<DynamicLabelList>& getInterpolationCells() const
        {
            return interpolationCells_;
        };
        
        //Update body field
        void updateBodyField
        ( 
            volScalarField& body,
            volVectorField& f
        );
        
        //Update imposed vector field
        void updateVectorField
        (
            volVectorField& VS,
            word VName,
            volScalarField& body
        );
        
        // reset IB (it will be recreated upon next call to updateBodyField)
        void resetIB(volScalarField& body);
        
        // detect contact with other particles
        void detectContact(volScalarField& body);
        
        // update force and torque acting on the body
        void updateFAndT
        (
            vector FIn,
            vector TIn
        )
        {
            F_+=FIn;
            T_+=TIn;
        }
        
        // move the body post-contact
        void postContactBodyUpdate(volScalarField& body);
        
        // function to compute local particle radii w.r.t. CoM_
        DynamicScalarList getLocPartRad(DynamicLabelList& cellsOfInt);
    };
}
#endif
