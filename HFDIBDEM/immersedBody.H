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
Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Šourek (2019-*)
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

#include "geomModel.H"
#include "stlBased.H"
#include "convexBody.H"
#include "nonConvexBody.H"


namespace Foam
{   
//     // custom field filtration 
//     template<class Type>
//     Foam::tmp<Foam::Field<Type>>
//     filterField
//     (
//         const Field<Type>& field,
//         const labelList& cellIDs
//     ) 
//     {
//         if (isNull(cellIDs))
//         {
//             return field;
//         }
//         else
//         {
//             return tmp<Field<Type>>(new Field<Type>(field, cellIDs));
//         }
//     }

    class immersedBody
    {
    
        private:
        
        // Private data types
        
        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;
        
        // name of the STL to base the body on
        word    bodyName_;
        
        //is the particle in contact with other particles
        bool    isInPrtContact_;
        
        //is the particle in contact with wall
        bool    isInWallContact_;
        
        //is the particle active in the simulation
        bool    isActive_;
        
        //did the body move enough that it needs to be reprojected?
        bool    recomputeProjection_;
        
        //Dictionary
        dictionary&             immersedDict_;
        
        //Mesh
        const  dynamicFvMesh&   mesh_;
        
        //Dictionary
        dictionary&             transportProperties_;
        
        autoPtr<geomModel> geomModel_;
        //Precalculate cellPoints for better performace
        List<pointField>& cellPoints_;
        
        
        //Body mass
        scalar M_;
        scalar M0_;                                                     //initial mass (at addition)
        
        //Location of center of mass
        vector CoM_;
        
        //Rotation axis
        vector Axis_;
        vector AxisOld_;
        
        //Angular velocity
        scalar omega_;
        scalar omegaOld_;
        
        //Body velocity
        vector Vel_;
        vector VelOld_;
        
        //Body linear acceleration
        vector a_;
        
        //Body angular acceleration
        vector alpha_;
        
        //Total Rodrigues angle of the body
        vector totalAngle_;
        
        //Moment of intertia tensor
        symmTensor  I_;
        
        //Forces acting on the particle
        vector  F_;
        
        //Total torque acting on the particle
        vector  T_;
        
        //spring stiffness for normal contact force
        scalar  kN_;
        
        //spring stiffness for tangential contact force
        scalar  kt_;
        
        //damping coefficient for normal contact force
        scalar gammaN_;
        
        //damping coefficient for tangential contact force
        scalar gammat_;
        
        //friction coefficient
        scalar mu_;
        
        //adhesion coefficient for normal contact force
        scalar adhN_;
        
        //adhesion coefficient for normal contact force
        scalar adhEqui_;
        
        //maximal body Courant number based on its linear velocity
        scalar CoNum_;
        
        //mean body Courant number based on its linear velocity
        scalar meanCoNum_;
        
        //fluid density
        dimensionedScalar rhoF_;
        
        //body density
        dimensionedScalar rhoS_;
        
        //characteristic diameter of the body
        scalar  dC_;
        
        //integer to identify current body
        label   bodyId_;
        
        //Auxiliary switch - do I want to update the torques?
        bool    updateTorque_;
        
        //Operation to perform
        int     bodyOperation_;
        
        //fraction of dC_to move by in one DEM loop
        scalar  maxDistInDEMloop_;
        
        //index holder for particle contact list
        List<List<label>> boundIndList_;
        
        //Bounding points of the surface
        point minBoundPoint_;
        point maxBoundPoint_;
        
        //surface cells
        List<DynamicLabelList>  surfCells_;
        
        //internal cells
        List<DynamicLabelList>  intCells_;
        
        //faces in contact with walls
        List<DynamicLabelList>  wallContactFaces_;
        
        // helper list to hold faces that were already checked for wall contact
        DynamicLabelList wallContactFacesHelp_;
        
        //processor num that owns the IB
        label owner_;
        
        //Relaxation factor for velocity
        scalar velRelaxFac_;
        
        //Update immersed body info
        void  updateImmersedBody
        ( 
            volScalarField& body,
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
        
        vector addToM
        (
            volScalarField& body,
            DynamicLabelList& labelCellLst,
            vector tmpCom
        );
        
        void addToI
        (
            volScalarField& body,
            DynamicLabelList& labelCellLst
        );
        
        void calculateGeometricalProperties(volScalarField& body);
        void calculateGeometricalProperties2(volScalarField& body, scalar delatT);
        void updateCoM(volScalarField& body, scalar deltaT);
        void updateI(volScalarField& body);
        
        
        //Update Fluid-Solid fluid
        void updateCoupling
        ( 
            volScalarField& body,
            volVectorField& f
        );
                        
        // history variables for contact correction . They are assigned at the beggining of the main (Pimple) time step
        //Forces acting on the particle
        vector  historyCouplingF_;
        
        //Total torque acting on the particle
        vector  historyCouplingT_;
        
        //Rotation axis
        vector historyAxis_;
        
        //Angular velocity
        scalar historyOmega_;
        
        //Body velocity
        vector historyVel_;
        
        //Body linear acceleration
        vector historya_;
        
        //Body angular acceleration
        vector historyAlpha_;
        
        //Total Rodrigues angle of the body
        vector historyTotalAngle_;
        
        // time step required for a contact resolution (Ib should move maximally by the set of fraction of its dC_)
        scalar      contactTimeRes_;
        
        //Tangential force from last timestep
        //First label -1 indicates wall other for bodyID 
        DynamicList<Tuple2<label,Tuple2<label,vector>>> historyFt_;
        
        // variables for re-using the same body several times
        // Note: all the copies have to have false as re-use
        //~ label   useNTimes_;
        //  -1 -> never repeat
        //   0 -> use indefinitly
        //   N -> use N times
        //~ word    startPos_;                                              //NOT IMPLEMENTED
        //~ scalar  timeBetweenUsage_;
        // approximate time before re-using the particles
        
        // auxiliary octree variables
        Field<label> octreeField_;
        label cellToStartInCreateIB_;
        
        bool    startSynced_;
        
        // just auxiliaries
        tensor  totRotMatrix_;
        bool    sdBasedLambda_;
        scalar  intSpan_;
        scalar  dS_;//how much did the body move since last projection
        scalar  charCellSize_;//characteristic cell size over surface
        label  refineBuffers_;//number of refine buffers around body
        bool  useInterpolation_;//number of refine buffers around body
        label recomputeM0_;//how many times rocompute M0_ because of refinement
        vector geometricD_;//geometric directions
        bool created_;
        bool inContactWithStatic_;
        label timeStepsInContWStatic_;
        label timesToSetStatic_;
        
        void initializeIB();
        
        protected:
        
        
        public:
        
        // Constructors ------------------------------------------------
        // - constructor loading STL from file
        immersedBody
        (
            word fileName,
            const Foam::dynamicFvMesh& mesh,
            dictionary& HFDIBDEMDict,
            dictionary& transportProperties,
            label bodyId,
            label recomputeM0,
            vector geometricD,
            geomModel* bodyGeomModel,
            List<pointField>& cellPoints
        );
        
        // Destructors -------------------------------------------------
        ~immersedBody();
        
        // Public subclasses -------------------------------------------
        //Class to keep interpolation data
        class interpolationInfo{
        public:
            label surfCell_;
            DynamicPointList intPoints_;
            label order_;
            DynamicLabelList intCells_;
            DynamicLabelList procWithIntCells_;
            DynamicVectorList intVec_;
        };
        
        //Class to ask for interpolation data on other processor
        class intVecRequest
        {
            public:
                label requestLabel_;
                label vecLabel_;
                point intPoint_;
                label intCell_;
        };
        
        // Public variables --------------------------------------------
        // Note (MI): why is this here? There is a getter bellow (20201126)
        //            -> catch XXII caused by the necessity to have public subclass
        //interpolation info list
        List<DynamicList<immersedBody::interpolationInfo>>  interpolationInfo_;
        
        //interpolation vec request list
        List<DynamicList<immersedBody::intVecRequest>>      interpolationVecReqs_;
        
        // Public methods ----------------------------------------------
        
        //Rotate immersed body   
        void moveImmersedBody
        (
            scalar deltaT = -1.0
        );
        
        //Reset body field
        void resetBody(volScalarField& body, bool resethistoryFt = true);
        
        //update body movement based on F_ and T_
        void updateMovement();
        void updateMovement
        (
            scalar deltaT
        );
        void updateMovement
        (
            scalar deltaT,
            vector Vel,
            vector Axis,
            scalar omega
        );
        void updateMovement
        (
            scalar deltaT,
            vector Vel,
            vector Axis,
            scalar omega,
            scalar velRelaxFac
        );
        void updateMovement
        (
            vector Vel,
            vector Axis,
            scalar omega
        );
        void updateMovement
        (
            vector Vel,
            vector Axis,
            scalar omega,
            scalar velRelaxFac
        );
        void updateMovementComp
        (
            scalar deltaT,
            vector Vel,
            vector Axis,
            scalar omega,
            scalar velRelaxFac
        );
        
        //Create immersed body info
        void createImmersedBody
        (
            volScalarField& body,
            volScalarField& refineF,
            bool synchCreation = true
        );
        
        void synchCreateImmersedBody
        (
            volScalarField& body,
            volScalarField& refineF
        );

        void constructRefineField
        (
            volScalarField& body,
            volScalarField& refineF,
            DynamicLabelList cellsToIterate,
            DynamicLabelList startLevel 
        );
        
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
            return M0_;
        }
        
        const dimensionedScalar& getRhoS()    const
        {
            return rhoS_;
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
        
        const scalar& getKt()   const
        {
            return kt_;
        }
        
        const scalar& getGammat()   const
        {
            return gammat_;
        }
        
        const scalar& getmu()   const
        {
            return mu_;
        }
        
        const scalar& getadhN()   const
        {
            return adhN_;
        }
        
        const scalar& getadhEqui()   const
        {
            return adhEqui_;
        }
        
        const scalar& getCoNum()   const
        {
            return CoNum_;
        }
        
        const scalar& getMeanCoNum()   const
        {
            return meanCoNum_;
        }
        
        const scalar& getDC()   const
        {
            return dC_;
        }
        
        const scalar& getcontactTimeRes()   const
        {
            return contactTimeRes_;
        }
        
        const int& getbodyOperation()   const
        {
            return bodyOperation_;
        }
        
        void updateContactTimeRes
        (
            scalar newContactTime
        )
        {
            if (newContactTime < contactTimeRes_)
            {
                contactTimeRes_ = newContactTime;
            }
        }
        
        const DynamicList<Tuple2<label,Tuple2<label,vector>>>& getHistoryhistoryFt()   const
        {
            return historyFt_;
        }
        
        const vector& getHistoryCouplingF()   const
        {
            return historyCouplingF_;
        }
        
        const vector& getHistoryCouplingT()   const
        {
            return historyCouplingT_;
        }
        
        const List<DynamicLabelList>& getSurfaceCellList()        const 
        {
            return surfCells_;
        };
        
        const List<DynamicLabelList>& getInternalCellList()       const 
        {
            return intCells_;
        };
        
        const List<DynamicLabelList>& getWallContactFaceList()    const 
        {
            return wallContactFaces_;
        };
        
        const List<DynamicList<immersedBody::intVecRequest>>& getinterpolationVecReqs() const
        {
            return interpolationVecReqs_;
        };
        
        List<DynamicList<immersedBody::interpolationInfo>>& getInterpolationInfo()
        {
            return interpolationInfo_;
        };
        
        List<List<label>>& getBoundIndList()
        {
            return boundIndList_;
        };
        
        const point& getMinBoundPoint() const
        {
            return minBoundPoint_;
        };
        
        const point& getMaxBoundPoint() const
        {
            return maxBoundPoint_;
        };
        
        // pre-contact update of body field
        void preContactUpdateBodyField
        (
            volScalarField& body,
            volVectorField& f
        );
        
        // move the body post-contact
        void postContactUpdateBodyField(volScalarField& body, volScalarField& refineF);
        
        void recreateBodyField(volScalarField& body, volScalarField& refineF);
        
        //Update imposed vector field
        void updateVectorField
        (
            volVectorField& VS,
            word VName,
            volScalarField& body,
            vectorField surfNorm
        );
        
        // detect contact with other particles
        void detectWallContact(volScalarField& body);
        
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
        
        void preContactUpdateImmersedBody
        (
            volScalarField& body,
            volVectorField& f
        );
        
        void switchPrtContact (bool newBool){isInPrtContact_ =newBool;}
        void switchWallContact(bool newBool){isInWallContact_=newBool;}
        
        void postContactUpdateImmersedBody(volScalarField& body,volVectorField& f);
        void postContactUpdateImmersedBody(scalar deltaT);
        
        // compute the mean and maximal body courant number (based on
        // the body linear velocity)
        void computeBodyCoNumber();
        
        // print out body linear and angular momentum for control
        void printMomentum();
        
        // print out body statistics for control
        void printStats();
        
        // print out "eulerian" forces acting on the body for control
        void printForcesAndTorques();
        
        // function to compute local particle radii w.r.t. CoM_
        DynamicScalarList getLocPartRad(DynamicLabelList& cellsOfInt);
        
        // auxiliary switch to turn off/on torque update
        bool& retUpdateTorque()
        {
            return updateTorque_;
        }
        
        // return to history position
        void returnPosition();
        
        // initialize variables base to history
        void initializeVarHistory( bool setHistory);
        
        void solveWallContact
        (
            scalar kWN,
            scalar gammaWN,
            scalar kWt,
            scalar gammaWt,
            scalar muW,
            scalar adhWN,
            scalar deltaT
        );
        
        void switchActiveOff(volScalarField& body);
        const bool& getIsActive()
        {
            return isActive_;
        }
                
        bool checkContactMovement(scalar deltaT);
        
        void assignFullHistory();
        
        const label& getOwner() const
        {
            return owner_;
        }
        
        const bool& getStartSynced() const
        {
            return startSynced_;
        }
        
        Tuple2<vector,Tuple2<label,label>> findCellCustom
        (
            vector& prevPoint,
            label& startCell,
            label& startProc,
            vector& gradToBody,
            scalar& intDist
        );
        
        void initSyncWithFlow(const volVectorField& U);
        
        void pimpleUpdate
        (
            volScalarField& body,
            volVectorField& f
        );
        
        const bool& getSDBasedLambda() const
        {
            return sdBasedLambda_;
        }
        
        const scalar& getIntSpan() const
        {
            return intSpan_;
        }
        
        const bool& getUseInterpolation() const
        {
            return useInterpolation_;
        }
        
        void computeBodyCharPars();
        
        void checkIfInDomain(volScalarField& body);
        
        void setRecomputeProjection(const bool newBool)
        {
            recomputeProjection_ = newBool;
        }
        
        const label& getrecomputeM0() const
        {
            return recomputeM0_;
        }
        void recomputedM0()
        {
            recomputeM0_--;
        }
        
        void inCotactWithStatic(bool inContact)
        {
            inContactWithStatic_ = inContact;
        }
        
        void chceckBodyOp();
        
        void setRestartSim(vector vel, scalar angVel, vector axisRot, bool setStatic, label timesInContact);
        
        geomModel& getGeomModel()
        {
            return geomModel_();
        }
        
        void recordBodyInfo(IOdictionary& dict, word path)
        {
            dict.add("bodyId",bodyId_);
            dict.add("bodyName",bodyName_);
            dict.add("Vel",Vel_);
            dict.add("omega",omega_);
            dict.add("Axis",Axis_);
            dict.add("static",bodyOperation_ == 0);
            dict.add("timeStepsInContWStatic", timeStepsInContWStatic_);
            geomModel_->writeBodyInfo(bodyId_,dict, path);
        }        
        
        //update old movement variables (these are stored for pimple loop)
        void updateOldMovementVars()
        {
            VelOld_   = Vel_;
            AxisOld_  = Axis_;
            omegaOld_ = omega_;
        }
    };
}
#endif
