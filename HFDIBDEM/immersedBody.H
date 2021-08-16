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

#include "contactInfo.H"

namespace Foam
{

    class immersedBody
    {

        private:

        // Private data types ------------------------------------------

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        // name of the body
        word    bodyName_;

        // is the particle active in the simulation
        bool    isActive_;

        // body set-up
        dictionary&             immersedDict_;

        // mesh reference
        const  dynamicFvMesh&   mesh_;

        // fluid properties
        dictionary&             transportProperties_;

        // body geometry
        autoPtr<geomModel> geomModel_;

        // contact models
        autoPtr<contactInfo> contactInfo_;

        // precalculate cellPoints for better performace
        List<pointField>& cellPoints_;

        // rotation axis
        vector Axis_;
        vector AxisOld_;

        // angular velocity (scalar)
        scalar omega_;
        scalar omegaOld_;

        // linear velocity
        vector Vel_;
        vector VelOld_;

        // linear acceleration
        vector a_;

        // angular acceleration
        vector alpha_;

        // total Rodrigues angle of the body
        vector totalAngle_;

        // forces acting on the particle
        vector  F_;

        // total torque acting on the particle
        vector  T_;

        // maximal body Courant number based on its linear velocity
        scalar CoNum_;

        // mean body Courant number based on its linear velocity
        scalar meanCoNum_;

        // fluid density
        dimensionedScalar rhoF_;

        // integer to identify current body
        label   bodyId_;

        // auxiliary switch - do I want to update the torques?
        bool    updateTorque_;

        // operation to perform
        int     bodyOperation_;

        // fraction of dC_to move by in one DEM loop
        scalar  maxDistInDEMloop_;

        // index holder for particle contact list
        List<List<label>> boundIndList_;

        // bounding points of the surface
        point minBoundPoint_;
        point maxBoundPoint_;

        // body surface cells
        List<DynamicLabelList>  surfCells_;

        // body internal cells
        List<DynamicLabelList>  intCells_;

        // processor num that owns the IB
        label owner_;

        // relaxation factor for velocity
        scalar velRelaxFac_;

        // update immersed body info
        void  updateImmersedBody
        (
            volScalarField& body,
            volVectorField& f
        );

        // create interpolation points
        void  calculateInterpolationPoints(volScalarField& body);

        // update fluid-solid coupling
        void updateCoupling
        (
            volScalarField& body,
            volVectorField& f
        );

        // history variables for contact correction
        // -- forces acting on the particle
        vector  historyCouplingF_;

        // -- total torque acting on the particle
        vector  historyCouplingT_;

        // -- rotation axis
        vector historyAxis_;

        // -- angular velocity
        scalar historyOmega_;

        // -- linear velocity
        vector historyVel_;

        // -- linear acceleration
        vector historya_;

        // -- angular acceleration
        vector historyAlpha_;

        // -- total Rodrigues angle of the body
        vector historyTotalAngle_;

        // -- tangential force from last timestep
        DynamicList<Tuple2<label,Tuple2<label,vector>>> historyFt_;

        // auxiliary octree variables
        Field<label> octreeField_;
        label cellToStartInCreateIB_;

        // start the body synchronized with the flow (approximate)
        bool    startSynced_;
        
        // auxiliary variables
        tensor  totRotMatrix_;
        bool    sdBasedLambda_;
        scalar  intSpan_;
        scalar  charCellSize_;                                          //characteristic cell size over surface
        label   refineBuffers_;                                         //number of refine buffers around body
        bool    useInterpolation_;                                      //number of refine buffers around body
        label   recomputeM0_;                                           //how many times rocompute M0_ because of refinement
        vector  geometricD_;                                            //geometric directions
        bool    created_;
        label   timesToSetStatic_;

        // immersed body initialization
        void initializeIB();

        protected:


        public:

        // Constructors ------------------------------------------------
        // -- general constructor for plug-in geomModel
        immersedBody
        (
            word bodyName,
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
        // -- class to keep interpolation data
        class interpolationInfo{
        public:
            label surfCell_;
            DynamicPointList intPoints_;
            label order_;
            DynamicLabelList intCells_;
            DynamicLabelList procWithIntCells_;
            DynamicVectorList intVec_;
        };

        // -- class to ask for interpolation data on other processor
        class intVecRequest
        {
            public:
                label requestLabel_;
                label vecLabel_;
                point intPoint_;
                label intCell_;
        };

        // Public variables --------------------------------------------
        // -- interpolation info list
        List<DynamicList<immersedBody::interpolationInfo>>  interpolationInfo_;

        // -- interpolation vec request list
        List<DynamicList<immersedBody::intVecRequest>>      interpolationVecReqs_;

        // Public methods ----------------------------------------------

        // move immersed body
        void moveImmersedBody
        (
            scalar deltaT = -1.0
        );

        // reset body field
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

        // create immersed body
        void createImmersedBody
        (
            volScalarField& body,
            volScalarField& refineF,
            bool synchCreation = true
        );

        void syncCreateImmersedBody
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

        // getters

        const bool& checkWallContact() const
        {
            return contactInfo_->checkWallContact();
        }

        const label& getBodyId() const
        {
            return bodyId_;
        }

        const scalar& getM0()    const
        {
            return geomModel_->getM0();
        }

        const scalar& getM()    const
        {
            return geomModel_->getM();
        }

        const dimensionedScalar& getRhoS()    const
        {
            return geomModel_->getRhoS();
        }

        const vector& getVel()  const
        {
            return Vel_;
        }

        const vector& getCoM()  const
        {
            return geomModel_->getCoM();
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
            return geomModel_->getI();
        }

        const scalar& getKN()   const
        {
            return contactInfo_->getkN();
        }

        const scalar& getGammaN()   const
        {
            return contactInfo_->getgammaN();
        }

        const scalar& getKt()   const
        {
            return contactInfo_->getkt();
        }

        const scalar& getGammat()   const
        {
            return contactInfo_->getgammat();
        }

        const scalar& getmu()   const
        {
            return contactInfo_->getmu();
        }

        const scalar& getadhN()   const
        {
            return contactInfo_->getAdhN();
        }

        const scalar& getadhEqui()   const
        {
            return contactInfo_->getAdhEqui();
        }

        const scalar& getCoNum()   const
        {
            return CoNum_;
        }

        const scalar& getMeanCoNum()   const
        {
            return meanCoNum_;
        }

        scalar getDC()
        {
            return geomModel_->getDC();
        }

        const int& getbodyOperation()   const
        {
            return bodyOperation_;
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

        const label& getTimeStepsInContWStatic() const
        {
            return contactInfo_->getTimeStepsInContWStatic();
        };

        // move the body post-contact
        void postContactUpdateBodyField(volScalarField& body, volScalarField& refineF);

        void recreateBodyField(volScalarField& body, volScalarField& refineF);

        // update imposed vector field
        void updateVectorField
        (
            volVectorField& VS,
            word VName,
            volScalarField& body,
            vectorField surfNorm
        );

        // should we detect contact with wall
        bool shouldDetectWallContact();

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
        void setWallContact(bool inContact){contactInfo_->setWallContact(inContact);;}

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

        void computeBodyCharPars(){geomModel_->computeBodyCharPars(surfCells_);}

        void checkIfInDomain(volScalarField& body);

        const label& getrecomputeM0() const
        {
            return recomputeM0_;
        }
        void recomputedM0()
        {
            recomputeM0_--;
        }

        void inContactWithStatic(bool inContact)
        {
            contactInfo_->inContactWithStatic(inContact);
        }

        void chceckBodyOp();

        void setRestartSim(vector vel, scalar angVel, vector axisRot, bool setStatic, label timesInContact);

        geomModel& getGeomModel()
        {
            return geomModel_();
        }

        contactInfo& getContactInfo()
        {
            return contactInfo_();
        }

        contactVars* getContactVars()
        {
            contactVars* cVars = new contactVars(
                bodyId_,
                Vel_,
                omega_,
                Axis_,
                getM0(),
                getM(),
                geomModel_->getRhoS()
            );
            return cVars;
        }

        void recordBodyInfo(IOdictionary& dict, word path)
        {
            dict.add("bodyId",bodyId_);
            dict.add("bodyName",bodyName_);
            dict.add("Vel",Vel_);
            dict.add("omega",omega_);
            dict.add("Axis",Axis_);
            dict.add("static",bodyOperation_ == 0);
            dict.add("timeStepsInContWStatic", getTimeStepsInContWStatic());
            geomModel_->writeBodyInfo(bodyId_,dict, path);
        }

        // update old movement variables for fluid-solid coupling
        void updateOldMovementVars()
        {
            VelOld_   = Vel_;
            AxisOld_  = Axis_;
            omegaOld_ = omega_;
        }
    };
}
#endif
