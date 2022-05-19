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
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*)
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
// #include "fvMesh.H"
#include "triSurfaceMesh.H"

#include "geomModel.H"
#include "stlBased.H"
#include "convexBody.H"
#include "nonConvexBody.H"

#include "ibContactClass.H"
#include "wallContactInfo.H"
#include "interpolationInfo.H"
#include "ibInterpolation.H"

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
        const  fvMesh&   mesh_;

        // fluid properties
        dictionary&             transportProperties_;

        // body geometry
        autoPtr<geomModel> geomModel_;

        // Interpolation Info
        autoPtr<interpolationInfo> intpInfo_;

        // contact models
        autoPtr<ibContactClass> ibContactClass_;

        autoPtr<ibContactVars> ibCntVars_;

        autoPtr<wallContactInfo> wallCntInfo_;

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

        // index holder for particle contact list
        List<List<label>> boundIndList_;

        // bounding points of the surface
        point minBoundPoint_;
        point maxBoundPoint_;

        // processor num that owns the IB
        label owner_;

        // update immersed body info
        void  updateImmersedBody
        (
            volScalarField& body,
            volVectorField& f
        );

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
        label   recomputeM0_;                                           //how many times rocompute M0_ because of refinement
        label   timesToSetStatic_;
        vector  staticContactPost_;

        protected:


        public:

        // Constructors ------------------------------------------------
        // -- general constructor for plug-in geomModel
        immersedBody
        (
            word bodyName,
            const Foam::fvMesh& mesh,
            const volScalarField& body,
            dictionary& HFDIBDEMDict,
            dictionary& transportProperties,
            label bodyId,
            label recomputeM0,
            geomModel* bodyGeomModel,
            autoPtr<ibInterpolation>& ibIntp,
            List<pointField>& cellPoints,
            HashTable<materialInfo,string,Hash<string>>& materialInfos,
            HashTable<materialInfo,string,Hash<string>>& wInfos,
            HashTable<scalar,string,Hash<string>>& matInterAdh
        );

        // Destructors -------------------------------------------------
        ~immersedBody();

        // Public subclasses -------------------------------------------
        // -- class to keep interpolation data
        class interpolationInfoOld{
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
        List<DynamicList<immersedBody::interpolationInfoOld>>  interpolationInfoOld_;

        // -- interpolation vec request list
        List<DynamicList<immersedBody::intVecRequest>>      interpolationVecReqs_;

        // Public methods ----------------------------------------------

        // move immersed body
        void moveImmersedBody
        (
            scalar deltaT = -1.0
        );

        // reset body field
        void resetBody(volScalarField& body);

        //update body movement based on F_ and T_
        void updateMovement
        (
            scalar deltaT
        );
        void updateMovement
        (
            vector Vel,
            vector Axis,
            scalar omega
        );

        void updateMovementComp
        (
            scalar deltaT,
            vector Vel,
            vector Axis,
            scalar omega
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

        interpolationInfo& getIntpInfo()
        {
            return intpInfo_();
        }

        void clearIntpInfo()
        {
            if(intpInfo_.valid())
            {
                intpInfo_().clearIntpInfo();
            }
        }

        const bool& checkWallContact() const
        {
            return ibContactClass_->checkWallContact();
        }

        const label& getBodyId() const
        {
            return bodyId_;
        }

        const vector& getVel()  const
        {
            return Vel_;
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
            return geomModel_->getSurfaceCellList();
        };

        const List<DynamicLabelList>& getInternalCellList()       const
        {
            return geomModel_->getInternalCellList();
        };

        const List<DynamicList<immersedBody::intVecRequest>>& getinterpolationVecReqs() const
        {
            return interpolationVecReqs_;
        };

        List<DynamicList<immersedBody::interpolationInfoOld>>& getInterpolationInfo()
        {
            return interpolationInfoOld_;
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

        void postContactUpdateBodyField(volScalarField& body, volScalarField& refineF);

        void recreateBodyField(volScalarField& body, volScalarField& refineF);

        // update imposed vector field
        void updateVectorField
        (
            volVectorField& VS,
            word VName,
            volScalarField& body
        );

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

        void setWallContact(bool inContact){ibContactClass_->setWallContact(inContact);;}

        void postPimpleUpdateImmersedBody(volScalarField& body,volVectorField& f);

        // compute the mean and maximal body courant number (based on
        // the body linear velocity)
        void computeBodyCoNumber();

        // print out body linear and angular momentum for control
        void printMomentum();

        // print out body statistics for control
        void printStats();

        // print out "eulerian" forces acting on the body for control
        void printForcesAndTorques();

        // auxiliary switch to turn off/on torque update
        bool& retUpdateTorque()
        {
            return updateTorque_;
        }

        // initialize variables base to history
        void resetForces()
        {
            F_*=0.0;
            T_*=0.0;
        }

        void switchActiveOff(volScalarField& body);
        const bool& getIsActive()
        {
            return isActive_;
        }

        const label& getOwner() const
        {
            return owner_;
        }

        const bool& getStartSynced() const
        {
            return startSynced_;
        }

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

        void computeBodyCharPars(){geomModel_->computeBodyCharPars();}

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
            ibContactClass_->inContactWithStatic(inContact);
        }

        void chceckBodyOp();

        void setRestartSim(vector vel, scalar angVel, vector axisRot, bool setStatic, label timesInContact);

        geomModel& getGeomModel()
        {
            return geomModel_();
        }

        autoPtr<geomModel>& getGeomModelPtr()
        {
            return geomModel_;
        }

        ibContactClass& getibContactClass()
        {
            return ibContactClass_();
        }

        ibContactVars& getContactVars()
        {
            return ibCntVars_();
        }

        wallContactInfo& getWallCntInfo()
        {
            return wallCntInfo_();
        }

        void recordBodyInfo(IOdictionary& dict, word path)
        {
            dict.add("bodyId",bodyId_);
            dict.add("bodyName",bodyName_);
            dict.add("Vel",Vel_);
            dict.add("omega",omega_);
            dict.add("Axis",Axis_);
            dict.add("static",bodyOperation_ == 0);
            dict.add("timeStepsInContWStatic", ibContactClass_->getTimeStepsInContWStatic());
            geomModel_->writeBodyInfo(bodyId_,dict, path);
        }

        // update old movement variables for fluid-solid coupling
        void updateOldMovementVars()
        {
            VelOld_   = Vel_;
            AxisOld_  = Axis_;
            omegaOld_ = omega_;
        }

        void writeIbInfo()
        {
            InfoH << iB_Info;
            InfoH << "-- body " << bodyId_ << " center of mass       : " << geomModel_->getCoM() << endl;
            InfoH << "-- body " << bodyId_ << " linear velocity      : " << Vel_ << endl;
            InfoH << "-- body " << bodyId_ << " angluar velocity     : " << omega_ << endl;
            InfoH << "-- body " << bodyId_ << " axis of rotation     : " << Axis_ << endl;
        }

        vectorField getUatIbPoints();
    };
}
#endif
