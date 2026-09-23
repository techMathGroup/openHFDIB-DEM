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
    implementation of the HFDIB method (Municchi and Radl, 2018) in OpenFOAM
    extended by connection with (level-set) DEM (Isoz and Sourek, 2022)


SourceFiles
    openHFDIBDEM.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef openHFDIBDEM_H
#define openHFDIBDEM_H

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
#include "triSurfaceMesh.H"
#include "immersedBody.H"

#include "addModel.H"
#include "addModelOnce.H"
#include "addModelOnceFromFile.H"
#include "addModelOnceScatter.H"
#include "addModelRepeatSamePosition.H"
#include "addModelRepeatRandomPosition.H"
#include "addModelDistribution.H"

#include "leastSquaresInt.H"
#include "lineInt.H"

#include "verletList.H"

//OS Time efficiency Testing
#include "clockTime.H"


namespace Foam
{

 class openHFDIBDEM
 {

    private:

    // Private data types ------------------------------------------

    typedef DynamicList<label>  DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point>  DynamicPointList;

    //Basic members
    const  fvMesh&          mesh_;
    IOdictionary            HFDIBDEMDict_;
    dictionary              HFDIBinterpDict_;
    IOdictionary            transportProperties_;

    const wordList          bodyNames_;

    // immersed bodies objects
    PtrList<immersedBody> immersedBodies_;

    // Interpolation class
    autoPtr<ibInterpolation> ibInterp_;

    // addModels (how to add the baseObjects)
    PtrList<addModel> addModels_;

    // precalculate cellPoints for better performace
    List<labelList> cellPoints_;

    HashTable<autoPtr<prtContactInfo>,Tuple2<label, label>,Hash<Tuple2<label, label>>>
        prtcInfoTable_;

    // DEM step
    scalar stepDEM_;

    // switch for the adaptive DEM time step (contact-potential test)
    bool adaptiveDEM_;

    // safety factors for the conservative sweep distance used when
    // sorting bodies into the sub-cycled / single-stepped sets:
    //  - translation: bounds |Vel|*dt (+ |a|*dt^2 for free bodies) with
    //    a margin for the force state drifting over the CFD step
    //  - rotation: bounds omega*rMax*dt with a larger margin - rotation
    //    moves contact points tangentially, which the tight-bbox tests
    //    of the live Verlet list only catch after the fact
    scalar sweepSafetyTrans_;
    scalar sweepSafetyRot_;

    // switch for the deltaTDEM estimator inside adaptiveStepDEM
    // (false restores the legacy stepDEM fraction of the CFD step)
    bool useDEMdtEstimator_;

    // pre-contact area estimate model: conservative uses areaCoeff*h^2,
    // velocityAware caps it by the hertz peak contact area
    // pi*R*deltaMax(v)
    bool dtEstimatorVelocityAware_;

    // c_A multiplier on h^2 in the pre-contact area estimate
    scalar dtAreaCoeff_;

    // multiplier on the tangential stiffness when sizing the step
    // (1 = full conservatism; the implemented tangential law carries
    // an empirical 200x tuning factor that makes it govern)
    scalar dtTangentialFactor_;

    // multiplier on the tangential stiffness of the rotational bound
    // only (escape hatch if the rotational estimate is too tight)
    scalar dtRotationalFactor_;

    // floor on the estimated deltaTDEM (0 = off). binding means the
    // estimated stability bound is being exceeded - warned in the log
    scalar dtMinDEMdt_;

    // cap on the DEM subcycle count (0 = unlimited). binding means the
    // estimated stability bound is being exceeded - warned in the log
    label dtMaxDEMSubCycles_;

    // report the used/hertz/rayleigh step line every N CFD steps
    label dtReportInterval_;

    // per-CFD-step DEM step estimate (stepDEM_ * min pair stability
    // step); GREAT when no pair constrains
    scalar deltaTDEM_;

    // hertz contact-duration and rayleigh surface-wave steps of the
    // same pair set (diagnostics only, never used for sizing)
    scalar dtDiagHertz_;
    scalar dtDiagRayleigh_;

    // how many times rocompute M0_ because of mesh refinement
    label recomputeM0_;

    // mesh geometric directions
    vector geometricD_;

    // Verlet list object
    verletList verletList_;

    //record simulation for start from lastTime possibility
    bool recordSimulation_;
    bool recordFirstTimeStep_ = false;
    word recordOutDir_;

    // counters, reset once per time step, tracking addRemoveBodies activity
    label nBodiesAddedLastStep_ = 0;
    label nBodiesRemovedLastStep_ = 0;

    public:

    // Constructors ----------------------------------------------------
    openHFDIBDEM(const Foam::fvMesh& mesh);

    // Destructors -----------------------------------------------------
    ~openHFDIBDEM();

    // Public methods --------------------------------------------------

    // scalar wallContactTime_ = 0.0;
    // scalar wallContactParallelTime_ = 0.0;
    // scalar wallContactReduceTime_ = 0.0;
    // scalar prtContactTime_ = 0.0;
    // scalar prtContactParallelTime_ = 0.0;
    // scalar prtContactReduceTime_ = 0.0;
    // scalar demItegrationTime_ = 0.0;

    void initialize(
        volScalarField& body,
        volVectorField& U,
        volScalarField& refineF,
        label recomputeM0, 
        word runTime
    );

    void createBodies
    (
        volScalarField& body,
        volScalarField& refineF
    );
    void recreateBodies
    (
        volScalarField& body,
        volScalarField& refineF
    );
    void createBodiesComputeDynamicsVars
    (
        volScalarField& body,
        volScalarField& refineF
    );

    void updateSurface
    (
        const scalar surfaceThreshold,
        const volScalarField& body,
        volScalarField& surface
    );

    void updateSurface
    (
        const scalar surfaceThreshold,
        const volScalarField& body,
        volScalarField& surface,
        surfaceScalarField& surfaceF
    );

    void updateGlobalFluidDensity(
        const volScalarField& body,
        volScalarField& rho
    );
    void updateGlobalFluidDensity(
        const volScalarField& body,
        volScalarField& rho,
        volScalarField& rhoS
    );

    void preUpdateBodies
    (
        volScalarField& body
    );
    void postUpdateBodies
    (
        const volScalarField& body,
        const volVectorField& f
    );
    void postUpdateBodies
    (
        const volScalarField& body,
        const volVectorField& f,
        const bool kinematicForce,
        const bool applyAddedMass
    );

    void postUpdateBodies
    (
        const volScalarField& body,
        const volVectorField& f,
        const volScalarField& rho,
        const bool kinematicForce,
        const bool applyAddedMass
    );

    void updateFSCoupling
    (
        volScalarField& body,
        volVectorField& f
    );
    void updateFSCoupling
    (
        volScalarField& body,
        volVectorField& f,
        const bool kinematicForce,
        const bool applyAddedMass
    );

    // void for DEM inner-loops
    void updateDEM(volScalarField& body,volScalarField& refineF);

    // adaptiveStepDEM: sort bodies into the sub-cycled set A (contact
    // potential with another body or a wall during the coming CFD time
    // step) and the single-stepped set B (no potential). Both sets
    // exclude static bodies and are rank-uniform (reduced over all
    // processors). hardFlags/sweepDist report the internal
    // classification for the DEM dt estimator.
    void computeAdaptiveSets
    (
        HashSet<label,Hash<label>>& subCycle,
        HashSet<label,Hash<label>>& singleStep,
        HashSet<label,Hash<label>>& hardFlags,
        HashTable<scalar,label,Hash<label>>& sweepDist
    );

    // compute the system DEM step from two pair classes: contact
    // (ongoing contact, cached contact state) and pre-contact
    // (sweep-inflated pairs, estimated geometry). hertz/rayleigh
    // steps are diagnostics only. the min and its governing pair are
    // reduced across ranks
    void computeDEMdtEstimate
    (
        const HashSet<label,Hash<label>>& subCycle,
        const HashSet<label,Hash<label>>& hardFlags,
        const HashTable<scalar,label,Hash<label>>& sweepDist
    );

    // upper bound on the contact-point speed magnitude of a body
    // (|Vel| + |omega|*rMax), used for the pre-contact approach speed
    scalar pairApproachSpeed(immersedBody& ib) const;

    // upper bound on the lever arm of any contact point of the body
    // (bbox corner distance from the CoM); static bodies
    // (bodyOperation_ == 0) return 0 - they are not integrated
    // rotationally
    scalar bodyRMax(immersedBody& ib) const;

    // smallest principal moment of inertia (lambda_min of I); GREAT
    // for static bodies so their rotational bound never governs
    scalar bodyIeff(immersedBody& ib) const;

    // drop all contact bookkeeping of  body bodyId: the
    // prtcInfoTable_ entries of every pair containing it (the
    // Verlet-driven erase paths no longer see the body) and its
    // Verlet-list entry.
    // Note (MI): used in clean-up of switched-off bodies
    void removeBodyContacts(label bodyId);

    prtContactInfo& getPrtcInfo(Tuple2<label,label> cPair);

    // interpolate vector field V using Vs for the surface values
    void interpolateIB
    (
        volVectorField& V,
        volVectorField& Vs,
        volScalarField& body
    );

    // save the moved STL files
    void writeBodiesInfo();

    // add or remove bodies from the simulation
    void addRemoveBodies(volScalarField& body,volVectorField& U,volScalarField& refineF);

    void restartSimulation
    (
        volScalarField& body,
        volScalarField& refineF,
        word runTime
    );

    void preCalculateCellPoints();

    bool getRecordFirstTime()
    {
        return recordFirstTimeStep_;
    }

    void setRecordFirstTime( bool setUp)
    {
        recordFirstTimeStep_ = setUp;
    }

    // number of bodies added by addRemoveBodies during the last time step
    label nBodiesAddedLastStep() const
    {
        return nBodiesAddedLastStep_;
    }

    // number of bodies switched off (left the domain) during the last time step
    label nBodiesRemovedLastStep() const
    {
        return nBodiesRemovedLastStep_;
    }

    void writeFirtsTimeBodiesInfo();

    void setSolverInfo();
    
    void updateBodiesRhoF(const scalar rho);
    void updateBodiesRhoF(const volScalarField& rho, const volScalarField& body);

    scalar computeBodiesLinCourantNo();
 };
}

#endif
