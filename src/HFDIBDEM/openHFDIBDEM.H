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
    implementation of the HFDIB method (Municchi and Radl, 2016) in OpenFOAM
    extended by connection with (level-set) DEM (Isoz and Sourek, 2020)


SourceFiles
    openHFDIBDEM.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*)
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
//// #include "fvMesh.H"
#include "triSurfaceMesh.H"
#include "immersedBody.H"

#include "addModel.H"
#include "addModelOnce.H"
#include "addModelOnceFromFile.H"
#include "addModelOnceScatter.H"
#include "addModelRepeatSamePosition.H"
#include "addModelRepeatRandomPosition.H"
#include "addModelDistribution.H"

// #include "contactModel.H"
#include "leastSquaresInt.H"
#include "lineInt.H"
// #include "IOstream.H"
// #include "virtualMesh.H"

#include "verletList.H"
//OS Time effitiency Testing
#include "clockTime.H"
//OS Time effitiency Testing
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

    void initialize(volScalarField& body,volVectorField& U,volScalarField& refineF,label recomputeM0, word runTime);

    void createBodies
    (
        volScalarField& body,
        volScalarField& refineF
    );

    void preUpdateBodies
    (
        volScalarField& body,
        volVectorField& f
    );
    void postUpdateBodies(volScalarField& body,volVectorField& f);

    void updateFSCoupling
    (
        volScalarField& body,
        volVectorField& f
    );

    void recreateBodies(volScalarField& body,volScalarField& refineF);
    // void for DEM inner-loops
    void updateDEM(volScalarField& body,volScalarField& refineF);

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

    void writeFirtsTimeBodiesInfo();
    
    void setSolverInfo();

 };
}

#endif
