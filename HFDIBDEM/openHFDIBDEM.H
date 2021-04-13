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
    Martin Isoz (2019-*), Martin Šourek (2019-*)
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
#include "dynamicFvMesh.H"
#include "triSurfaceMesh.H"
#include "immersedBody.H"

#include "addModel.H"
#include "addModelOnce.H"
#include "addModelOnceScatter.H"
#include "addModelRepeatSamePosition.H"
#include "addModelRepeatRandomPosition.H"
#include "addModelDistribution.H"

namespace Foam
{
    
 class openHFDIBDEM
 {

    private:
  
    // Private data types

    typedef DynamicList<label>  DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point>  DynamicPointList;

    //Basic members
    const  dynamicFvMesh&   mesh_;
    IOdictionary            HFDIBDEMDict_;
    dictionary              HFDIBinterpDict_;
    IOdictionary            transportProperties_;
    
    const wordList          stlNames_;

    //Immersed bodies objects
    PtrList<immersedBody> immersedBodies_;
    
    //addModels (how to add the baseObjects)
    PtrList<addModel> addModels_;
    
    //Contact variables (global)
    //spring stiffness for normal contact force (walls)
    scalar  kWN_;
    
    //damping coefficient for normal contact force (walls)
    scalar gammaWN_;
    
    //spring stiffness for normal contact force (walls)
    scalar  kWt_;
    
    //damping coefficient for normal contact force (walls)
    scalar gammaWt_;
    
    //friction coefficient
    scalar muW_;
    
    //adhesion coefficient
    scalar adhWN_;
    
    DynamicLabelList ibContactList_;
    
    //Minimal number of DEM inside loops
    label minDEMloops_;
    
    //Minimal timestep in DEM inside loops
    scalar minDEMtimeStep_;
    
    //how many times rocompute M0_ because of refinement
    label recomputeM0_;
    
    //solve body-wall contact
    //~ void solveWallContact(const DynamicLabelList& wallContactLst);
    // Note (MI): moved to immersedBody class by MS (as it should be there)
    
    //Neighbour list variables - First list is related to dimensions so the size is always 3
    // List for all bounding points positions
    DynamicList<DynamicScalarList> boundValueNeighbourList_;
    // List for labels holding information about bounding points in boundValueNeighbourList_
    // This is only list in which bodyID is increased by 1
    // Positive label indicates maximal bounding point of body Label-1
    // Negative label indicates minimal bounding point of body -label-1
    DynamicList<DynamicLabelList> boundLabelNeighbourList_;
    // List for pairs of intersected Bodies in specific coordinate
    DynamicList<DynamicList<Tuple2<label, label>>> contactInCoordNeighbourList_;
    // List of possible contact du to intersection of bounding boxes.
    // Pair is added to this list when it is present in all 3 contactInCoordNeighbourList_.
    DynamicList<Tuple2<label, label>> possibleContactNeighbourList_;
    // List of pairs that are in contact due to intersection of surfaces
    DynamicList<Tuple2<label, label>> prtContactIBList_;
    // Center of the contact area
    DynamicVectorList prtContactCenter_;
    //  Penetration volume for prt-prt contact
    DynamicScalarList prtContactVolume_;
    //  Normal direction for prt-prt contact
    DynamicVectorList prtContactNormal_;
    //  Contact area for prt-prt contact
    DynamicScalarList prtContactArea_;
    // number of DEM loops
    DynamicScalarList numberOfDEMloops_;
    // Method for iteration over list and find if the pair is present in given list
    label findIndexOfPairInNeighbourList(DynamicList<Tuple2<label, label>>& listToSearch, Tuple2<label, label> pair);
    // Update positions of bounding points
    // sort the list
    // Solve potential intersections
    void updateNeighbourLists();
    // detect prt-prt contact
    void detectPrtContact();
    
    protected:
    
    class prtPrtContactInfo{
    public:
        Tuple2<label, label> prtsInContact_;
        vector contactCenter_;
        scalar contactVolume_;
        vector contactNormal_;
        scalar contactArea_;
        bool inContact_;
    };
    
    public:
    
    openHFDIBDEM(const Foam::dynamicFvMesh& mesh);
    openHFDIBDEM(const Foam::dynamicFvMesh& mesh, const wordList stlNames);
    ~openHFDIBDEM();
    
    //~ void initialize(volScalarField& body);
    void initialize(volScalarField& body,volVectorField& U,volScalarField& refineF,label recomputeM0);
    
    void preUpdateBodies
    (
        volScalarField& body,
        volVectorField& f
    );
    scalar postUpdateBodies(volScalarField& body,volVectorField& f);
    // Note (MI): this function returns the maximal CoNum encountered
    //            between particles (wanted to maintain a stable integration)
    
    void updateFSCoupling
    (
        volScalarField& body,
        volVectorField& f
    );
    
    void moveBodies(volScalarField& body,volScalarField& refineF);
    
    void recreateBodies(volScalarField& body,volScalarField& refineF);
    // Void for DEM inner-loops
    void correctContact(volScalarField& body,volScalarField& refineF);
    // Detect wall contact at the end of the time step and return to original positions
    void getContactListAndReturnPositions(volScalarField& body);
    // Sort list of bounding points
    void sortBoundingListPrtContact();
    // Solve possible intersection of bounding boxes
    void swapInBoundingListPrtContact(label coord, label j);
    
    //Interpolate vector V using Vs for the surface values
    void interpolateIB
    (
        volVectorField& V,
        volVectorField& Vs,
        volScalarField& body
    );
    
    //Interpolate scalar S using Ss for the surface values
    void interpolateIB
    (
        volScalarField& S,
        volScalarField& Ss,
        volScalarField& body
    )
    {
        FatalError << "not implemented" << exit(FatalError);
    }
    
    // save the moved STL files
    void writeBodySurfMeshes();
    
    // add or remove bodies from the simulation
    //~ void addRemoveBodies(volScalarField& body);
    void addRemoveBodies(volScalarField& body,volVectorField& U,volScalarField& refineF);
    
    // function to find normal vector
    vector getContactNormal
    (
        DynamicLabelList commonCells,
        vector normalVector,vector contactCenter
    );
    
    // function to find contact area
    scalar getContactArea(DynamicLabelList commonCells);
    
    // function to return variables needed for prt-prt resolution
    openHFDIBDEM::prtPrtContactInfo getPrtContactInfo(Tuple2<label, label> pairToTest);
    
    // Method for iteration over list and find if pairs that contains at least some ib
    void findIndexesOfPairWithSomeIb
    (
        DynamicList<Tuple2<label, label>>& listToSearch,
        Tuple2<label, label> pair,
        DynamicList<Tuple2<label,
        label>>& listToAppend
    );    
        
    //solve body-body contact
    void solvePrtContact(openHFDIBDEM::prtPrtContactInfo contactInfo, scalar deltaT);
    
    //Resolve possible 3D contact info
    Tuple2<scalar,vector> get3DcontactInfo
    (
        DynamicLabelList commonCells,
        vector normalVector,
        vector contactCenter,
        label owner
    );

 };

}

#endif
