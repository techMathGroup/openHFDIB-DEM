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
    class for immersed bodies representation.
SourceFiles
    immersedBodies.C
Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Šourek (2019-*)
\*---------------------------------------------------------------------------*/
#include "immersedBody.H"
#include "polyMesh.H"
#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"

#include "interpolationCellPoint.H"
#include "interpolationCell.H"
#include "meshSearch.H"
#include "List.H"
#include "ListOps.H"

#include "OFstream.H"

#include "cellSet.H"
#include "addToRunTimeSelectionTable.H"
#include "PstreamReduceOps.H"

#include "fvcSmooth.H"
#include "fvMeshSubset.H"

#define ORDER 2

using namespace Foam;

//---------------------------------------------------------------------------//
immersedBody::immersedBody
(
    word fileName,
    const Foam::dynamicFvMesh& mesh,
    dictionary& HFDIBDEMDict,
    dictionary& transportProperties,
    label bodyId,
    label recomputeM0,
    vector geometricD
)
:
stlName_(fileName),
isInPrtContact_(false),
isInWallContact_(false),
isActive_(true),
recomputeProjection_(true),
immersedDict_(HFDIBDEMDict.subDict(stlName_)),
mesh_(mesh),
transportProperties_(transportProperties),
M_(0.0),
M0_(0.0),
CoM_(vector::zero),
Axis_(vector::one),
AxisOld_(vector::one),
omega_(0.0),
omegaOld_(0.0),
Vel_(vector::zero),
VelOld_(vector::zero),
a_(vector::zero),
alpha_(vector::zero),
totalAngle_(vector::zero),
I_(symmTensor::zero),
F_(vector::zero),
T_(vector::zero),
kN_(readScalar(immersedDict_.lookup("kN"))),
kt_(readScalar(immersedDict_.lookup("kt"))),
gammaN_(readScalar(immersedDict_.lookup("gammaN"))),
gammat_(readScalar(immersedDict_.lookup("gammat"))),
mu_(readScalar(immersedDict_.lookup("mu"))),
adhN_(readScalar(immersedDict_.lookup("adhN"))),
CoNum_(0.0),
rhoF_(transportProperties_.lookup("rho")),
rhoS_(immersedDict_.lookup("rho")),
bodyConvex_(false),
dC_(0.0),
thrSurf_(readScalar(HFDIBDEMDict.lookup("surfaceThreshold"))),
bodyId_(bodyId),
updateTorque_(false),
bodyOperation_(0),
maxDistInDEMloop_(readScalar(HFDIBDEMDict.lookup("maxDistInDEMloop"))),
boundIndList_(3),
bodySurfMesh_
(
    IOobject
    (
        stlName_ +".stl",
        "constant",
        "triSurface",
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
writeBodySurfMesh_(readBool(immersedDict_.lookup("writeBodySTL"))),
owner_(false),
historyCouplingF_(vector::zero),
historyCouplingT_(vector::zero),
historyAxis_(vector::one/mag(vector::one)),
historyOmega_(0.0),
historyVel_(vector::zero),
historya_(vector::zero),
historyAlpha_(vector::zero),
historyTotalAngle_(vector::zero),
contactTimeRes_(1),
octreeField_(mesh_.nCells(), 0),
cellToStartInCreateIB_(0),
totRotMatrix_(tensor::I),
sdBasedLambda_(false),
intSpan_(2.0),
dS_(0.0),
charCellSize_(1e3),
refineBuffers_(0),
useInterpolation_(true),
recomputeM0_(recomputeM0),
geometricD_(geometricD)
{
    initializeIB();    
}

immersedBody::immersedBody
(
    triSurface triSurf,
    word    surfName,
    const Foam::dynamicFvMesh& mesh,
    dictionary& HFDIBDEMDict,
    dictionary& transportProperties,
    label bodyId,
    label recomputeM0,
    vector geometricD
)
:
stlName_(surfName),
isInPrtContact_(false),
isInWallContact_(false),
isActive_(true),
recomputeProjection_(true),
immersedDict_(HFDIBDEMDict.subDict(stlName_)),
mesh_(mesh),
transportProperties_(transportProperties),
M_(0.0),
M0_(0.0),
CoM_(vector::zero),
Axis_(vector::one),
AxisOld_(vector::one),
omega_(0.0),
omegaOld_(0.0),
Vel_(vector::zero),
VelOld_(vector::zero),
a_(vector::zero),
alpha_(vector::zero),
totalAngle_(vector::zero),
I_(symmTensor::zero),
F_(vector::zero),
T_(vector::zero),
kN_(readScalar(immersedDict_.lookup("kN"))),
kt_(readScalar(immersedDict_.lookup("kt"))),
gammaN_(readScalar(immersedDict_.lookup("gammaN"))),
gammat_(readScalar(immersedDict_.lookup("gammat"))),
mu_(readScalar(immersedDict_.lookup("mu"))),
adhN_(readScalar(immersedDict_.lookup("adhN"))),
adhEqui_(readScalar(immersedDict_.lookup("adhEqui"))),
CoNum_(0.0),
rhoF_(transportProperties_.lookup("rho")),
rhoS_(immersedDict_.lookup("rho")),
bodyConvex_(false),
dC_(0.0),
thrSurf_(readScalar(HFDIBDEMDict.lookup("surfaceThreshold"))),
bodyId_(bodyId),
updateTorque_(false),
bodyOperation_(0),
maxDistInDEMloop_(readScalar(HFDIBDEMDict.lookup("maxDistInDEMloop"))),
boundIndList_(3),
bodySurfMesh_
(
    IOobject
    (
        stlName_ +".stl",
        "constant",
        "triSurface",
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    triSurf
),
writeBodySurfMesh_(readBool(immersedDict_.lookup("writeBodySTL"))),
owner_(false),
velRelaxFac_(1.0),
historyCouplingF_(vector::zero),
historyCouplingT_(vector::zero),
historyAxis_(vector::one/mag(vector::one)),
historyOmega_(0.0),
historyVel_(vector::zero),
historya_(vector::zero),
historyAlpha_(vector::zero),
historyTotalAngle_(vector::zero),
contactTimeRes_(1),
octreeField_(mesh_.nCells(), 0),
cellToStartInCreateIB_(0),
startSynced_(false),
totRotMatrix_(tensor::I),
sdBasedLambda_(false),
intSpan_(2.0),
dS_(0.0),
charCellSize_(1e3),
refineBuffers_(0),
useInterpolation_(true),
recomputeM0_(recomputeM0),
geometricD_(geometricD)
{
    initializeIB();    
}
//---------------------------------------------------------------------------//
immersedBody::~immersedBody()
{
    bodySurfMesh_.clearOut();
}
//---------------------------------------------------------------------------//
void immersedBody::initializeIB()
{
    Info<< "Data on immersed boundary object triSurface" << endl;
    bodySurfMesh_.writeStats(Info);
    
    if (mesh_.nGeometricD() < 3)
    {
        const vector validDirs = (mesh_.geometricD() + Vector<label>::one)/2;
        Axis_ -= validDirs;
    }
    Axis_ /= mag(Axis_);
    
    // read the body operation name from immersedDict_
    if (immersedDict_.found("staticBody"))
    {
        bodyOperation_ = 0;
        Info << stlName_ << " is a static body." << endl;
    }
    else if (immersedDict_.found("prescribedTransBody"))
    {
        bodyOperation_ = 1;
        
        Vel_   = immersedDict_.subDict("prescribedTransBody").lookup("velocity");
        
        Info << stlName_ << " is a freely rotating body with prescribed linear velocity." << endl;
    }
    else if (immersedDict_.found("prescribedRotBody"))
    {
        bodyOperation_ = 2;
        
        Axis_ = immersedDict_.subDict("prescribedRotBody").lookup("axis");       
        omega_  = readScalar(
                    immersedDict_.subDict("prescribedRotBody").lookup("omega")
                );
    
        Info << stlName_ << " is a freely moving body with prescribed rotation." << endl;
        
        // Note: at the moment, the body always rotates around its center of mass
    }
    else if (immersedDict_.found("prescribedTransRotBody"))
    {
        bodyOperation_ = 3;
        
        Vel_   = immersedDict_.subDict("prescribedTransRotBody").lookup("velocity");
        Axis_ = immersedDict_.subDict("prescribedTransRotBody").lookup("axis");       
        omega_  = readScalar(
                    immersedDict_.subDict("prescribedTransRotBody").lookup("omega")
                );
    
        Info << stlName_ << " has prescribed both movement and rotation." << endl;
        
        // Note: at the moment, the body always rotates around its center of mass
    }
    else if (immersedDict_.found("prescribedTransFixedAxisRotBody"))
    {
        bodyOperation_ = 4;
        
        Vel_   = immersedDict_.subDict("prescribedTransFixedAxisRotBody").lookup("velocity");
        Axis_ = immersedDict_.subDict("prescribedTransFixedAxisRotBody").lookup("axis");
        
        Info << stlName_ << " has prescribed movement and axis of rotation." << endl;
        
        // Note: at the moment, the body always rotates around its center of mass
    }
    else if (immersedDict_.found("fullyCoupledBody"))
    {
        bodyOperation_ = 5;
        
        Info << stlName_ << " is fully coupled with fluid phase." << endl;
    }
    else
    {
        Info << "No body operation was found for " << stlName_ << endl
          << "Assuming static body.";
    }
    
    // check if the immersedDict_ contains switch for updateTorque_
    if (immersedDict_.found("updateTorque"))
    {
        updateTorque_ = readBool(immersedDict_.lookup("updateTorque"));
        Info << "Found updateTorque: " << updateTorque_ << endl;
    }
    else
    {
        Info << "Did not find updateTorque, using updateTorque: " << updateTorque_ << endl;
    }
    // check if the immersedDict_ contains switch for bodyConvex_
    if (immersedDict_.found("bodyConvex"))
    {
        bodyConvex_ = readBool(immersedDict_.lookup("bodyConvex"));
        Info << "Found bodyConvex, the body is convex?: " << bodyConvex_ << endl;
    }
    else
    {
        Info << "Did not find bodyConvex, using bodyConvex: " << bodyConvex_ << endl;
    }
    
    // do I want the IB to start as in sync with the flow?
    if (immersedDict_.found("startSynced"))
    {
        startSynced_ = readBool(immersedDict_.lookup("startSynced"));
        if (startSynced_)
        {
            Info << "Will try to sync the body with the flow upon creation" << endl;
        }
        else
        {
            Info << "The body will be created as static" << endl;
        }
    }
    else
    {
        Info << "startSynced was not specified, using startSynced: " << startSynced_ << endl;
    }
    
    // just auxiliaries due to MOR
    if (immersedDict_.found("sdBasedLambda"))
    {
        sdBasedLambda_ = readBool(immersedDict_.lookup("sdBasedLambda"));
    }
    if (immersedDict_.found("interfaceSpan"))
    {
        intSpan_ = readScalar(immersedDict_.lookup("interfaceSpan"));
    }
    if (immersedDict_.found("velRelaxFac"))
    {
        velRelaxFac_ = readScalar(immersedDict_.lookup("velRelaxFac"));
    }
    if (immersedDict_.found("refineBuffers"))
    {
        refineBuffers_ = readLabel(immersedDict_.lookup("refineBuffers"));
    }
    if (immersedDict_.found("useInterpolation"))
    {
        useInterpolation_ = readBool(immersedDict_.lookup("useInterpolation"));
    }

    // Set sizes for parallel runs
    surfCells_.setSize(Pstream::nProcs());
    intCells_.setSize(Pstream::nProcs());
    ibPartialVolume_.setSize(Pstream::nProcs());
    wallContactFaces_.setSize(Pstream::nProcs());
    interpolationInfo_.setSize(Pstream::nProcs());
    interpolationVecReqs_.setSize(Pstream::nProcs());
    
    // Initialize bounding points
    boundBox bound(bodySurfMesh_.bounds());
    minBoundPoint_ = bound.min();
    maxBoundPoint_ = bound.max();    
    
    //Initialize list of labels for bounding points based on bodyID
    //This list hold index of both bounding points for each dimension
    //Thanks to this information the update of bounding points positions is more efficient
    for (label coord = 0; coord < 3; coord = coord + 1)
    {
        List<label> labeli(2);
        labeli[0] = 2 * bodyId_;
        labeli[1] = 2 * bodyId_ + 1;
        boundIndList_[coord] = labeli;
    }
    
    Info << "Finished body initialization" << endl;
    Info << "New bodyID: " << bodyId_ << " name: " << stlName_ << " rhoS: " << rhoS_ << " dC: " << dC_ << endl;
}
//---------------------------------------------------------------------------//
//Update immersed body (pre-contact)
void immersedBody::preContactUpdateBodyField
(
    volScalarField& body,
    volVectorField& f
)
{
    preContactUpdateImmersedBody( body, f );
}
//---------------------------------------------------------------------------//
//Create immersed body info
void immersedBody::createImmersedBody(volScalarField& body, volScalarField& refineF, bool createHistory, scalar deltaT)
{
    if (bodyConvex_)
    {
        //~ Info << "Creating body as convex" << endl;
        createImmersedBodyConvex(body, deltaT);
    }
    else
    {
        //~ Info << "Creating body as NOT convex" << endl;
        createImmersedBodyConcave(body, deltaT);
    }
    
    DynamicLabelList zeroList(surfCells_[Pstream::myProcNo()].size(), 0);
    constructRefineField(body, refineF, surfCells_[Pstream::myProcNo()], zeroList);
    
    scalarList charCellSizeL(Pstream::nProcs(),1e4);
    forAll (surfCells_[Pstream::myProcNo()],sCellI)
    {
        label cellI = surfCells_[Pstream::myProcNo()][sCellI];
        charCellSizeL[Pstream::myProcNo()] = min(charCellSizeL[Pstream::myProcNo()],Foam::pow(mesh_.V()[cellI],0.3333));
    }
    
    forAll(charCellSizeL,indl)
    {
        if(charCellSizeL[indl] > 5e3)
            charCellSizeL[indl] = -1.0;
    }
    
    charCellSize_ = gMax(charCellSizeL);
    Info << "Body characteristic cell size: " << charCellSize_ << endl;
    
    // I do NOT want to recompute the body again
    recomputeProjection_ = false;
}
//---------------------------------------------------------------------------//
//Auxiliary for createImmersedBody
labelList immersedBody::getBBoxCellsByOctTree
(
    label cellToCheck,
    bool& insideBB,
    vector& bBoxMin,
    vector& bBoxMax,
    List<DynamicLabelList>& bBoxCells
)
{
    labelList retList;
    
    if (octreeField_[cellToCheck] ==0)
    {
        octreeField_[cellToCheck] = 1;
        vector cCenter = mesh_.C()[cellToCheck];
        label   partCheck(0);
        forAll (bBoxMin,vecI)
        {
            if (cCenter[vecI] >= bBoxMin[vecI] and cCenter[vecI] <= bBoxMax[vecI])
            {
                partCheck++;
            }
        }
        bool cellInside = (partCheck == 3) ? true : false;
        if (cellInside)
        {
            bBoxCells[Pstream::myProcNo()].append(cellToCheck);
            insideBB = true;
        }
        if (not insideBB or cellInside)
        {
            retList = mesh_.cellCells()[cellToCheck];
        }
    }
    return retList;
}
//---------------------------------------------------------------------------//
//Create immersed body for concave body
void immersedBody::createImmersedBodyConcave
(
    volScalarField& body,
    bool createHistory,
    scalar deltaT
)
{        
    // reduce computational domain to the body bounding box
    scalar inflFact(2*sqrt(mesh_.magSf()[0]));
    
    vector expMinBBox = minBoundPoint_ - vector::one*inflFact;
    vector expMaxBBox = maxBoundPoint_ + vector::one*inflFact;
    
    octreeField_ *= 0;
    List<DynamicLabelList> bBoxCells(Pstream::nProcs());
    
    bool isInsideBB(false);
    labelList nextToCheck(1,0);
    label iterCount(0);label iterMax(mesh_.nCells());
    while ((nextToCheck.size() > 0 or not isInsideBB) and iterCount < iterMax)
    {
        iterCount++;        
        DynamicLabelList auxToCheck;
        
        forAll (nextToCheck,cellToCheck)
        {
            auxToCheck.append(
                getBBoxCellsByOctTree(
                    nextToCheck[cellToCheck],
                    isInsideBB,
                    expMinBBox,expMaxBBox,bBoxCells
                )
            );
        }
        nextToCheck = auxToCheck;
    }
                
    // define the search methods
    const triSurface ibTemp( bodySurfMesh_);
    triSurfaceSearch ibTriSurfSearch( ibTemp );
        
    // get cell centers inside the body bounding box
    const pointField& cp = mesh_.C();
    //~ boolList centersInside = ibTriSurfSearch.calcInside(cp);
    const pointField fCp = filterField(cp,bBoxCells[Pstream::myProcNo()]);
    //~ const pointField fCp(cp,bBoxCells[Pstream::myProcNo()]);
    // Note (MI): for some reason, custom filterField is slightly faster
    boolList fCentersInside = ibTriSurfSearch.calcInside(fCp);
    // Note (MI): there is almost no speedup using filtering.
    //            it seems as calcInside does filtering by itself
    //            also, due to the filtering, I need to make a copy
    //            of cp as fCp
    
    // get reference to mesh points
    const pointField& pp = mesh_.points();
    
    // clear old list contents
    intCells_[Pstream::myProcNo()].clear();
    surfCells_[Pstream::myProcNo()].clear();
    
    //Find the processor with most of this IB inside
    ibPartialVolume_[Pstream::myProcNo()] = 0;
    
    //
    vector sDSpan(4.0*(mesh_.bounds().max()-mesh_.bounds().min()));
    //
        
    // first loop, construction of body field and identification of 
    // the number of inside and surface cells
    forAll (bBoxCells[Pstream::myProcNo()],bCellI)                       //go only through bBox
    {
        label cellI(bBoxCells[Pstream::myProcNo()][bCellI]);
        
        //Check if partially or completely inside
        const labelList& vertexLabels = mesh_.cellPoints()[cellI];
        //~ const pointField vertexPoints(pp,vertexLabels);
        const pointField vertexPoints = filterField(pp,vertexLabels);
        boolList vertexesInside = ibTriSurfSearch.calcInside( vertexPoints );
        bool centerInside(fCentersInside[bCellI]);
        scalar rVInSize(0.5/vertexesInside.size());
        // Note: weight of a single vertex in the cell
        
        scalar cBody(0);
        forAll (vertexesInside, verIn)
        {
            if (vertexesInside[verIn]==true)
            {
                cBody  += rVInSize; //fraction of cell covered                
            }
        }
        
        // Note: this is needed for correct definition of internal and
        //       surface cells of the body
        if (centerInside)//consistency with Blais 2016
        {
            cBody+=0.5;
        }
        if (cBody > thrSurf_)
        {
            if (cBody > (1.0-thrSurf_))
            {
                intCells_[Pstream::myProcNo()].append(cellI);
                cellToStartInCreateIB_ = cellI;
            }
            else if (cBody  <= (1.0-thrSurf_))
            {
                surfCells_[Pstream::myProcNo()].append(cellI);
                if (sdBasedLambda_)
                {
                    pointIndexHit pointHit(
                        ibTriSurfSearch.nearest(
                            mesh_.C()[cellI],
                            sDSpan
                        )
                    );
                    scalar signedDist(0.0);
                    if (pointHit.hit())
                    {
                        signedDist = mag(pointHit.hitPoint()-cp[cellI]);
                    }
                    else
                    {// this is due to robustness, not much more can be done here
                        Info << "Missed point in signedDist computation !!" << endl;
                    }
                    if (centerInside)
                    {
                        cBody = 0.5*(Foam::tanh(intSpan_*signedDist/Foam::pow(mesh_.V()[cellI],0.333))+1.0);
                    }
                    else
                    {
                        cBody = 0.5*(-1.0*Foam::tanh(intSpan_*signedDist/Foam::pow(mesh_.V()[cellI],0.333))+1.0);
                    }
                }
            }
            ibPartialVolume_[Pstream::myProcNo()] += 1;
        }
        body[cellI]+= cBody;
        // clip the body field values
        body[cellI] = min(max(0.0,body[cellI]),1.0);
        // Note (MI): max should be useless, min is for overlaps
    }
    
    //gather partial volume from other processors
    Pstream::gatherList(ibPartialVolume_, 0);
    Pstream::scatter(ibPartialVolume_, 0);
    
    for (label i = 0; i < ibPartialVolume_.size(); i++)
    {
        if (ibPartialVolume_[i] == max(ibPartialVolume_))//MI 20201119 OR GMAX?
        {
        //Set owner of the IB which will move this IB
            owner_ = i;
            break;
        }
    }
    
    Info << "body: " << bodyId_ << " owner: " << owner_ << endl;
    
    //refine body as stated in the dictionary
    if(useInterpolation_)
    {
        Info << "Computing interpolation points" << endl;
        calculateInterpolationPoints(body);
    }
    Info << "Computing geometrical properties" << endl;
    calculateGeometricalProperties(body);
    //~ calculateGeometricalProperties2(body, deltaT);
}
//---------------------------------------------------------------------------//
//Create immersed body for convex body
void immersedBody::createImmersedBodyConvex
(
    volScalarField& body,
    bool createHistory,
    scalar deltaT
)
{        
    // Note (MI): I did NOT really modified this function.
    // -> there is almost no speed up filtering cp before calling
    //    calcInside
    // -> I guess that for pp the situation will be similar
    // -> to include filtering actually meens running octree on the
    //    bounding box, which is awfully similar to the actuall body
    //    creation MS: you can test this and modify, if you deem
    //    necessary
    boundBox ibBound(minBoundPoint_,maxBoundPoint_);
    boundBox meshBound(mesh_.bounds());
    bool ibInsideMesh(false);
    forAll(geometricD_,dir)
    {
        if(geometricD_[dir] == 1)
        {
            if(meshBound.contains(minBoundPoint_) || meshBound.contains(maxBoundPoint_))
            {
                ibInsideMesh = true;
                break;
            }
        }
    }
    
    // clear old list contents
    intCells_[Pstream::myProcNo()].clear();
    surfCells_[Pstream::myProcNo()].clear();
    //Find the processor with most of this IB inside
    ibPartialVolume_[Pstream::myProcNo()] = 0;
    octreeField_ *= 0;
    
    if(ibInsideMesh)
    {
        const triSurface ibTemp( bodySurfMesh_);
        triSurfaceSearch ibTriSurfSearch( ibTemp );
        const pointField& pp = mesh_.points();            
        // get the list of cell centroids
        const pointField& cp = mesh_.C();
        
        // see which cell centroids are inside the current body
//         const boolList centersInside = ibTriSurfSearch.calcInside(cp);
//         const boolList pointsInside = ibTriSurfSearch.calcInside(pp);

        bool insideIB(false);
        
        if(cellToStartInCreateIB_ >= octreeField_.size())
            cellToStartInCreateIB_ = 0;
        
        labelList nextToCheck(1,cellToStartInCreateIB_);
        label iterCount(0);label iterMax(mesh_.nCells());
        labelList vertexLabels;
        boolList vertexesInside;
        pointField pointPos;
        bool centerInside;
        DynamicLabelList auxToCheck;
        while (nextToCheck.size() > 0 and iterCount < iterMax)
        {
            iterCount++;        
            auxToCheck.clear();
            
            forAll (nextToCheck,cellToCheck)
            {
                if (octreeField_[nextToCheck[cellToCheck]] == 0)
                {
                    octreeField_[nextToCheck[cellToCheck]] = 1;
                    
                    vertexLabels = mesh_.cellPoints()[nextToCheck[cellToCheck]];
                    pointPos = filterField(pp,vertexLabels);
                    bool cellInsideBB(false);
                    forAll(pointPos,pos)
                    {
                        if(ibBound.contains(pointPos[pos]))
                        {
                            cellInsideBB = true;
                            break;
                        }
                    }
                    
                    if(cellInsideBB)
                    {
//                         vertexesInside = new boolList(pointsInside, vertexLabels);
                        vertexesInside = ibTriSurfSearch.calcInside(pointPos);
                        pointField centerPoint(1,cp[nextToCheck[cellToCheck]]);
                        centerInside = (ibTriSurfSearch.calcInside(centerPoint))[0];
//                         centerInside = centersInside[nextToCheck[cellToCheck]];
                        
                        if(std::any_of(vertexesInside.begin(),vertexesInside.end(),[](bool b){return b;}) || centerInside || !insideIB)
                        {
                            auxToCheck.append(
                                createImmersedBodyByOctTree(
                                    nextToCheck[cellToCheck],
                                    insideIB,
                                    ibPartialVolume_,centerInside,vertexesInside,body
                                )
                            );
                        }
                    }
                    else if(!insideIB)
                    {
                        auxToCheck.append(mesh_.cellCells()[nextToCheck[cellToCheck]]);
                    }
                }
            }
            nextToCheck = auxToCheck;
        }
        
    //     createImmersedBodyByOctTree(cellToStartInCreateIB_, insideIB, ibPartialVolume_, centersInside, pointsInside, body);
        cellToStartInCreateIB_ = min(intCells_[Pstream::myProcNo()]);
    }
    //gather partial volume from other processors
    Pstream::gatherList(ibPartialVolume_, 0);
    Pstream::scatter(ibPartialVolume_, 0);
    for (label i = 0; i < ibPartialVolume_.size(); i++)
    {
        if (ibPartialVolume_[i] == max(ibPartialVolume_))
        {
            //Set owner of the IB which will move this IB
            owner_ = i;
            break;
        }
    }
    
    Info << "body: " << bodyId_ << " owner: " << owner_ << endl;
    
    if(useInterpolation_)
    {
        Info << "Computing interpolation points" << endl;
        calculateInterpolationPoints(body);
    }
    Info << "Computing geometrical properties" << endl;
    calculateGeometricalProperties(body);
    //~ calculateGeometricalProperties2(body, deltaT);
}
//Create immersed body info
labelList immersedBody::createImmersedBodyByOctTree
(
    label cellToCheck,
    bool& insideIB,
    DynamicLabelList& ibPartialVolumei,
    bool& centerInside,
    boolList& vertexesInside,
    volScalarField& body
)
{
    labelList retList;    

    //Check if partially or completely inside
//         const labelList& vertexLabels = mesh_.cellPoints()[cellToCheck];
    //~ const pointField vertexPoints(pp,vertexLabels);
//         boolList vertexesInside(pointsInside, vertexLabels);
//         bool centerInside(centersInside[cellToCheck]);
    scalar rVInSize(0.5/vertexesInside.size());
    // Note: weight of a single vertex in the cell
    
    scalar cBody(0);
    forAll (vertexesInside, verIn)
    {
        if (vertexesInside[verIn]==true)
        {
            cBody  += rVInSize; //fraction of cell covered                
        }
    }
    
    // Note: this is needed for correct definition of internal and
    //       surface cells of the body
    if (centerInside)//consistency with Blais 2016
    {
        cBody+=0.5;
    }
    bool cellInside(false);
    if (cBody > thrSurf_)
    {
        if (cBody > (1.0-thrSurf_))
        {
            intCells_[Pstream::myProcNo()].append(cellToCheck);
        }
        else if (cBody  <= (1.0-thrSurf_))
        {
            surfCells_[Pstream::myProcNo()].append(cellToCheck);
        }
        ibPartialVolume_[Pstream::myProcNo()] += 1;
        cellInside = true;
        insideIB = true;
    }
    body[cellToCheck]+= cBody;
    // clip the body field values
    body[cellToCheck] = min(max(0.0,body[cellToCheck]),1.0);
    // Note (MI): max should be useless, min is for overlaps
    
    if (!insideIB || cellInside)
    {
        retList = mesh_.cellCells()[cellToCheck];
//             labelList cellNb(mesh_.cellCells()[cellToCheck]);//list of neighbours
//             forAll (cellNb,nbCellI)
//             {
//                 createImmersedBodyByOctTree(cellNb[nbCellI], insideIB, ibPartialVolumei, centersInside, pointsInside, body);
//             }
    }
    
    return retList;
}
//---------------------------------------------------------------------------//
//Create immersed body for concave body
void immersedBody::constructRefineField
(
    volScalarField& body,
    volScalarField& refineF,
    DynamicLabelList cellsToIterate,
    DynamicLabelList startLevel 
)
{
    DynamicLabelList cellsToIterateC;
    DynamicLabelList cellsToIterateF;
    
    List<DynamicLabelList> cellsToSendToProcs;
    cellsToSendToProcs.setSize(Pstream::nProcs());
    List<DynamicLabelList> cellsToSendToProcsLevel;
    cellsToSendToProcsLevel.setSize(Pstream::nProcs());
    
    for(label i = 0; i < refineBuffers_; i++)
    {
        forAll(cellsToIterate, cellI)
        {
            if(startLevel[cellI] == i)
                cellsToIterateC.append(cellsToIterate[cellI]);
        }
        
        forAll(cellsToIterateC, cellI)
        {
            labelList cellFaces(mesh_.cells()[cellsToIterateC[cellI]]);
            forAll(cellFaces, faceI)
            {
                if (mesh_.isInternalFace(cellFaces[faceI]))
                {
                    label nCell(mesh_.owner()[cellFaces[faceI]]);
                    if(nCell == cellsToIterateC[cellI])
                    {
                        nCell = mesh_.neighbour()[cellFaces[faceI]];
                    }
                    
                    if(refineF[nCell] == 0)
                    {
                        if(i > 0)
                        {
                            if(body[nCell] < SMALL)
                            {
                                refineF[nCell] = 1;
                                cellsToIterateF.append(nCell);
                            }
                        }
                        else
                        {
                            refineF[nCell] = 1;
                            cellsToIterateF.append(nCell);
                        }
                    }
                }
                else
                {
                    label facePatchId(mesh_.boundaryMesh().whichPatch(cellFaces[faceI]));
                    const polyPatch& cPatch = mesh_.boundaryMesh()[facePatchId];
                    if (cPatch.type() == "processor")
                    {
                        const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(cPatch);
                        if (procPatch.myProcNo() == Pstream::myProcNo())
                        {
                            cellsToSendToProcs[procPatch.neighbProcNo()].append(cPatch.whichFace(cellFaces[faceI]));
                            cellsToSendToProcsLevel[procPatch.neighbProcNo()].append(i+1);
                        }
                        else
                        {
                            cellsToSendToProcs[procPatch.myProcNo()].append(cPatch.whichFace(cellFaces[faceI]));
                            cellsToSendToProcsLevel[procPatch.myProcNo()].append(i+1);
                        }
                    }
                }
            }
        }
        cellsToIterateC = cellsToIterateF;
        cellsToIterateF.clear();
    }
    
    List<DynamicLabelList> facesReceivedFromProcs;
    List<DynamicLabelList> cellsReceivedFromProcsLevel;
    //Send points that are not on this proc to other proc
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UOPstream send(proci, pBufs);
            send << cellsToSendToProcs[proci];
        }
    }
    pBufs.finishedSends();
    
    //Recieve points from other procs
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recv(proci, pBufs);
            DynamicLabelList recList (recv);
            facesReceivedFromProcs.append(recList);
        }
        else
        {
            DynamicLabelList recList;
            facesReceivedFromProcs.append(recList);
        }
    }
    
    pBufs.clear();
    
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            
            UOPstream send(proci, pBufs);
            send << cellsToSendToProcsLevel[proci];
        }
    }
    
    pBufs.finishedSends();
    
    //Recieve points from other procs
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recv(proci, pBufs);
            DynamicLabelList recList (recv);
            cellsReceivedFromProcsLevel.append(recList);
        }
        else
        {
            DynamicLabelList recList;
            cellsReceivedFromProcsLevel.append(recList);
        }
    }
    
    pBufs.clear();
    
    
    DynamicLabelList newCellsToIterate;
    DynamicLabelList newCellsToIterateStartLevel;
    
    //Check if some point from other proc is on this processor
    for (label otherProci = 0; otherProci < facesReceivedFromProcs.size(); otherProci++)
    {
        for (label faceI = 0; faceI < facesReceivedFromProcs[otherProci].size(); faceI++)
        {
            label cellProcI(0);
            forAll (mesh_.boundaryMesh(), patchi)
            {
                const polyPatch& cPatch = mesh_.boundaryMesh()[patchi];
                if (cPatch.type() == "processor")
                {
                    const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(cPatch);
                    if (procPatch.myProcNo() == Pstream::myProcNo() && procPatch.neighbProcNo() == otherProci)
                    {
                        cellProcI = mesh_.faceOwner()[cPatch.start() + facesReceivedFromProcs[otherProci][faceI]];
                        if(refineF[cellProcI] == 0)
                        {
                            newCellsToIterate.append(cellProcI);
                            newCellsToIterateStartLevel.append(cellsReceivedFromProcsLevel[otherProci][faceI]);
                        }
                        break;
                    }
                    else if (procPatch.myProcNo() == otherProci && procPatch.neighbProcNo() == Pstream::myProcNo())
                    {
                        cellProcI = mesh_.faceOwner()[cPatch.start() + facesReceivedFromProcs[otherProci][faceI]];
                        if(refineF[cellProcI] == 0)
                        {
                            newCellsToIterate.append(cellProcI);
                            newCellsToIterateStartLevel.append(cellsReceivedFromProcsLevel[otherProci][faceI]);
                        }
                        break;
                    }
                }
            }
        }
    }
    
    bool contBool(false);
    if(newCellsToIterate.size() > 0)
    {
        contBool = true;
    }

    reduce(contBool, orOp<bool>());
    if(contBool)
    {
        constructRefineField(body, refineF, newCellsToIterate, newCellsToIterateStartLevel);
    }
}
//---------------------------------------------------------------------------//
//detect contact with walls
void immersedBody::detectWallContact(volScalarField& body)
{
    if((M0_-M_) < 0)
        return;
    
    // clear all list required for wall contact
    wallContactFaces_[Pstream::myProcNo()].clear();
    wallContactFacesHelp_.clear();
    
    //Go through all surfCells and check if there is any surfCell whose face is a boundary face        
    forAll (surfCells_[Pstream::myProcNo()],sCellI)
    {
        label cCell(surfCells_[Pstream::myProcNo()][sCellI]);
        
        const labelList& cFaces = mesh_.cells()[cCell];
        
        forAll (cFaces,faceI)
        {
            if (!mesh_.isInternalFace(cFaces[faceI]))
            {
                // Get reference to the patch which is in contact with IB. There is contact only if the patch is marked as a wall
                label facePatchId(-1);
                facePatchId = mesh_.boundaryMesh().whichPatch(cFaces[faceI]);
                const polyPatch& cPatch = mesh_.boundaryMesh()[facePatchId];
                if (cPatch.type()=="wall")
                {          
                    wallContactFaces_[Pstream::myProcNo()].append(cFaces[faceI]);
                    switchWallContact(true);
                    
                    // check the wall for possible contact only if it wasnt already checked
//                     if (findIndex(wallContactFacesHelp_,cFaces[faceI]) == -1)
//                     {
//                         // Append the checked face to the list to prevent cumulative checking
//                         wallContactFacesHelp_.append(cFaces[faceI]);
//                         wallContactFaces_[Pstream::myProcNo()].append(cFaces[faceI]);
//                         switchWallContact(true);
//                         // Estimate penetration depth for current face and if there is intersection with the surface
//                         // Append the face to contactList together with the penetration depth and estimate contactTimeRes
// //                         scalar depth(getPenetrationDepth(cFaces[faceI]));
// //                         if (depth > 0)
// //                         {
// //                             updateContactTimeRes((dC_ * maxDistInDEMloop_)/(mag(Vel_)+SMALL));
// //                             switchWallContact(true);
// //                             wallContactFaces_[Pstream::myProcNo()].append(cFaces[faceI]);
// //                         }
// //                         // Find neighbour faces that are also in contact with the IB
// //                         getContactFaces(cFaces[faceI], body);
//                     }
                }
            }
        }
    }
    
    reduce(isInWallContact_, orOp<bool>());
    
    if(isInWallContact_)
    {
        forAll(intCells_[Pstream::myProcNo()],iCellI)
        {
            label cCell(intCells_[Pstream::myProcNo()][iCellI]);
        
            const labelList& cFaces = mesh_.cells()[cCell];
            
            forAll (cFaces,faceI)
            {
                if (!mesh_.isInternalFace(cFaces[faceI]))
                {
                    // Get reference to the patch which is in contact with IB. There is contact only if the patch is marked as a wall
                    label facePatchId(-1);
                    facePatchId = mesh_.boundaryMesh().whichPatch(cFaces[faceI]);
                    const polyPatch& cPatch = mesh_.boundaryMesh()[facePatchId];
                    if (cPatch.type()=="wall")
                    {          
                        wallContactFaces_[Pstream::myProcNo()].append(cFaces[faceI]);
                    }
                }
            }
            }
    }
    
    
    //Distribute wallContactfaces over processors
    Pstream::gatherList(wallContactFaces_, 0);
    Pstream::scatter(wallContactFaces_, 0);
    
    // Note: - I go through all the surface cells.
    //       - for each surface cell I check the cell faces
    //       - if a face is not internal, I look up the patch in which
    //         it is
    //       - if the face is on "wall" type patch, I check the penetration depth
    //       - If there is intersection include it in current wall contact faces
}
//---------------------------------------------------------------------------//
//Get all surrounding contact faces
void immersedBody::getContactFaces
(
    label faceIndex,
    volScalarField& body
)
{
    // get reference to the patch which is in contact with IB
    label facePatchId = mesh_.boundaryMesh().whichPatch(faceIndex);
    const polyPatch& pp = mesh_.boundaryMesh()[facePatchId];
    
    // Go through neighbour faces and check if they are in contact with IB
    const labelList& faceFacessList = pp.faceFaces()[pp.whichFace(faceIndex)];
    forAll (faceFacessList,faceI)
    {
        // Get global index of the face
        label globalFaceInd = pp.start() + faceFacessList[faceI];
        // check the wall for possible contact only if it wasnt already checked
        if (findIndex(wallContactFacesHelp_, globalFaceInd) == -1)
        {
            // get cell index
            label cellId = pp.faceCells()[faceFacessList[faceI]];
            // Check the current face only of lambda for its owner is not zero
            if (body[cellId] > SMALL)
            {
                wallContactFacesHelp_.append(globalFaceInd);   
                
                scalar depth(getPenetrationDepth(globalFaceInd));
                if (depth > 0)
                {
                    updateContactTimeRes((dC_ * maxDistInDEMloop_)/(mag(Vel_)+SMALL));
                    switchWallContact(true);
                    wallContactFaces_[Pstream::myProcNo()].append(globalFaceInd);
                }
                getContactFaces(globalFaceInd, body);
            }
        }
    }
}

//---------------------------------------------------------------------------//
//Get penetration depth
scalar immersedBody::getPenetrationDepth
(
    label faceInd
)
{
    // Create triSurfaceSearch for the IB surface
    const triSurface ibTemp( bodySurfMesh_);
    triSurfaceSearch ibTriSurfSearch( ibTemp );
    scalar depth(0);
    // Estimate the normal unit vector for the contact 
    vector vecN(mesh_.Sf()[faceInd]/mag(mesh_.Sf()[faceInd]));
    // Move from the contact point further from the IB.
    // note (MS): The nearest search works only from outside of the IB surface!
    vector pointSearch(mesh_.Cf()[faceInd] + dC_ * vecN / 2);
    // Find nearest point on the IB surface
    pointIndexHit ibPointIndexHit = ibTriSurfSearch.nearest(pointSearch, -vecN * dC_);
    //Evaluate the penetration depth only if the nearest point was found
    if (ibPointIndexHit.hit())
    {
        scalar diff(mag(pointSearch - mesh_.Cf()[faceInd]) - mag(pointSearch - ibPointIndexHit.hitPoint()));
        if (diff > 0)
        {
            depth = diff;
        }
    }
    return depth;
}

//---------------------------------------------------------------------------//
//Update immersed body info (pre-contact, ends with contact detection)
void immersedBody::preContactUpdateImmersedBody
(
    volScalarField& body,
    volVectorField& f
)
{
    F_*=0.0;
    T_*=0.0;
    switchWallContact(false);
    
    updateOldMovementVars();
    
    // Assigned variables for potential contact correction
    historyAxis_ = Axis_;
    historyOmega_ = omega_;
    historyVel_ = Vel_;
    historya_ = a_;
    historyAlpha_ = alpha_;
    historyTotalAngle_ = totalAngle_;
    historyBodyPoints_ = bodySurfMesh_.points();
    
    // Note: at the moment, I work only with FLUIDCOUPLING body operation
    //~ updateCoupling(body,f);
    
    detectWallContact(body);
}
//---------------------------------------------------------------------------//
//Update immersed body info (post-contact, contact forces need to be included)
void immersedBody::postContactUpdateImmersedBody
(
    volScalarField& body,
    volVectorField& f
)
{    
    //update Vel_, Axis_ and omega_
    updateCoupling(body,f);
    updateMovement(VelOld_,AxisOld_,omegaOld_);
    
    //update body courant number
    computeBodyCoNumber();
    
    Info << "-- body: " << bodyId_ << " current center of mass position: " << CoM_ << endl;
}
void immersedBody::postContactUpdateImmersedBody(scalar deltaT)
{    
    //update Vel_, Axis_ and omega_
    updateMovement(deltaT);
    
    //update body courant number
    computeBodyCoNumber();
    
    Info << "-- body: " << bodyId_ << " current center of mass position: " << CoM_ << endl;
}
//---------------------------------------------------------------------------//
void immersedBody::updateCoupling
(
    volScalarField& body,
    volVectorField& f
)
{
  const uniformDimensionedVectorField& g =
        mesh_.lookupObject<uniformDimensionedVectorField>("g");

  vector FV(vector::zero);
  vector FG(M_*(1.0-rhoF_.value()/rhoS_.value())*g.value());
  vector TA(vector::zero);
  
//   Pstream::gatherList(surfCells_, 0);
//   forAll(surfCells_, proci)
//   {
//       Info << "Info proc: " << proci << " Computing FV: " << surfCells_[proci] << endl;
//   }
    
  //Calcualate viscous force and torque
  forAll (surfCells_[Pstream::myProcNo()],sCellI)
  {
     label cellI = surfCells_[Pstream::myProcNo()][sCellI];

     // Note (MI): what about the sign of torque?
     //~ F_ -=  f[cellI] * mesh_.V()[cellI] * rhof.value();//viscosity?
     //~ F_ +=  body[cellI]*mesh_.V()[cellI]*(rho_.value()-rhof.value())*g.value();//gravity/buyoancy?
     //~ FV +=  f[cellI] * mesh_.V()[cellI] * rhoF_.value();//viscosity?
     //~ T_ +=  ( (mesh_.C()[cellI]-CoM_) ^ f[cellI] )
                 //~ *mesh_.V()[cellI]* rhoF_.value();
     FV -=  f[cellI]*mesh_.V()[cellI];//viscosity?
     //~ T_ +=  ((mesh_.C()[cellI] - CoM_)^f[cellI])*mesh_.V()[cellI];
     TA +=  ((mesh_.C()[cellI] - CoM_)^f[cellI])*mesh_.V()[cellI];
     // Note (MI): shouldn't this be multiplied by (1-body)?
     //~ FV *=  (scalar(1)-body[cellI]);
     //~ T_ *=  (scalar(1)-body[cellI]);
     // Note (MI): thats what I am trying to do here
  }
//   forAll (intCells_[Pstream::myProcNo()],iCellI)
//   {
//      label cellI = intCells_[Pstream::myProcNo()][iCellI];
// 
//      FV -=  f[cellI]*mesh_.V()[cellI];//viscosity?
//      TA +=  ((mesh_.C()[cellI] - CoM_)^f[cellI])*mesh_.V()[cellI];
//   }
  
  reduce(FV, sumOp<vector>());
  reduce(TA, sumOp<vector>()); 
  
  FV *= rhoF_.value();
  TA *= rhoF_.value();
  
  //~ Info << "-- body " << bodyId_ << " viscous   force: " << FV << endl;
  //~ Info << "-- body " << bodyId_ << " grav/buoy force: " << FG << endl;
  //~ Info << "-- body " << bodyId_ << "          torque: " << T_ << endl;
  
  historyCouplingF_ = FV+FG;
  historyCouplingT_ = TA;
  
  F_ += FV+FG;
  T_ -= TA;
  
  printForcesAndTorques();

}
//---------------------------------------------------------------------------//
//update movement variables of the body
void immersedBody::updateMovement()
{         
    scalar deltaT = mesh_.time().deltaT().value();
    
    updateMovementComp(deltaT,Vel_,Axis_,omega_,1);
    
    //~ Info << "new Vel_  : " << Vel_ << endl;
    //~ Info << "new Axis_ : " << Axis_ << endl;
    //~ Info << "new omega_: " << omega_ << endl;
}
void immersedBody::updateMovement
(
    scalar deltaT
)
{
    updateMovementComp(deltaT,Vel_,Axis_,omega_,1);
    
    //~ Info << "new Vel_  : " << Vel_ << endl;
    //~ Info << "new Axis_ : " << Axis_ << endl;
    //~ Info << "new omega_: " << omega_ << endl;
}
void immersedBody::updateMovement
(
    vector Vel,
    vector Axis,
    scalar omega
)
{
    scalar deltaT = mesh_.time().deltaT().value();
    updateMovementComp(deltaT,Vel,Axis,omega,1);
    
    //~ Info << "new Vel_  : " << Vel_ << endl;
    //~ Info << "new Axis_ : " << Axis_ << endl;
    //~ Info << "new omega_: " << omega_ << endl;
}
void immersedBody::updateMovement
(
    vector Vel,
    vector Axis,
    scalar omega,
    scalar velRelaxFac
)
{
    scalar deltaT = mesh_.time().deltaT().value();
    updateMovementComp(deltaT,Vel,Axis,omega,velRelaxFac);
    
    //~ Info << "new Vel_  : " << Vel_ << endl;
    //~ Info << "new Axis_ : " << Axis_ << endl;
    //~ Info << "new omega_: " << omega_ << endl;
}
void immersedBody::updateMovement
(
    scalar deltaT,
    vector Vel,
    vector Axis,
    scalar omega
)
{
    updateMovementComp(deltaT,Vel,Axis,omega,1);
    
    //~ Info << "new Vel_  : " << Vel_ << endl;
    //~ Info << "new Axis_ : " << Axis_ << endl;
    //~ Info << "new omega_: " << omega_ << endl;
}
void immersedBody::updateMovement
(
    scalar deltaT,
    vector Vel,
    vector Axis,
    scalar omega,
    scalar velRelaxFac
)
{
    updateMovementComp(deltaT,Vel,Axis,omega,velRelaxFac);
    
    //~ Info << "new Vel_  : " << Vel_ << endl;
    //~ Info << "new Axis_ : " << Axis_ << endl;
    //~ Info << "new omega_: " << omega_ << endl;
}
void immersedBody::updateMovementComp
(
    scalar deltaT,
    vector Vel,
    vector Axis,
    scalar omega,
    scalar velRelaxFac
)
{    
    //~ Info << "-- body " << bodyId_ << "!! movement update working with F_ = " << F_
         //~ << " and with T_ = " << T_ << endl;
    
    //~ Info << "old Vel_  : " << Vel_ << endl;
    //~ Info << "old Axis_ : " << Axis_ << endl;
    //~ Info << "old omega_: " << omega_ << endl;
    
    // auxiliary (nested) functions
    // Note: due to the different body operations, it is better to split
    //       the different parts of the movement update in separate steps
    auto updateTranslation = [&]()
    {
        //compute current acceleration (assume constant over timeStep)
        a_  = F_/(M0_+SMALL);
        
        //Update body linear velocity
        Vel_ = (Vel + deltaT*a_) * velRelaxFac;
    };
    
    auto updateRotation = [&]()
    {
        //Update body angular acceleration
        alpha_ = inv(I_) & T_;
        
        //Update body angular velocity
        vector Omega((Axis*omega + deltaT*alpha_) * velRelaxFac);
        
        //Split Omega into Axis_ and omega_
        omega_ = mag(Omega);
        
        if (omega_ < SMALL)
        {
            Axis_ = vector::one;
            if (mesh_.nGeometricD() < 3)
            {
                const vector validDirs = (mesh_.geometricD() + Vector<label>::one)/2;
                Axis_ -= validDirs;
            }
        }
        else
        {
            //~ vector oldAxis = Axis_;
            Axis_ =  Omega/(omega_+SMALL);
            // Note (MI): below is kind of a hack for 2D, in 3D, the
            //            rotation may pose problems - it is still a
            //            little bit of an open issue
            if (mesh_.nGeometricD() < 3)
            {// in 2D, I need to keep only the correct part of the rotation axis
                const vector validDirs = (mesh_.geometricD() + Vector<label>::one)/2;
                Axis_ = cmptMultiply(vector::one-validDirs,Axis_);
            }
            //~ forAll (Axis_,axElI)
            //~ {
                //~ if (mag(Axis_[axElI]) < 1.0e-08) Axis_[axElI] = 0.0;
            //~ }
        }
        Axis_ /= mag(Axis_);
        // Note (MI): I am cutting of small elements of Axis_ (robustness)
    };
    
    auto updateRotationFixedAxis = [&]()
    {
        //Update body angular velocity
        vector Omega((Axis*omega + deltaT * ( inv(I_) & T_ )) * velRelaxFac);
        
        //Split Omega into Axis_ and omega_
        omega_ = mag(Omega);
        
        vector newAxis = Omega/(omega_+SMALL);
        if ((newAxis & Axis_) < 0) Axis_ *= (-1.0);;
        // Note (MI): can this take care of rotation direction?
    };
    
    if (bodyOperation_ == 0 or bodyOperation_ == 3)
    {
        F_ *= 0.0;
        T_ *= 0.0;
        return;
    }
    else if (bodyOperation_ == 1)
    {
        updateRotation();
        
        F_ *= 0.0;
        T_ *= 0.0;
        return;
    }
    else if (bodyOperation_ == 2)
    {
        updateTranslation();
        
        F_ *= 0.0;
        T_ *= 0.0;
        return;
    }
    else if (bodyOperation_ == 4)
    {
        updateRotationFixedAxis();
        
        F_ *= 0.0;
        T_ *= 0.0;
        return;
    }
         
    updateTranslation();
    if (updateTorque_) updateRotation();
    
    F_ *= 0.0;
    T_ *= 0.0;
    
    return;
    
    // Note (MI): after the body movement update, I should discard the
    //            forces used for the update
}

//---------------------------------------------------------------------------//
//Create interpolation points
void immersedBody::calculateInterpolationPoints
(
    volScalarField& body
)
{    
    meshSearch search(mesh_);
    
    // stabilisation for normalisation of the interface normal
    const dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/pow(average(mesh_.V()), 1.0/3.0)
    );
    // Note (MI): this was copied from "interfaceProperties"
    //Create temporary unit surface normals
    vectorField surfNorm(-fvc::grad(body));
    surfNorm /= (mag(surfNorm)+deltaN.value());
    // Note (MI): I am interested only in the normals to the current
    //            body. However, I re-compute surfNorm for all the
    //            body fields.
    
    // get interpolation distances for the body
    //~ scalarField intDists = 1.0*sqrtThree*Foam::pow(mesh_.V(),1.0/3.0);
    // Note (MI): doesn't this work only for hexahedral cells? and even
    //            on such not exactly?
    
    // scale the surfNorm (prepare it for the future computations)
    //~ surfNorm = intDists*surfNorm;

    // clear the old interpolation data
    interpolationInfo_[Pstream::myProcNo()].clear();
    interpolationVecReqs_.clear();
    interpolationVecReqs_.setSize(Pstream::nProcs());
    // Note (MI): OK, in this cycle, I go through all the surface cells
    //            in each cell, I
    //            1. approximate a position of surface point
    //            2. 
    
    //Variables to find in other processors
    List<DynamicLabelList> surfCellRef;
    surfCellRef.setSize(Pstream::nProcs());
    List<DynamicLabelList> orderRef;
    orderRef.setSize(Pstream::nProcs());
    List<DynamicLabelList> surfCellPointRef;
    surfCellPointRef.setSize(Pstream::nProcs());
    List<DynamicLabelList> surfCellLabelRef;
    surfCellLabelRef.setSize(Pstream::nProcs());
    List<DynamicPointList> surfCellLabelRefP;
    surfCellLabelRefP.setSize(Pstream::nProcs());
    
    forAll (surfCells_[Pstream::myProcNo()],cell)
    {               
        //get surface cell label
        label scell = surfCells_[Pstream::myProcNo()][cell];
            
        // estimate intDist
        scalar intDist(0);
        label cellI(scell);
        point surfPoint(mesh_.C()[scell]);
        
        if (sdBasedLambda_)
        {
            scalar minMaxBody(max(min(body[scell],1.0-SMALL),SMALL));
            // NoteMI: this is necessary for robustness - Foam does not like body=1 in atanh
            intDist = Foam::atanh(2.0*minMaxBody - 1.0)*Foam::pow(mesh_.V()[scell],0.333)/intSpan_;
            surfPoint += surfNorm[scell]*intDist;
        }
        else
        {
            scalar dotProd(-GREAT);        
            //~ const labelList& cellNb(mesh_.cellCells()[cellI]);//list of neighbours
            labelList cellNb(mesh_.cellCells()[scell]);//list of neighbours
            forAll (cellNb,nbCellI)
            {
                vector rVec(mesh_.C()[cellNb[nbCellI]] - mesh_.C()[scell]);
                intDist += mag(rVec);
                scalar auxDotProd(rVec & surfNorm[scell]);
                if (auxDotProd > dotProd)
                {
                    dotProd = auxDotProd;
                    cellI   = cellNb[nbCellI];
                }
                //~ Info << "Dot product: " << (rVec & surfNorm[scell]) << endl;
            }
            intDist /= (scalar(cellNb.size())+SMALL);
            surfPoint -= intDist*surfNorm[scell]*(0.5-body[scell]);
        }
        
        //create vector for points and cells and add to main vectors
        DynamicPointList intPoints;
        DynamicLabelList intCells;
        DynamicLabelList procOfCells;
        // Note (MI): these things get extended for each surfCell
        // Note (MI): what happens during mesh update? does the number of
        //            interface cells change? what should I do then?
        // -> destroy old lists and create another ones? (probably)
        // -> I might need a method to call after each mesh update 
        // Note (MI): I should set the size of an element of the global variable
        //            before 
        
        // Note (MI): I need to rewrite intDist and switch it from ds
        //            to res (see HFDIBDEM::interpolateIB)
        intDist = Foam::pow(mesh_.V()[scell],0.333);
        //~ intDist*= 1.5;// move by 1.5 cellScale 
        intDist*= 0.5;// move by 1.5 cellScale 
        
        //Add to list
        intPoints.append(surfPoint);
        vector surfNormToSend(vector::zero);
        if (mag(surfNorm[scell]) > SMALL)
        {
            surfNormToSend = surfNorm[scell]/mag(surfNorm[scell]);
            Tuple2<label,label> helpTup(scell,-1);
            Tuple2<vector,Tuple2<label,label>> startCell(vector::zero,helpTup);
            //~ label cellI = scell;
            //Add other interpolation points
            for (int order=0;order<ORDER;order++)
            {
                //~ surfPoint = surfPoint + surfNorm[scell];
    //             surfPoint = surfPoint + intDist*surfNorm[scell];
    //             cellI = search.findCell(surfPoint,cellI,true);
                startCell = findCellCustom(startCell.first(),startCell.second().first(),startCell.second().second(),surfNormToSend, intDist);
                surfPoint = startCell.first();
                cellI = startCell.second().first();
                
                if (startCell.second().second() == -1)
                {
                    if (startCell.second().first() != -1)
                    {
//                         if (body[cellI] > SMALL)
//                             cellI = -1;
                        if (body[cellI] > SMALL)
                        {
                            order--;
                            continue;
                        }
                    }
                    intPoints.append(surfPoint);
                    intCells.append(cellI);
                    procOfCells.append(Pstream::myProcNo());
                }
                else
                {
                    intPoints.append(surfPoint);
                    intCells.append(-1);
                    procOfCells.append(startCell.second().second());
                    //Info << "Missed intPoint" << endl;
                    // Point is not on this processor. Assign to check if it is on other proc
                    surfCellRef[startCell.second().second()].append(cell);
                    orderRef[startCell.second().second()].append(order);
                    surfCellLabelRef[startCell.second().second()].append(cellI);
                    surfCellPointRef[startCell.second().second()].append(order);
                    surfCellLabelRefP[startCell.second().second()].append(surfPoint);
                }
            }
        }
        else
        {
            for (int order=0;order<ORDER;order++)
            {
                intPoints.append(surfPoint);
                intCells.append(-1);
                procOfCells.append(Pstream::myProcNo());
            }
        }
        //~ Info << intCells << endl;
        // Note (MI): is there a way to estimate ORDER before? (to know the
        //            size of the interpolationPoints_[cell] and the size of 
        //            interpolationCells_[cell]? (I do not like these dynamic
        //            lists (vectors) very much
        // Note (MI): ORDER is hardcoded (at the moment?) ORDER = 2
        // -> size(intPoints) = 2 + 1 (surfPoint) (ORDER + 2)
        // -> size(intCells)  = 2 (ORDER + 1)

        // Note (MI): it seems that there is a problem with List of
        //            std::vector< point > (different field sizes?)
        //            I might do the same trick as in createImmersedBody and
        //            first compute the object sizes and after that allocate
        //            the necessary memory -> Yes, I did it
        
        //Get cells
        //~ Info << "interpolation cells" << endl;
        //~ for (int order=0;order<ORDER;order++)
        //~ {
            //~ //Use findNearestCell()...it is faster
            //~ //Check if outside the domain (if not, then set to -1)
            //~ label cellI = search_.findCell(intPoints[order+1],scell,false);
            //~ intCells.append(cellI);
            //~ Info << cellI << endl;
        //~ }
        //~ // Note (MI): the code should be able to handle outside of domain pts
        //~ Info << intCells << endl;
        // assign to global variables
        immersedBody::interpolationInfo intInfo;
        intInfo.surfCell_ = scell;
        intInfo.intPoints_ = intPoints;
        intInfo.intCells_ = intCells;
        intInfo.procWithIntCells_ = procOfCells;
        intInfo.intVec_.setSize(2);
        interpolationInfo_[Pstream::myProcNo()].append(intInfo);
        
        //~ Info << intPoints << endl;
    }
    //Send points that are not on this proc to other proc
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    PstreamBuffers pBufs2(Pstream::commsTypes::nonBlocking);
    
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UOPstream send(proci, pBufs);
            send << surfCellLabelRef[proci];
            
            UOPstream send2(proci, pBufs2);
            send2 << surfCellPointRef[proci];
        }
    }
    pBufs.finishedSends();
    pBufs2.finishedSends();
    
    List<DynamicLabelList> surfPointRefFromProcs;
    List<DynamicLabelList> surfPointLabelRefFromProcs;
    List<DynamicLabelList> surfCellRefToReturn;
    //Recieve points from other procs
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recv(proci, pBufs);
            DynamicLabelList recSurfPointLabelRef (recv);
            surfPointLabelRefFromProcs.append(recSurfPointLabelRef);
            UIPstream recv2(proci, pBufs2);
            DynamicLabelList recSurfPointRef (recv2);
            surfPointRefFromProcs.append(recSurfPointRef);
        }
        else
        {
            DynamicLabelList recSurfPointLabelRef;
            surfPointLabelRefFromProcs.append(recSurfPointLabelRef);
            DynamicLabelList recSurfPointRef;
            surfPointRefFromProcs.append(recSurfPointRef);
        }
    }
    
    //Check if some point from other proc is on this processor
    for (label otherProci = 0; otherProci < surfPointRefFromProcs.size(); otherProci++)
    {
        DynamicLabelList surfCellRefFromProc;
        for (label surfPointI = 0; surfPointI < surfPointRefFromProcs[otherProci].size(); surfPointI++)
        {
            label cellProcI(0);
            forAll (mesh_.boundaryMesh(), patchi)
            {
                const polyPatch& cPatch = mesh_.boundaryMesh()[patchi];
                if (cPatch.type() == "processor")
                {
                    const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(cPatch);
                    if (procPatch.myProcNo() == Pstream::myProcNo() && procPatch.neighbProcNo() == otherProci)
                    {
                        cellProcI = mesh_.faceOwner()[cPatch.start() + surfPointLabelRefFromProcs[otherProci][surfPointI]];
                        if (surfPointRefFromProcs[otherProci][surfPointI] == 1)
                        {
                            vector sfUnit(mesh_.Sf()[cPatch.start() + surfPointLabelRefFromProcs[otherProci][surfPointI]]/mag(mesh_.Sf()[cPatch.start() + surfPointLabelRefFromProcs[otherProci][surfPointI]]));
                            
                            label bestCell(0);
    
                            scalar dotProd(-GREAT);
                            labelList cellNb(mesh_.cellCells()[cellProcI]);//list of neighbours
                            forAll (cellNb,nbCellI)
                            {
                                vector vecI(mesh_.C()[cellNb[nbCellI]] - mesh_.C()[cellProcI]);
                                vecI /= mag(vecI);
                                scalar auxDotProd(vecI & sfUnit);
                                if (auxDotProd > dotProd)
                                {
                                    dotProd = auxDotProd;
                                    bestCell = cellNb[nbCellI];
                                }
                            }
                            cellProcI = bestCell;
                        }
                        break;
                    }
                    else if (procPatch.myProcNo() == otherProci && procPatch.neighbProcNo() == Pstream::myProcNo())
                    {
                        cellProcI = mesh_.faceOwner()[cPatch.start() + surfPointLabelRefFromProcs[otherProci][surfPointI]];
                        if (surfPointRefFromProcs[otherProci][surfPointI] == 1)
                        {
                            vector sfUnit(mesh_.Sf()[cPatch.start() + surfPointLabelRefFromProcs[otherProci][surfPointI]]/mag(mesh_.Sf()[cPatch.start() + surfPointLabelRefFromProcs[otherProci][surfPointI]]));
                            
                            label bestCell(0);
    
                            scalar dotProd(-GREAT);
                            labelList cellNb(mesh_.cellCells()[cellProcI]);//list of neighbours
                            forAll (cellNb,nbCellI)
                            {
                                vector vecI(mesh_.C()[cellNb[nbCellI]] - mesh_.C()[cellProcI]);
                                vecI /= mag(vecI);
                                scalar auxDotProd(vecI & sfUnit);
                                if (auxDotProd > dotProd)
                                {
                                    dotProd = auxDotProd;
                                    bestCell = cellNb[nbCellI];
                                }
                            }
                            cellProcI = bestCell;
                        }
                        break;
                    }
                }
            }
            
            if (cellProcI >= 0)
            {
                if (body[cellProcI] > SMALL)
                    cellProcI = -1;
            }
            surfCellRefFromProc.append(cellProcI);
        }
        //Assing found points to send them back
        surfCellRefToReturn.append(surfCellRefFromProc);
    }
    
    pBufs.clear();
    //Send points back
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UOPstream send(proci, pBufs);
            send << surfCellRefToReturn[proci];
        }
    }
    
    pBufs.finishedSends();
    
    List<DynamicLabelList> surfCellRefRcv;
    //Receive points from other processors
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recv(proci, pBufs);
            DynamicLabelList recsurfCellRef (recv);
            surfCellRefRcv.append(recsurfCellRef);
        }
        else
        {
            DynamicLabelList recsurfCellRef;
            surfCellRefRcv.append(recsurfCellRef);
        }
    }

    //Check if some points were found and assign them to info
    for (label otherProci = 0; otherProci < surfCellRefRcv.size(); otherProci++)
    {
        for (label intCellI = 0; intCellI < surfCellRefRcv[otherProci].size(); intCellI++)
        {
            label cellProcI(surfCellRefRcv[otherProci][intCellI]);
            if (cellProcI > -1)
            {
                interpolationInfo_[Pstream::myProcNo()][surfCellRef[otherProci][intCellI]].intCells_[orderRef[otherProci][intCellI]] = cellProcI;
            }
        }
    }
    //Decide which order should be used
    for (label infoI = 0; infoI < interpolationInfo_[Pstream::myProcNo()].size(); infoI++)
    {
        List<bool> allowedOrder;
        allowedOrder.setSize(ORDER);
        // Note (MI): Ok, i will create a list of bools of size ORDER
        
        for (int intPoint=0;intPoint<ORDER;intPoint++)
        {
            if (interpolationInfo_[Pstream::myProcNo()][infoI].intCells_[intPoint] == -1)
            {
                allowedOrder[intPoint] = false;
            }
            else
            {
                allowedOrder[intPoint] = true;
            }
        }
        
        interpolationInfo_[Pstream::myProcNo()][infoI].order_ = 2;
        if ( allowedOrder[1] == false)
        {
            interpolationInfo_[Pstream::myProcNo()][infoI].order_ = 1;
        }
        
        //Check if first order is possible
        if ( allowedOrder[0] == false)
        {
            interpolationInfo_[Pstream::myProcNo()][infoI].order_ = 0;
        }
        
        if (interpolationInfo_[Pstream::myProcNo()][infoI].order_ == 2)
        {
            if (interpolationInfo_[Pstream::myProcNo()][infoI].intCells_[0] == interpolationInfo_[Pstream::myProcNo()][infoI].intCells_[1])
                interpolationInfo_[Pstream::myProcNo()][infoI].order_ = 1;
        }
        //If some cell on other proc prepare requst for interpolation
        switch(interpolationInfo_[Pstream::myProcNo()][infoI].order_)
        {
            case 1:
                if (interpolationInfo_[Pstream::myProcNo()][infoI].procWithIntCells_[0] != Pstream::myProcNo())
                {
                    immersedBody::intVecRequest vecReq;
                    vecReq.requestLabel_ = infoI;
                    vecReq.vecLabel_ = 0;
                    vecReq.intPoint_ = interpolationInfo_[Pstream::myProcNo()][infoI].intPoints_[1];
                    vecReq.intCell_ = interpolationInfo_[Pstream::myProcNo()][infoI].intCells_[0];
                    interpolationVecReqs_[interpolationInfo_[Pstream::myProcNo()][infoI].procWithIntCells_[0]].append(vecReq);
                }
                break;
            case 2:
                if (interpolationInfo_[Pstream::myProcNo()][infoI].procWithIntCells_[0] != Pstream::myProcNo())
                {
                    immersedBody::intVecRequest vecReq;
                    vecReq.requestLabel_ = infoI;
                    vecReq.vecLabel_ = 0;
                    vecReq.intPoint_ = interpolationInfo_[Pstream::myProcNo()][infoI].intPoints_[1];
                    vecReq.intCell_ = interpolationInfo_[Pstream::myProcNo()][infoI].intCells_[0];
                    interpolationVecReqs_[interpolationInfo_[Pstream::myProcNo()][infoI].procWithIntCells_[0]].append(vecReq);
                }
                
                if (interpolationInfo_[Pstream::myProcNo()][infoI].procWithIntCells_[1] != Pstream::myProcNo())
                {
                    immersedBody::intVecRequest vecReq;
                    vecReq.requestLabel_ = infoI;
                    vecReq.vecLabel_ = 1;
                    vecReq.intPoint_ = interpolationInfo_[Pstream::myProcNo()][infoI].intPoints_[2];
                    vecReq.intCell_ = interpolationInfo_[Pstream::myProcNo()][infoI].intCells_[1];
                    interpolationVecReqs_[interpolationInfo_[Pstream::myProcNo()][infoI].procWithIntCells_[1]].append(vecReq);
                }
                break;
        }
    }
}
//---------------------------------------------------------------------------//
// Custom function to find cell containing point
Tuple2<vector,Tuple2<label,label>> immersedBody::findCellCustom
(
    vector& prevPoint,
    label& startCell,
    label& startProc,
    vector& gradToBody,
    scalar& intDist
)
{
    if (startProc != -1)
    {
        
        Tuple2<label,label> helpTup(startCell,startProc);
        Tuple2<vector,Tuple2<label,label>> tupleToReturn(prevPoint + intDist*gradToBody,helpTup);
        return tupleToReturn;
    }
    else if (startCell == -1)
    {
        Tuple2<label,label> helpTup(-1,-1);
        Tuple2<vector,Tuple2<label,label>> tupleToReturn(prevPoint + intDist*gradToBody,helpTup);
        return tupleToReturn;
    }
    
    labelList cellFaces(mesh_.cells()[startCell]);//list of face indices
    label bestFace(0);
    
    scalar dotProd(-GREAT);
    forAll (cellFaces,faceI)
    {
        vector vecI(mesh_.Cf()[cellFaces[faceI]] - mesh_.C()[startCell]);
        vecI /= mag(vecI);
        scalar auxDotProd(vecI & gradToBody);
        if (auxDotProd > dotProd)
        {
            dotProd = auxDotProd;
            bestFace = cellFaces[faceI];
        }
    }
    labelList cellPoints(mesh_.faces()[bestFace]);//list of vertex indicies
    DynamicList<Tuple2<vector,Tuple2<label,label>>> pointsToCheck;
    forAll (cellPoints, pointI)
    {
        labelList pointFaces(mesh_.pointFaces()[cellPoints[pointI]]);
        forAll (pointFaces, faceI)
        {
            if (mesh_.isInternalFace(pointFaces[faceI]))
            {
                label owner(mesh_.owner()[pointFaces[faceI]]);
                label neighbour(mesh_.neighbour()[pointFaces[faceI]]);
                if (owner != startCell)
                {
                    bool add(true);
                    forAll (pointsToCheck, i)
                    {
                        if (pointsToCheck[i].second().first() == owner)
                        {
                            add = false;
                            break;
                        }
                    }
                    if (add)
                    {
                        Tuple2<label,label> helpTup(owner,-1);
                        Tuple2<vector,Tuple2<label,label>> tupleToAdd(mesh_.C()[owner],helpTup);
                        pointsToCheck.append(tupleToAdd);
                    }
                }
                if (neighbour != startCell)
                {
                    bool add(true);
                    forAll (pointsToCheck, i)
                    {
                        if (pointsToCheck[i].second().first() == neighbour)
                        {
                            add = false;
                            break;
                        }
                    }
                    if (add)
                    {
                        Tuple2<label,label> helpTup(neighbour,-1);
                        Tuple2<vector,Tuple2<label,label>> tupleToAdd(mesh_.C()[neighbour],helpTup);
                        pointsToCheck.append(tupleToAdd);
                    }
                }
            }
            else
            {
                label owner(mesh_.faceOwner()[pointFaces[faceI]]);
                vector distToFace(mesh_.Cf()[pointFaces[faceI]] - mesh_.C()[owner]);
                vector sfUnit(mesh_.Sf()[pointFaces[faceI]]/mag(mesh_.Sf()[pointFaces[faceI]]));
//                 vector sfUnit(mesh_.Sf()[bestFace]/mag(mesh_.Sf()[bestFace]));
                vector pointToAppend(mesh_.C()[owner] + mag(distToFace)*sfUnit);
                
                label facePatchId(mesh_.boundaryMesh().whichPatch(pointFaces[faceI]));
                const polyPatch& cPatch = mesh_.boundaryMesh()[facePatchId];
                if (cPatch.type() == "processor")
                {
                    const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(cPatch);
                    if (procPatch.myProcNo() == Pstream::myProcNo())
                    {
                        Tuple2<label,label> helpTup(cPatch.whichFace(pointFaces[faceI]),procPatch.neighbProcNo());
                        Tuple2<vector,Tuple2<label,label>> tupleToAdd(pointToAppend,helpTup);
                        pointsToCheck.append(tupleToAdd);
                    }
                    else
                    {
                        Tuple2<label,label> helpTup(cPatch.whichFace(pointFaces[faceI]),procPatch.myProcNo());
                        Tuple2<vector,Tuple2<label,label>> tupleToAdd(pointToAppend,helpTup);
                        pointsToCheck.append(tupleToAdd);
                    }
                }
                else
                {
                    bool add(true);
                    forAll (pointsToCheck, i)
                    {
                        if (pointsToCheck[i].first() == pointToAppend)
                        {
                            add = false;
                            break;
                        }
                    }
                    if (add)
                    {
                        Tuple2<label,label> helpTup(-1,-1);
                        Tuple2<vector,Tuple2<label,label>> tupleToAdd(pointToAppend,helpTup);
                        pointsToCheck.append(tupleToAdd);
                    }
                }
            }
        }
    }
    
    dotProd = -GREAT;
    Tuple2<label,label> helpTup(-1,-1);
    Tuple2<vector,Tuple2<label,label>> tupleToReturn(vector::zero,helpTup);
    forAll (pointsToCheck,pointI)
    {
        vector vecI(pointsToCheck[pointI].first() - mesh_.C()[startCell]);
        vecI /= (mag(vecI)+SMALL);
        scalar auxDotProd(vecI & gradToBody);
        if (auxDotProd > dotProd)
        {
            dotProd = auxDotProd;
            tupleToReturn   = pointsToCheck[pointI];
        }
    }

    return tupleToReturn;
}
//---------------------------------------------------------------------------//
// Auxiliary function to add to M_ and I_
vector immersedBody::addToMAndI
(
    volScalarField& body,
    DynamicLabelList& labelCellLst,
    vector tmpCom
)
{
    //Get density
    dimensionedScalar rho=dimensionedScalar(immersedDict_.lookup("rho"));
    
    forAll (labelCellLst,cell)
    {
        label cellI = labelCellLst[cell];
        M_ += body[cellI] * rho.value() * mesh_.V()[cellI];
        tmpCom  += body[cellI] * rho.value() * mesh_.V()[cellI] * mesh_.C()[cellI];
        
        I_.xx() += body[cellI]*rho.value()*mesh_.V()[cellI]
                * (
                     mesh_.C()[cellI].y()*mesh_.C()[cellI].y()
                   + mesh_.C()[cellI].z()*mesh_.C()[cellI].z()
               );
        
        I_.yy() += body[cellI]*rho.value()*mesh_.V()[cellI]
                * (
                     mesh_.C()[cellI].x()*mesh_.C()[cellI].x()
                   + mesh_.C()[cellI].z()*mesh_.C()[cellI].z()
               );
        I_.zz() += body[cellI]*rho.value()*mesh_.V()[cellI]
                * (
                     mesh_.C()[cellI].y()*mesh_.C()[cellI].y()
                   + mesh_.C()[cellI].x()*mesh_.C()[cellI].x()
               );
        
        I_.xy() -= body[cellI]*rho.value()*mesh_.V()[cellI]
                * (
                     mesh_.C()[cellI].x()*mesh_.C()[cellI].y()
                 );
        
        I_.xz() -= body[cellI]*rho.value()*mesh_.V()[cellI]
                * (
                     mesh_.C()[cellI].x()*mesh_.C()[cellI].z()
                 );
        
        I_.yz() -= body[cellI]*rho.value()*mesh_.V()[cellI]
                * (
                     mesh_.C()[cellI].y()*mesh_.C()[cellI].z()
                 );
    }
    
    return tmpCom;
}
// Note (MI): tmpCom vector is shared between surfCells and intCells
//            -> I should make it an input argument and return it
//            from the function (20190206)

//---------------------------------------------------------------------------//
// auxiliary function just to compute the center of mass and the total
// mass of the body
vector immersedBody::addToM
(
    volScalarField& body,
    DynamicLabelList& labelCellLst,
    vector tmpCoM
)
{
    forAll (labelCellLst,cell)
    {
        label cellI  = labelCellLst[cell];
        scalar Mi    = body[cellI] * rhoS_.value() * mesh_.V()[cellI];
        M_          += Mi;
        tmpCoM      += Mi * mesh_.C()[cellI];
    }
    return tmpCoM;
}
//---------------------------------------------------------------------------//
// auxiliary function just to inertial tensor of the body
void immersedBody::addToI
(
    volScalarField& body,
    DynamicLabelList& labelCellLst
)
{
    forAll (labelCellLst,cell)
    {
        label cellI  = labelCellLst[cell];
        
        scalar Mi    = body[cellI] * rhoS_.value() * mesh_.V()[cellI];
        //~ const scalar& xLoc = mesh_.C()[cellI].x();
        //~ const scalar& yLoc = mesh_.C()[cellI].y();
        //~ const scalar& zLoc = mesh_.C()[cellI].z();
        scalar xLoc = mesh_.C()[cellI].x() - CoM_.x();
        scalar yLoc = mesh_.C()[cellI].y() - CoM_.y();
        scalar zLoc = mesh_.C()[cellI].z() - CoM_.z();
        // Note (MI): this is just for better code readability
        
        I_.xx() += Mi*(yLoc*yLoc + zLoc*zLoc);
        I_.yy() += Mi*(xLoc*xLoc + zLoc*zLoc);
        I_.zz() += Mi*(xLoc*xLoc + yLoc*yLoc);
        
        I_.xy() -= Mi*(xLoc*yLoc);
        I_.xz() -= Mi*(xLoc*zLoc);
        I_.yz() -= Mi*(yLoc*zLoc);
    }
}

void immersedBody::updateCoM(volScalarField& body, scalar deltaT)
{
    // go through cells and compute local mass and tmpCoM
    vector tmpCoM(vector::zero);
    //~ M_ = 0;
    tmpCoM = addToM(body,surfCells_[Pstream::myProcNo()],tmpCoM);
    tmpCoM = addToM(body,intCells_[Pstream::myProcNo()],tmpCoM);
    
    // collect from processors
    reduce(M_, sumOp<scalar>());
    reduce(tmpCoM,  sumOp<vector>());
    
    // final computation of the body center of mass
    CoM_ = tmpCoM / (M_+SMALL);
    
    //~ Info << "Recomputed M" << endl;
}

void immersedBody::updateI(volScalarField& body)
{
    I_      = symmTensor::zero;
    
    addToI(body,surfCells_[Pstream::myProcNo()]);
    addToI(body,intCells_[Pstream::myProcNo()]);
    
    reduce(I_,  sumOp<symmTensor>());
}

void immersedBody::calculateGeometricalProperties( volScalarField& body )
{
    //Evaluate center of mass
    M_      = scalar(0);
    CoM_    = vector::zero;
    I_      = symmTensor::zero;
    vector tmpCom(vector::zero);
    // Note (MI): M_ ... body mass

    // Note (MI): why is the below loop coded this way?
    //            - forAll intCells and surfCells and cellI is updated
    //            based on which index I am using at the moment?
    // Q:         could I make the coding cleaner?
    // Q (MI)   : is there a common body field for all the immersed
    //            bodies? (I guess it is) -> yes, it seems that the
    //            body field is common to all the bodies and intCells_
    //            and surfCells_ are proper to each body
    // =>         I need to go only through intCells and surfCells in here
    //            and going conditionally through body field is not OK
    tmpCom = addToMAndI(body,surfCells_[Pstream::myProcNo()],tmpCom);
    tmpCom = addToMAndI(body,intCells_[Pstream::myProcNo()],tmpCom);
    // Note (MI): clean up of the body geometricalProperties computation
    //            fixed the issue with improper rotation (at least in
    //            2D on static mesh), which I probably broke while
    //            tweaking the code
    // Note (MI): the rotation seems to be working on 2D dynamically
    //            refined meshes as well
    // Note (MI): velocity field with rotation is still broken (20190201)
    // =>         the imposed force is for some reason applied inside
    //            the body
    
  
    // Note (MI): I splitted the original combined for loop (see v8) in
    //            two distinct loops - one for internal cells and one
    //            for surface cells. for the reasoning behind, go through
    //            the notes in the function createImmersedBody
    // Note (MI): it would be good to put addToM_ and addToI_ in separate
    //            functions - done (20190201)
    
    //Collect from processors
    reduce(M_, sumOp<scalar>());
    reduce(tmpCom,  sumOp<vector>());
    reduce(I_,  sumOp<symmTensor>());
    
//     CoM_ = tmpCom / (M_+SMALL);
    // Note (MI): added +SMALL for robustness, effect on result not tested
    
    //~ if (isFirstUpdate_)
    //~ {
        //~ forAll (surfCells_[Pstream::myProcNo()],sCellI)
        //~ {
            //~ label cellI = surfCells_[Pstream::myProcNo()][sCellI];
            //~ dC_ = max(dC_,mag(CoM_-mesh_.C()[cellI]));
        //~ }
        //~ Info << M_ << " " << M0_ << endl;
        //~ M0_ = M_;
        //~ Info << "asigning M0_" << endl;
    //~ }
    
//     M_ = M0_;
    
    CoM_ = tmpCom / (M_+SMALL);
    
    Info << "-- body " << bodyId_ << " current center of mass position: " << CoM_ << endl;
    // Note (MI): compute characteristic diameter as maximum distance
    //            of surface cell to the body center of mass
    
    //~ //if only 1% of the initial particle mass remains in the domain, switch it off
    //~ Info << "-- body " << bodyId_ << " current M/M0: " << M_/M0_ << endl;
    //~ if (M_/(M0_+SMALL) < 1e-2) switchActiveOff(body);
}

void immersedBody::calculateGeometricalProperties2( volScalarField& body, scalar deltaT )
{
    //Evaluate center of mass
    updateCoM(body, deltaT);
    
    //Evaluate tensor of inertia
    updateI(body);
    
    //~ if (isFirstUpdate_)
    //~ {
        //~ forAll (surfCells_[Pstream::myProcNo()],sCellI)
        //~ {
            //~ label cellI = surfCells_[Pstream::myProcNo()][sCellI];
            //~ dC_ = max(dC_,mag(CoM_-mesh_.C()[cellI]));
        //~ }
        //~ M0_ = M_;
        //~ Info << "Recomputed M0" << endl;
    //~ }
    Info << "-- body " << bodyId_ << " current center of mass position: " << CoM_ << endl;
    // Note (MI): compute characteristic diameter as maximum distance
    //            of surface cell to the body center of mass
    
    //~ if (isNewBody_)
    //~ {
        //~ M0_ = M_;
        //~ isNewBody_ = false;
        //~ Info << "Recomputed M0" << endl;
    //~ }

    //~ Info << "-- body " << bodyId_ << " current M/M0: " << M_/M0_ << endl;
    //~ //if only 1% of the initial particle mass remains in the domain, switch it off
    //~ if (M_/(M0_+SMALL) < 1e-2) switchActiveOff(body);
}
//---------------------------------------------------------------------------//
//Move immersed body according to body operation
void immersedBody::moveImmersedBody
(
    scalar deltaT
)
{
    if (bodyOperation_ == 0) return;
    
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    if (owner_ == Pstream::myProcNo())
    {
        if (mag(deltaT + 1.0) < SMALL) deltaT = mesh_.time().deltaT().value();

        //Rotation angle
        //~ scalar angle     = omega_*deltaT;//Forward Euler
        scalar angle     = omega_*deltaT - 0.5*mag(alpha_)*deltaT*deltaT;
        //~ vector transIncr = (Vel_-a_*deltaT)*deltaT + 0.5*a_*deltaT*deltaT;
        //~ vector transIncr = Vel_*deltaT;//Forward Euler
        vector transIncr = Vel_*deltaT - 0.5*a_*deltaT*deltaT;
        // Note: (Vel_-a_*deltaT) => Vel_OLD
        // Note: this is midpoint 2nd order time integration
        // Note (MI): expanded the parantheses and simplified the relation
        scalar rotIncr(0.5*omega_*deltaT*dC_);
        
        dS_ += Foam::sqrt(
            Foam::magSqr(transIncr) + Foam::pow(rotIncr,2.0)
        );
        
//         Pout << "current ds/charCellSize = " << dS_/(charCellSize_+SMALL) << endl;
        //~ if (dS_ > 5.0*charCellSize_)
        if (dS_ > 0.0*charCellSize_)
        //~ if (dS_ > 0.0*charCellSize_)
        {
            recomputeProjection_ = true;
            dS_  = 0.0;
//             Pout << "!! The body will be re-projected into the mesh" << endl;
        }
        else
        {
            CoM_ += transIncr;
        }
        
        // standard rotation
        //~ tensor rotMatrix(Foam::cos(angle)*tensor::I);
        //~ rotMatrix += Foam::sin(angle)*tensor(0,-Axis_.z(),Axis_.y(),Axis_.z(),0,-Axis_.x(),-Axis_.y(),Axis_.x(),0);
        //~ rotMatrix += (1.0-Foam::cos(angle))*(Axis_ * Axis_);
        tensor rotMatrix(Foam::cos(angle)*tensor::I);
        rotMatrix += Foam::sin(angle)*tensor(
            0.0,      -Axis_.z(),  Axis_.y(),
            Axis_.z(), 0.0,       -Axis_.x(),
            -Axis_.y(), Axis_.x(),  0.0
        ); // original
        //~ rotMatrix += Foam::sin(angle)*tensor(
            //~ 0.0,       Axis_.z(), -Axis_.y(),
            //~ -Axis_.z(), 0.0,        Axis_.x(),
            //~ Axis_.y(),-Axis_.x(),  0.0
        //~ ); // transpose
        rotMatrix += (1.0-Foam::cos(angle))*(Axis_ * Axis_);
        
        //update total rotation matrix
        totRotMatrix_ = rotMatrix & totRotMatrix_;
        vector eulerAngles;
        scalar sy = Foam::sqrt(totRotMatrix_.xx()*totRotMatrix_.xx() + totRotMatrix_.yy()*totRotMatrix_.yy());
        if (sy > SMALL)
        {
            eulerAngles.x() = Foam::atan2(totRotMatrix_.zy(),totRotMatrix_.zz());
            eulerAngles.y() = Foam::atan2(-totRotMatrix_.zx(),sy);
            eulerAngles.z() = Foam::atan2(totRotMatrix_.yx(),totRotMatrix_.xx());
        }
        else
        {
            eulerAngles.x() = Foam::atan2(-totRotMatrix_.yz(),totRotMatrix_.yy());
            eulerAngles.y() = Foam::atan2(-totRotMatrix_.zx(),sy);
            eulerAngles.z() = 0.0;
        }
        
        //~ vector locAxis(Axis_+CoM_);
        //~ locAxis /= mag(locAxis);
        //~ tensor rotMatrix(Foam::cos(angle)*tensor::I);
        //~ rotMatrix += Foam::sin(angle)*tensor(0,-locAxis.z(),locAxis.y(),locAxis.z(),0,-locAxis.x(),-locAxis.y(),locAxis.x(),0);
        //~ rotMatrix += (1.0-Foam::cos(angle))*(locAxis * locAxis);
        
        // Rodrigues rotation
        //~ vector angleVec  = angle*Axis_;
        //~ tensor skewAngleVec = tensor(0.0,angleVec.z(),-angleVec.y(),-angleVec.z(),0.0,angleVec.x(),angleVec.y(),-angleVec.x(),0.0);
        //~ tensor  rotMatrix(tensor::I);
        //~ rotMatrix += 4.0/(1.0+Foam::pow(angle,2.0))*(skewAngleVec+skewAngleVec & skewAngleVec);
        //~ totalAngle_ = 4.0/(4.0 - (totalAngle_ & angleVec)) * (totalAngle_ + angleVec - 0.5*(totalAngle_ ^ angleVec));
        
        //~ Info << "-- body " << bodyId_ << " rotation incr:    " << angleVec << endl;
        //~ Info << "-- body " << bodyId_ << " total rotation:   " << totalAngle_ << endl;
        
        pointField bodyPoints (bodySurfMesh_.points());
        
        Info << "-- body " << bodyId_ << " linear velocity      :  " << Vel_ << endl;
        Info << "-- body " << bodyId_ << " angluar velocity     : " << omega_ << endl;
        Info << "-- body " << bodyId_ << " axis of rotation     : " << Axis_ << endl;
        Info << "-- body " << bodyId_ << " total rotation matrix: " << totRotMatrix_ << endl;
        Info << "-- body " << bodyId_ << " total euler angles   : " << eulerAngles << endl;
        //~ Info << "-- body " << bodyId_ << " loc. axis of rotation: " << locAxis << endl;
                
            
        //Move points
        forAll (bodyPoints,p)
        {
            //vector normV    = (bodyPoints[p]-CoM_)^Axis_;
            //scalar axisMod  = ((bodyPoints[p]-CoM_)&Axis_);
            //scalar magDist  = mag((bodyPoints[p]-CoM_) - Axis_*axisMod/(mag(Axis_)+SMALL));
            
            //Move in tangential direction
            //~ bodyPoints[p] = bodyPoints[p] + angle*normV;
            // Note (MI): I believe that Municchi has incorect sign of normV
            //            (clock-wise direction, while standard is 
            //            counter-clockwise (for omega > 0))
            
            // move the point in a way that CoM_ becomes origin
            bodyPoints[p] -= CoM_;
            
            //Rotate the point
            //~ Info << bodyPoints[p] << endl;
            bodyPoints[p] = rotMatrix & bodyPoints[p];
            //~ bodyPoints[p] = Axis_*(Axis_ & bodyPoints[p]) + Foam::cos(angle)*((Axis_ ^ bodyPoints[p]) ^ Axis_) + Foam::sin(angle)*(Axis_ ^ bodyPoints[p]);
            //~ Info << bodyPoints[p] << endl;
            //~ Info << "----" << endl;
            
            // move the point back in its original place (after rotation)
            bodyPoints[p] += CoM_;
            
            //Rescale
            //~ bodyPoints[p] = (bodyPoints[p] - Axis_*((bodyPoints[p]-CoM_)&Axis_))
                        //~ * magDist /(SMALL+mag(
                                //~ bodyPoints[p] - Axis_*((bodyPoints[p]-CoM_)&Axis_)
                            //~ )) +
                        //~ Axis_*((bodyPoints[p]-CoM_)&Axis_);
            //~ bodyPoints[p] = 
            //~ (
                //~ (bodyPoints[p] - Axis_*axisMod)/(mag(bodyPoints[p] - Axis_*axisMod) + SMALL)*magDist
                //~ + Axis_*axisMod
            //~ );
            // Note (MI): I turned this off on June 21 2019. What is the point of
            //            rescaling the body (if it is a rigid one)
            // Note (MI): might this be connected to potential distortion of
            //            body with rotation? without that the body gets
            //            distorted (or it seems so), with it the solver
            //            crashes with segFault
            
            //Translate point
            bodyPoints[p] += transIncr;
        
        }
        //~ Info << bodyPoints << endl;
        // Note (MI): there has to be something rotten in here. the body
        //            movement on dynamic mesh is not clear
        //            !! test it on static mesh to identify the root cause !!
        // =>         it seems that the problem was in the computation of body
        //            geometrical properties (hopefully corrected on 20190201)
        
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            UOPstream send(proci, pBufs);
            send << bodyPoints;
        }
    }
    pBufs.finishedSends();
    //Move body to points calculated by owner_
    UIPstream recv(owner_, pBufs);
    pointField bodyPoints2 (recv);
    reduce(recomputeProjection_, orOp<bool>());
    
    //move mesh
    bodySurfMesh_.movePoints(bodyPoints2);
    
    // Update bounds of the body
    boundBox bound(bodyPoints2);
    minBoundPoint_ = bound.min();
    maxBoundPoint_ = bound.max();
}
//---------------------------------------------------------------------------//
//Update imposed vector field
//~ void immersedBody::updateVectorField(volVectorField& VS, word VName)
void immersedBody::updateVectorField(volVectorField& VS, word VName,volScalarField& body, vectorField surfNorm)
{
    //Check dictionary for parameters (only noSlip allowed)
    word BC = immersedDict_.subDict(VName).lookup("BC");
    
    if (BC=="noSlip")
    {
        //If STATICBODY set to zero
        if ( bodyOperation_==0)
        {
            forAll (surfCells_[Pstream::myProcNo()],cell)
            {
                label cellI = surfCells_[Pstream::myProcNo()][cell];
                VS[cellI]   = Vel_;
            }
            forAll (intCells_[Pstream::myProcNo()],cell)
            {
                label cellI = intCells_[Pstream::myProcNo()][cell];
                VS[cellI] = Vel_;
            }
        }
        else
        {
            label cellI;
            forAll (surfCells_[Pstream::myProcNo()],cell)
            {
                cellI=surfCells_[Pstream::myProcNo()][cell];
                
                // correction for estimated interface position
                scalar intDist(0);
                //~ labelList cellNb(mesh_.cellCells()[cellI]);//list of neighbours
                //~ forAll (cellNb,nbCellI)
                //~ {
                    //~ vector rVec(mesh_.C()[cellNb[nbCellI]] - mesh_.C()[cellI]);
                    //~ intDist += mag(rVec);
                //~ }
                //~ intDist /= (scalar(cellNb.size())+SMALL);
                point surfPoint(mesh_.C()[cellI]);
                
                if (sdBasedLambda_)
                {
                    scalar minMaxBody(max(min(body[cellI],1.0-SMALL),SMALL));
                    // NoteMI: this is necessary for robustness - Foam does not like body=1 in atanh
                    intDist = Foam::atanh(2.0*minMaxBody - 1.0)*Foam::pow(mesh_.V()[cellI],0.333)/intSpan_;
                    surfPoint += surfNorm[cellI]*intDist;
                }
                else
                {
                    labelList cellNb(mesh_.cellCells()[cellI]);//list of neighbours
                    forAll (cellNb,nbCellI)
                    {
                        vector rVec(mesh_.C()[cellNb[nbCellI]] - mesh_.C()[cellI]);
                        intDist += mag(rVec);
                    }
                    intDist /= (scalar(cellNb.size())+SMALL);
                    surfPoint -= intDist*surfNorm[cellI]*(0.5-body[cellI]);
                }
                
                //~ if (body[cellI] < 0.5) {surfPoint -= intDist*surfNorm[cellI]*(0.5-body[cellI]);}
                
                // end of interface position estimate
                
                
                
                //~ vector planarVec       =  mesh_.C()[cellI] - CoM_
                                     //~ - Axis_*(
                                          //~ (mesh_.C()[cellI]-CoM_)&Axis_
                                         //~ );
                vector planarVec       =  surfPoint - CoM_
                                     - Axis_*(
                                          (surfPoint-CoM_)&Axis_
                                         );
                
                vector VSvalue = -(planarVec^Axis_)*omega_ + Vel_;
                VS[cellI] = VSvalue;
            }
            forAll (intCells_[Pstream::myProcNo()],cell)
            {
                cellI=intCells_[Pstream::myProcNo()][cell];
                
                vector planarVec       =  mesh_.C()[cellI] - CoM_
                                     - Axis_*(
                                          (mesh_.C()[cellI]-CoM_)&Axis_
                                         );
                
                vector VSvalue = -(planarVec^Axis_)*omega_ + Vel_;
                VS[cellI] = VSvalue;
            }
            // Note (MI): the "else" part is enough even for static body,
            //            but this way, I am attempting to save a very small
            //            ammount of computing time
            // Note (MI): ^ ... cross product, & ... dot product
        }
    
    }

}
//---------------------------------------------------------------------------//
//Reset body field for this immersed object
void immersedBody::resetBody(volScalarField& body, bool resethistoryFt)
{
    forAll (intCells_[Pstream::myProcNo()],cellI) // modified version
    {
        body[intCells_[Pstream::myProcNo()][cellI]] = 0;
    }
    forAll (surfCells_[Pstream::myProcNo()],cellI)
    {
        body[surfCells_[Pstream::myProcNo()][cellI]] = 0;
    }
    switchWallContact(false);
    
    interpolationInfo_[Pstream::myProcNo()].clear();
    interpolationVecReqs_[Pstream::myProcNo()].clear();
    
    surfCells_[Pstream::myProcNo()].clear();
    intCells_[Pstream::myProcNo()].clear();
    ibPartialVolume_[Pstream::myProcNo()] = 0;
    
    //keep last Ft force only if it was updated in this time step
    if (resethistoryFt)
    {
        DynamicList<Tuple2<label,Tuple2<label,vector>>> historyFtNew;
        
        forAll (historyFt_,Fti)
        {
            if (historyFt_[Fti].second().first())
            {
                Tuple2<label,vector> help(0,historyFt_[Fti].second().second());
                Tuple2<label,Tuple2<label,vector>> help2(historyFt_[Fti].first(), help);
                
                historyFtNew.append(help2);
            }
        }
        
        historyFt_ = historyFtNew;
        contactTimeRes_ = 1;
    }
}
//---------------------------------------------------------------------------//
// function to compute local particle radii w.r.t. CoM_
DynamicList<scalar> immersedBody::getLocPartRad(DynamicLabelList& cellsOfInt)
{
    DynamicScalarList aLst;
    
    forAll (cellsOfInt,cellI)
    {
        label cCell(cellsOfInt[cellI]);
        aLst.append(mag(mesh_.C()[cCell] - CoM_));
    }
    return aLst;
}
//---------------------------------------------------------------------------//
// function to move the body after the contact
void immersedBody::postContactUpdateBodyField(volScalarField& body, volScalarField& refineF)
{
    // move the body due to the contact
    moveImmersedBody();
    if (recomputeProjection_)
    {
        resetBody(body);
        createImmersedBody(body,refineF);
        checkIfInDomain(body);
    }
    updateOldMovementVars();
}
//---------------------------------------------------------------------------//
void immersedBody::recreateBodyField(volScalarField& body, volScalarField& refineF)
{    
    octreeField_ = Field<label>(mesh_.nCells(), 0);
    interpolationInfo_[Pstream::myProcNo()].clear();
    interpolationVecReqs_[Pstream::myProcNo()].clear();
    
    surfCells_[Pstream::myProcNo()].clear();
    intCells_[Pstream::myProcNo()].clear();
    ibPartialVolume_[Pstream::myProcNo()] = 0;
    
    createImmersedBody(body,refineF);
    checkIfInDomain(body);
    Info << "-- body " << bodyId_ << " Re-created" << endl;
}
//---------------------------------------------------------------------------//
// function to compute maximal and mean courant number of the body
void immersedBody::computeBodyCoNumber()
{
    label auxCntr(0);
    scalar VelMag(mag(Vel_));
    
    meanCoNum_ = 0.0;
    
    // rotation body courant number
    scalar rotCoNumB(omega_*dC_*0.5*mesh_.time().deltaT().value());
    
    //Info << "-- rotCo*dCell: " << rotCoNumB << endl;
    
    forAll (surfCells_[Pstream::myProcNo()],sCellI)
    {
        label cellI(surfCells_[Pstream::myProcNo()][sCellI]);
        
        scalar dCell(Foam::pow(mesh_.V()[cellI],0.3333));
        scalar CoNumCell(VelMag*mesh_.time().deltaT().value()/dCell);
        
        CoNumCell+=rotCoNumB/dCell;
        // Note (MI): this formula overshoots the cell courant number
        //            (each cell is assumed to be on the cirvumventing
        //            circle)
        
        CoNum_      = max(CoNum_,CoNumCell);
        meanCoNum_ += CoNumCell;
        auxCntr    += 1;
    }
    forAll (intCells_[Pstream::myProcNo()],iCellI)
    {
        label cellI(intCells_[Pstream::myProcNo()][iCellI]);
        
        scalar dCell(Foam::pow(mesh_.V()[cellI],0.3333));
        scalar CoNumCell(VelMag*mesh_.time().deltaT().value()/dCell);
        
        CoNum_      = max(CoNum_,CoNumCell);
        meanCoNum_ += CoNumCell;
        auxCntr    += 1;
    }
    
    reduce(meanCoNum_, sumOp<scalar>());
    reduce(auxCntr, sumOp<scalar>());
    reduce(CoNum_, maxOp<scalar>());
    
    if(auxCntr > 0)
    {
        meanCoNum_ /= auxCntr;
    }
    
    Info << "-- body " << bodyId_ << " Courant Number mean: " << meanCoNum_
         << " max: " << CoNum_ << endl;
        
    // Note (MI): these two numbers will be different only on
    //            unstructured meshes (or structured meshes with non-uniform
    //            cells)
}

//---------------------------------------------------------------------------//
// print out body linear and angular momentum
void immersedBody::printMomentum()
{
    vector L(I_&(Axis_*omega_));
    vector p(M_*Vel_);
    
    Info << "-- body " << bodyId_ << "  linear momentum:" << p 
         << " magnitude: " << mag(p) <<endl;
    Info << "-- body " << bodyId_ << " angular momentum:" << L 
         << " magnitude: " << mag(L) <<endl;
}
//---------------------------------------------------------------------------//
// print out body statistics
void immersedBody::printStats()
{
    vector L(I_&(Axis_*omega_));
    vector p(M_*Vel_);
    
    Info << "-- body " << bodyId_ << "  linear momentum:" << p 
         << " magnitude: " << mag(p) <<endl;
    Info << "-- body " << bodyId_ << " angular momentum:" << L 
         << " magnitude: " << mag(L) <<endl;
    Info << "-- body " << bodyId_ << "  linear velocity:" << Vel_ 
         << " magnitude: " << mag(Vel_) <<endl;
    Info << "-- body " << bodyId_ << " angular velocity:" << omega_ 
         << " magnitude: " << mag(omega_) <<endl;
    Info << "-- body " << bodyId_ << "    rotation axis:" << Axis_ 
         << " magnitude: " << mag(Axis_) <<endl;
}
//---------------------------------------------------------------------------//
// print out eulerian forces acting on the body
void immersedBody::printForcesAndTorques()
{
    const uniformDimensionedVectorField& g =
        mesh_.lookupObject<uniformDimensionedVectorField>("g");
    vector FG(M_*(1.0-rhoF_.value()/rhoS_.value())*g.value());
    
    Info << "-- body " << bodyId_ << "     linear force:" << F_ 
         << " magnitude: " << mag(F_) <<endl;
    Info << "-- body " << bodyId_ << "   grav/buy force:" << FG 
         << " magnitude: " << mag(FG) <<endl;
    Info << "-- body " << bodyId_ << "    viscous force:" << F_ - FG 
         << " magnitude: " << mag(F_ - FG) <<endl;
    Info << "-- body " << bodyId_ << "           torque:" << T_ 
         << " magnitude: " << mag(T_) <<endl;
}

//---------------------------------------------------------------------------//
// write the bodySurfMesh_ as STL
void immersedBody::writeBodySurfMesh()
{
    if (writeBodySurfMesh_)
    {
        word tmPath(mesh_.time().timePath());
        triSurface triToRet(bodySurfMesh_);
        triToRet.write(tmPath+"/"+name(bodyId_) + ".stl",".stl");
    }
}
//---------------------------------------------------------------------------//
// return to history position when the particle gets in contact during the time step
void immersedBody::returnPosition()
{
    bodySurfMesh_.movePoints(historyBodyPoints_);
    boundBox bound(historyBodyPoints_);
    minBoundPoint_ = bound.min();
    maxBoundPoint_ = bound.max();
}
//---------------------------------------------------------------------------//
// initialize variables base to history. Set history variables only when needed
void immersedBody::initializeVarHistory(bool setHistory)
{
    if (setHistory)
    {
        Axis_ = historyAxis_;
        omega_ = historyOmega_;
        Vel_ = historyVel_;
        a_ = historya_;
        alpha_ = historyAlpha_;
        totalAngle_ = historyTotalAngle_;
    }
    F_*=0.0;
    T_*=0.0;
}
//---------------------------------------------------------------------------//
void immersedBody::solveWallContact
(
    scalar kWN,
    scalar gammaWN,
    scalar kWt,
    scalar gammaWt,
    scalar muW,
    scalar adhWN,
    scalar deltaT
)
{    
    Info << "-- Body " << bodyId_ << " is in contact with wall" << endl;
    
    // compute mean model parameters
    scalar aKN(0.5*(kN_+kWN));
    scalar aGammaN(0.5*(gammaN_+gammaWN));
    scalar aKt(0.5*(kt_+kWt));
    scalar aGammat(0.5*(gammat_+gammaWt));
    scalar amu(0.5*(mu_+muW));
    scalar aadhN(0.5*(adhN_+adhWN));
    
    label nContactFaces(wallContactFaces_[Pstream::myProcNo()].size());
    //Create placeholders for forces
    vector FN(vector::zero);
    vector Ft(vector::zero);
    vector cLVec(vector::zero);
    vector nVecF(vector::zero);
    scalar overallContactArea(0);
    
    //if the IB was in contact in previous DEM time step, find the information about tangential force and assigne it
    vector FtLast(vector::zero);
    bool FtLastFinded(false);
    forAll (historyFt_,Fti)
    {
        if (historyFt_[Fti].first() == -1)
        {
            FtLastFinded = true;
            FtLast = historyFt_[Fti].second().second();
            break;
        }
    }
    
    DynamicLabelList contactCells;
    vector cVel(vector::zero);
    
    // loop over all faces in contact with walls
    forAll (wallContactFaces_[Pstream::myProcNo()],faceI)
    {        
        label cFace(wallContactFaces_[Pstream::myProcNo()][faceI]);
    
        contactCells.append(mesh_.faceOwner()[cFace]);
        // get the local wall velocity
//             label cFacePatch(mesh_.boundaryMesh().whichPatch(cFace));
        //~ vector wVel(U.boundaryField()[cFacePatch][cFace]);
        vector wVel(vector::zero);
        // Note (MI): could I go around whichPatch?
        
        // compute normal to movement and relative velocitydeltaTDEM
        vector nVec(-mesh_.Sf()[cFace]/mag(mesh_.Sf()[cFace]));

        nVecF += nVec * mag(mesh_.Sf()[cFace]);
        overallContactArea += mag(mesh_.Sf()[cFace]);        
        
        //project last Ft to new tangential direction
        vector FtLastP(FtLast - (FtLast & nVec) * nVec);
        //scale the projected vector to remain the magnitude
        vector FtLastr(mag(FtLast) * (FtLastP/(mag(FtLastP)+SMALL)));
        //Evaluate tangential velocity
        vector planarVec       =  mesh_.Cf()[cFace] - CoM_
                                     - Axis_*(
                                          (mesh_.Cf()[cFace]-CoM_)&Axis_
                                         );
        vector cVeli(-(planarVec^Axis_)*omega_ + Vel_);
        vector Vt((cVeli-wVel) - ((cVeli-wVel) & nVec) * nVec);
        cVel += cVeli;
        //Compute tangential force
        Ft += (FtLastr - aKt*Vt*deltaT - aGammat*Vt);
        
        // compute the position vector of the current face and add it
        // to the mean
//         label cFacePatch(mesh_.boundaryMesh().whichPatch(cFace));
//         cLVec += mesh_.boundaryMesh()[cFacePatch].faceCentres()[cFace]-CoM_;mesh_.Cf()[faceInd]
        cLVec += mesh_.Cf()[cFace]-CoM_;
    }

    
    reduce(nContactFaces, sumOp<label>());
    reduce(overallContactArea, sumOp<scalar>());
    reduce(Ft, sumOp<vector>());
    reduce(cVel, sumOp<vector>());
    reduce(cLVec, sumOp<vector>());
    reduce(nVecF, sumOp<vector>());
    
    vector wVel(vector::zero);
    cLVec /= (nContactFaces+SMALL);
    cVel /= (nContactFaces+SMALL);
    nVecF /= (overallContactArea+SMALL);
    scalar VnF(-(cVel-wVel) & nVecF);
    Ft /= (nContactFaces+SMALL);
    Info << "Ft: " << Ft << endl;
    
    scalar intersectedVolume((M0_-M_)/(rhoS_.value() + SMALL));
    scalar Lc(4*mag(cLVec)*mag(cLVec)/(mag(cLVec)+mag(cLVec)));
    scalar reduceM(M0_*M0_/(M0_+M0_));
    
    FN = (aKN*intersectedVolume/(Lc+SMALL) + aGammaN*sqrt(aKN*reduceM/pow(Lc+SMALL,3))*(VnF*overallContactArea))*nVecF;
    
    if (mag(Ft) > amu * mag(FN))
    {
        Ft *= amu * mag(FN) / mag(Ft);
    }
    Info << "Ftamu: " << Ft << endl;
    
    scalar FAc(aadhN*overallContactArea);  
    scalar FAeq(aKN*((adhEqui_*M0_)/(rhoS_.value() + SMALL))/(Lc+SMALL));
    scalar partMul((M0_-M_)/(M0_+SMALL)/adhEqui_);
    if(partMul > 1)
    {
        partMul = 1;
    }
    vector FA((FAeq * partMul  + FAc * (1-partMul)) * nVecF);
    Info << "FN: " << FN << endl;
    Info << "FAeq: " << FAeq << endl;
    Info << "FAc: " << FAc << endl;
    Info << "FA: " << FA << endl;
    FN -= FA;
    
    // Update or add the history of tangential force
    if (FtLastFinded)
    {
        forAll (historyFt_,Fti)
        {
            if (historyFt_[Fti].first() == -1)
            {
                Tuple2<label,vector> help(1,Ft);
                historyFt_[Fti].second() = help;
                break;
            }
        }
    }
    else
    {
        Tuple2<label,vector> help(1,Ft);
        Tuple2<label,Tuple2<label,vector>> help2(-1, help);
                
        historyFt_.append(help2);
    }

    F_+=FN + Ft;
    T_+=-1*cLVec ^ FN; // MS: why 0.1??
}
//---------------------------------------------------------------------------//
// switch the particle off (remove it from the simulation)
void immersedBody::switchActiveOff
(
    volScalarField& body
)
{
    // turn of the particle
    isActive_ = false;
    
    // rewrite the body field
    resetBody(body);
}
//---------------------------------------------------------------------------//
//update movement variables of the body
bool immersedBody::checkContactMovement
(
    scalar deltaT
)
{
    if (mag(deltaT + 1.0) < SMALL) deltaT = mesh_.time().deltaT().value();
    
    vector ai(vector::zero);
    vector Veli(Vel_);
    vector Axisi = Axis_;
    scalar omegai = omega_;
    vector alphai = alpha_;
    
    // auxiliary (nested) functions
    // Note: due to the different body operations, it is better to split
    //       the different parts of the movement update in separate steps
    auto updateTranslation = [&]()
    {
        //compute current acceleration (assume constant over timeStep)
        ai  = F_/(M0_+SMALL);
        
        //Update body linear velocity
        Veli += deltaT*ai;
    };
    
    auto updateRotation = [&]()
    {
        //Update body angular acceleration
        alphai = inv(I_) & T_;
        
        //Update body angular velocity
        vector Omega(Axis_*omegai + deltaT*alphai);
        
        //Split Omega into Axis_ and omega_
        omegai = mag(Omega);
        
        if (omegai < SMALL)
        {
            Axisi = vector::one;
        }
        else
        {
            //~ vector oldAxis = Axis_;
            Axisi =  Omega/(omegai+SMALL);
            forAll (Axisi,axElI)
            {
                if (mag(Axisi[axElI]) < 1.0e-08) Axisi[axElI] = 0.0;
            }
        }
        Axisi /= mag(Axisi);
        // Note (MI): I am cutting of small elements of Axis_ (robustness)
    };
    
    auto updateRotationFixedAxis = [&]()
    {
        //Update body angular velocity
        vector Omega(Axisi*omegai + deltaT * ( inv(I_) & T_ ));
        
        //Split Omega into Axis_ and omega_
        omegai = mag(Omega);
        
        vector newAxis = Omega/(omegai+SMALL);
        if ((newAxis & Axisi) < 0) Axisi *= (-1.0);;
        // Note (MI): can this take care of rotation direction?
    };
    
    if (bodyOperation_ == 0 or bodyOperation_ == 3)
    {
        if (mag(Veli)/(mag(Vel_)+SMALL) > maxDistInDEMloop_ || (omegai)/(omega_+SMALL) > maxDistInDEMloop_) return false;
        return true;
    }
    else if (bodyOperation_ == 1)
    {
        updateRotation();

        if (mag(Veli)/(mag(Vel_)+SMALL) > maxDistInDEMloop_ || (omegai)/(omega_+SMALL) > maxDistInDEMloop_) return false;
        return true;
    }
    else if (bodyOperation_ == 2)
    {
        updateTranslation();
 
        if (mag(Veli)/(mag(Vel_)+SMALL) > maxDistInDEMloop_ || (omegai)/(omega_+SMALL) > maxDistInDEMloop_) return false;
        return true; 
    }
    else if (bodyOperation_ == 4)
    {
        updateRotationFixedAxis();

        if (mag(Veli)/(mag(Vel_)+SMALL) > maxDistInDEMloop_ || (omegai)/(omega_+SMALL) > maxDistInDEMloop_) return false;
        return true; 
    }
         
    updateTranslation();
    if (updateTorque_) updateRotation();
    
    if (mag(Veli)/(mag(Vel_)+SMALL) > maxDistInDEMloop_ || (omegai)/(omega_+SMALL) > maxDistInDEMloop_) return false;
    return true; 
    
    // Note (MI): after the body movement update, I should discard the
    //            forces used for the update
}
//---------------------------------------------------------------------------//
void immersedBody::assignFullHistory()
{
    // Assigned variables for potential contact correction
    historyAxis_ = Axis_;
    historyOmega_ = omega_;
    historyVel_ = Vel_;
    historya_ = a_;
    historyAlpha_ = alpha_;
    historyTotalAngle_ = totalAngle_;
    historyBodyPoints_ = bodySurfMesh_.points();
}
//---------------------------------------------------------------------------//
void immersedBody::initSyncWithFlow(const volVectorField& U)
{
    // estimate particle linear velocity
    //~ vector meanV(vector::zero);scalar totVol(0);
    //~ label  cellI;
    //~ forAll (intCells_[Pstream::myProcNo()],iCellI)
    //~ {
        //~ cellI   = intCells_[Pstream::myProcNo()][iCellI];
        //~ meanV  += U[cellI]*mesh_.V()[cellI];
        //~ totVol += mesh_.V()[cellI];
    //~ }
    //~ reduce(meanV, sumOp<vector>());
    //~ reduce(totVol, sumOp<scalar>());
    //~ Vel_ = meanV/totVol;
    
    // smarter and faster albeir broken implementation
    // -> oF does not want to see region0SubSet/fvSchemes.gradSchemes
    //~ fvMeshSubset intCellsMesh(mesh_);
    //~ intCellsMesh.setCellSubset(intCells_[Pstream::myProcNo()]);
    //~ volVectorField intU("inBodyU",intCellsMesh.interpolate(U));
    //~ volVectorField curlU = fvc::curl(intU);
    
    // auxiliary computation (unnecessarily expensive)
    volVectorField curlU(fvc::curl(U));
    // Note (MI): if this initialization proves OK, than this needs to
    //            be computed only ONCE for all the bodies and re-used
    
    // computation itself
    vector meanV(vector::zero);
    scalar totVol(0);
    vector meanC(vector::zero);
    label  cellI;
    forAll (intCells_[Pstream::myProcNo()],iCellI)
    {
        cellI   = intCells_[Pstream::myProcNo()][iCellI];
        meanV  += U[cellI]*mesh_.V()[cellI];
        meanC  += curlU[cellI]*mesh_.V()[cellI];
        totVol += mesh_.V()[cellI];
    }
    reduce(meanV, sumOp<vector>());
    reduce(meanC, sumOp<vector>());
    reduce(totVol, sumOp<scalar>());
    Vel_ = meanV/(totVol+SMALL);
    meanC/=(totVol+SMALL);
    vector Omega(0.5*meanC);
    if(updateTorque_)
    {
        omega_ = mag(Omega);
        if (omega_ < SMALL)
        {
            Axis_ = vector::one;
            if (mesh_.nGeometricD() < 3)
            {
                const vector validDirs = (mesh_.geometricD() + Vector<label>::one)/2;
                Axis_ -= validDirs;
            }
        }
        else
        {
            Axis_ =  Omega/(omega_+SMALL);
            if (mesh_.nGeometricD() < 3)
            {// in 2D, I need to keep only the correct part of the rotation axis
                const vector validDirs = (mesh_.geometricD() + Vector<label>::one)/2;
                Axis_ = cmptMultiply(vector::one-validDirs,Axis_);
            }
        }
        Axis_ /= mag(Axis_);
    }
    // update old storage
    VelOld_     = Vel_;
    omegaOld_   = omega_;
    AxisOld_    = Axis_;
    // print data:
    Info << "-- body " << bodyId_ << "initial movement variables:" << endl;
    printStats();
}
//---------------------------------------------------------------------------//
void immersedBody::pimpleUpdate
(
    volScalarField& body,
    volVectorField& f
)
{
    updateCoupling(body,f);
    updateMovement(VelOld_,AxisOld_,omegaOld_,velRelaxFac_);
}
//---------------------------------------------------------------------------//
void immersedBody::computeBodyCharPars()
{
    forAll (surfCells_[Pstream::myProcNo()],sCellI)
    {
        label cellI = surfCells_[Pstream::myProcNo()][sCellI];
        dC_ = max(dC_,mag(CoM_-mesh_.C()[cellI]));
    }
    M0_ = M_;
    reduce(dC_, maxOp<scalar>());
}
//---------------------------------------------------------------------------//
void immersedBody::checkIfInDomain(volScalarField& body)
{
    Info << "-- body " << bodyId_ << " current M/M0: " << M_/M0_ << endl;
    //if only 1% of the initial particle mass remains in the domain, switch it off
    if (M_/(M0_+SMALL) < 1e-2) switchActiveOff(body);
}
//---------------------------------------------------------------------------//
void immersedBody::setRestartSim(vector vel, scalar angVel, vector axisRot)
{
    Vel_ = vel;
    omega_ = angVel;
    Axis_ = axisRot;
}
