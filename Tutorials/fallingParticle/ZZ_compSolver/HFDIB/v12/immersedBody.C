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

Contributors
    Federico Municchi (2016)
    Martin Isoz (2019)

Notes on implementation:
 - would it be profitable to rewrite this using cellZone?
   (create cellZone based on STL, work in the cellZone, similarly to the
    porousZone in porousSimpleFoam) 20190201 -> postponed at the moment
 - there seems to be a problem with IB transRotMotion
   (the problem was replicated on both static and dynamic mesh)
\*---------------------------------------------------------------------------*/
#include "immersedBody.H"
//~ #include "fvMesh.H"
//~ #include "dynamicFvMesh.H"
#include "polyMesh.H"
#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
//~ #include <cmath>
//~ #include <algorithm>

#include "interpolationCellPoint.H"
#include "interpolationCell.H"
#include "meshSearch.H"
#include "List.H"

#define ORDER 2

using namespace Foam;

// Q (MI)   : why are there doubled includes?
// Note (MI): removed the doubled includes and STD includes

// Note (MI): I assume that no particles are generated in contact with
//            others or with walls

//---------------------------------------------------------------------------//
immersedBody::immersedBody
(
    word fileName,
    const Foam::dynamicFvMesh& mesh,
    dictionary& HFDIBDict,
    dictionary& transportProperties,
    label bodyId
)
:
isFirstUpdate_(true),
isInPrtContact_(false),
isInWallContact_(false),
immersedDict_(HFDIBDict.subDict(fileName)),
mesh_(mesh),
transportProperties_(transportProperties),
M_(0.0),
CoM_(vector::zero),
Axis_(vector::one),
omega_(0.0),
Vel_(vector::zero),
I_(symmTensor::zero),
F_(vector::zero),
T_(vector::zero),
thrSurf_(readScalar(HFDIBDict.lookup("surfaceThreshold"))),
bodyId_(bodyId),
bodySurfMesh_
(
    IOobject
    (
        fileName +".stl",
        "constant",
        "triSurface",
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
//~ bodyPoints_(bodySurfMesh_.points())
kN_(readScalar(immersedDict_.lookup("kN"))),
gammaN_(readScalar(immersedDict_.lookup("gammaN")))
{
    
    Info<< "Read Immersed Boundary triSurface" << endl;
    bodySurfMesh_.writeStats(Info);
    
    //Return if devlared static
    if(immersedDict_.found("staticBody") )
    {
        bodyOperation_=STATICBODY;
        Info << fileName << " is static body." << endl;
    }
    else if(immersedDict_.found("transRotatingBody"))
    {
         bodyOperation_=TRANSROTATINGBODY;
        
         //Get basic quantities from dict
         Axis_ = immersedDict_.subDict("transRotatingBody").lookup("axis");
         CoM_  = immersedDict_.subDict("transRotatingBody").lookup("center");
        
         omega_  = readScalar(
                        immersedDict_.subDict("transRotatingBody").lookup("omega")
                    );
        
         Vel_   = immersedDict_.subDict("transRotatingBody").lookup("velocity");
        
        
         Info << fileName << " has scripted-trans rotational motion." << endl;
    }
    else if(immersedDict_.found("fluidCoupling"))
    {
        bodyOperation_=FLUIDCOUPLING;
        Info << fileName << " is coupled with fluid phase." << endl;
    }
    else
    {
        Info << "No body operation was found for " << fileName << endl
          << "Assuming static body.";
        bodyOperation_=STATICBODY;
    }
}
//---------------------------------------------------------------------------//
immersedBody::~immersedBody()
{
    bodySurfMesh_.clearOut();
    //~ delete bodySurfMesh_;
}
//---------------------------------------------------------------------------//
//Update immersed body
void immersedBody::updateBodyField
(
    volScalarField& body,
    volVectorField& f
)
{
    if(isFirstUpdate_)
    {
        createImmersedBody( body );
        isFirstUpdate_ = false;
    }
    else
    {
        updateImmersedBody( body, f );
    }
}
//---------------------------------------------------------------------------//
//Create immersed body info
void immersedBody::createImmersedBody(volScalarField& body )
{
        
    const triSurface ibTemp( bodySurfMesh_);
    triSurfaceSearch ibTriSurfSearch( ibTemp );
    const pointField& pp = mesh_.points();

    Info << "clearing old list contents" << endl;
    intCells_.clear();
    surfCells_.clear();
    prtContactCells_.clear();
    wallContactCells_.clear();
    Info << "cleared old list contents" << endl;
    
    // Note: this function may be rewritten using 2 loops but no
    //       dynamically updated fields (MI)
    // Note: I think this is actually faster (better memory management)
    
    Info << "Defining the body field" << endl;
    // first loop, construction of body field and identification of 
    // the number of inside and surface cells
    //~ forAll(mesh_.C(),cellI)
    forAll(body,cellI)
    {
        //Check if partially or completely inside
        const labelList& vertexLabels = mesh_.cellPoints()[cellI];
        const pointField vertexPoints(pp,vertexLabels);
        boolList vertexesInside = ibTriSurfSearch.calcInside( vertexPoints );
        scalar rVInSize(1.0/vertexesInside.size());
        // Note: weight of a single vertex in the cell
        
        scalar nVertIn(0);
        forAll(vertexesInside, verIn)
        {
            if(vertexesInside[verIn]==true)
            {
                body[cellI] += rVInSize; //fraction of cell covered
                // Note: doesn't this work only for 8 vertex polyhedra?
                //       (for 8+ vertex polyhedra, I can have
                //       body[cellI] > 1 and still have less then 50% of
                //       vertices inside the body -> I should find out
                //       how many vertices does the current cell have
                //       and then update the weights accordingly
                //       (should be enough to just get the length of 
                //       vertexesInside)
                
                //~ body[cellI] += rVInSize + bodyId_; //fraction of cell covered
                // Note:
                nVertIn++;
            }
        }
        
        // test if the current body is in contact with other body
        scalar remBody(body[cellI] - rVInSize*nVertIn);
        if (remBody > SMALL and nVertIn > SMALL)
        {
            prtContactCells_.append(cellI);
            if (!isInPrtContact_) isInPrtContact_ = true;
        }
        // Note: if the body[cellI] is higher than it should be after
        //       addition of the representative number of vertices, some
        //       body had to be already present
        // Note: need to take into account only the cells in the current
        //       body
        
        scalar corrBody(body[cellI] - remBody);
        
        //Check if the current cell is in IB or surface
        if(corrBody>1.0-thrSurf_ and nVertIn>0)
        {
            intCells_.append(cellI);
        }
        else if (corrBody>thrSurf_ and nVertIn>0)
        {
            surfCells_.append(cellI);
        }
        // Note (MI&MS): I go through all the body field cells, even
        //               if they are NOT in the CURRENT body!!!
        // -> go through this... and CHECK IT
        // -> go through this... and CHECK IT
    }
    Info << "Defined the body field" << endl;
    
    Info << prtContactCells_ << endl;
        
    // Note (MI): in the internal cells (intCells), I have all the cells
    //            with body function value > 0.9
    // Note (MI): in the surface cells (surfCells), I have the cells
    //            with body function value \in [0.1,0.9)
    // => the internal and surface cells are different sets
    // Note (MI): to cut of all the cells with body function value
    //            \in (0,0.1) seems rather brutal => look into it
    // Note (MI): included scalar for surface cell threshold (thrSurf_),
    //            which might be later moved to HFDIBDict as user input
    
    //refine body as stated in the dictionary
    //~ refineBody(body,&ibTriSurfSearch, & pp);
    
    Info << "Computing interpolation points" << endl;
    calculateInterpolationPoints(body);
    Info << "Computed interpolation points" << endl;
    Info << "Computing geometrical properties" << endl;
    calculateGeometricalProperties(body);
    Info << "Computed geometrical properties" << endl;
}
//---------------------------------------------------------------------------//
//detect contact with other particles
// Note (MI): this function has a lot of structure similar with
//            createImmersedBody
void immersedBody::detectContact(volScalarField& body)
{
    const triSurface ibTemp( bodySurfMesh_);
    triSurfaceSearch ibTriSurfSearch( ibTemp );
    const pointField& pp = mesh_.points();
    
    prtContactCells_.clear();
    wallContactCells_.clear();
    
    forAll(body,cellI)
    {
        //Check if partially or completely inside
        const labelList& vertexLabels = mesh_.cellPoints()[cellI];
        const pointField vertexPoints(pp,vertexLabels);
        boolList vertexesInside = ibTriSurfSearch.calcInside( vertexPoints );
        scalar rVInSize(1.0/vertexesInside.size());
        // Note: weight of a single vertex in the cell
        
        scalar testBody(0);
        scalar nVertIn(0);
        
        forAll(vertexesInside, verIn)
        {
            if(vertexesInside[verIn]==true)
            {
                testBody += rVInSize; //fraction of cell covered
                nVertIn++;
            }
        }
        
        // test if the current body is in contact with other body
        scalar remBody(body[cellI] - testBody);
        if (remBody > SMALL and nVertIn > SMALL)
        {
            prtContactCells_.append(cellI);
            if (!isInPrtContact_) isInPrtContact_ = true;
        }
        // Note: if the body[cellI] is higher than it should be after
        //       addition of the representative number of vertices, some
        //       body had to be already present
        // Note: need to take into account only the cells in the current
        //       body
    }
    Info << prtContactCells_ << endl;
}
//---------------------------------------------------------------------------//
//Update immersed body info
void immersedBody::updateImmersedBody
(
    volScalarField& body,
    volVectorField& f
)
{
    Info << "WORKING ON BODY " << bodyId_ << endl;
    //~ volScalarField bodyCopy(body);
    
    isInPrtContact_ = false;
    isInWallContact_= false;
    F_*=0.0;
    T_*=0.0;
    
    //Check Operation to perform
    if(bodyOperation_==STATICBODY) return;
    else
    {        
        if(bodyOperation_==FLUIDCOUPLING)
        {
            Info << "Updating coupling" << endl;
            updateCoupling(body,f);
            Info << "Updated coupling" << endl;
        }
        
        detectContact(body);
        if (isInPrtContact_) return;
        
        Info << "Reseting body" << endl;
        resetBody(body);
        Info << "Reseted body" << endl;
        
        Info << "Moving body" << endl;
        moveImmersedBody();
        Info << "Moved body" << endl;
    }
    
    //TODO:Very inefficient find other algorithm
    Info << "Re-creating body" << endl;
    createImmersedBody(body);
    Info << "Re-created body" << endl;
    // Note (MI): is this so bad tho? the bad thing is that I need to
    //            always recreate IB (and I have discrepancies in the
    //            mesh refinement). But how slow is actually the create
    //            function?
    
    //~ if (isInPrtContact_)
    //~ {
        //~ bodySurfMesh_.movePoints(bodyPoints_);
        //~ body          = bodyCopy;
    //~ }
    //~ else
    //~ {
        //~ bodyPoints_ = bodySurfMesh_.points();
    //~ }
    // Note (MI): if the body is in contact, I will just throw the whole
    //            computed movement down the toilet and reverse it
    // Note (MI): if the body is in contact, I will first add the contact
    //            induced forces and move it afterwards with these
    //            included

}
//---------------------------------------------------------------------------//
void immersedBody::updateCoupling
(
    volScalarField& body,
    volVectorField& f
)
{
  //~ const uniformDimensionedVectorField g =
   //~ mesh_.lookupObject<uniformDimensionedVectorField>("g");
   
  const dimensionedScalar rhof =
   dimensionedScalar(
          transportProperties_.lookup("rho")
      );
      
  //~ dimensionedScalar rho_ =
       //~ dimensionedScalar(immersedDict_.lookup("rho"));

  //~ vector F(vector::zero);
  //~ vector T(vector::zero);
  
  Info << "Old viscous force: "<< F_ << endl;
  
  //Calcualate viscous force and torque
  forAll(surfCells_,cell)
  {
     label cellI = surfCells_[cell];

     //~ F_ +=  f[cellI] * mesh_.V()[cellI] * rhof.value();
     F_ -=  f[cellI] * mesh_.V()[cellI] * rhof.value();
     T_ +=  ( (mesh_.C()[cellI]-CoM_) ^ f[cellI] )
                 *mesh_.V()[cellI]* rhof.value();
     // Note (MI): shouldn't this be multiplied by (1-body)?
     //~ F *=  (scalar(1)-body[cellI]);
     //~ T *=  (scalar(1)-body[cellI]);
     // Note (MI): thats what I am trying to do here
  }
  // Note (MI): this loop will be probably crucial for FSI
  
  Info << "Viscous force: " << F_ << endl;

  updateMovement();

  //~ //Update body linear velocity
  //~ Vel_ += mesh_.time().deltaT().value()
            //~ * (
                  //~ F_ / (M_+SMALL)
                //~ + (1.0-rhof.value()/rho_.value())*g.value()
              //~ );
  //~ // Note (MI): added +SMALL for robustness, effect on result not tested
//~ 
  //~ //Update body angular velocity
  //~ vector Omega_(vector::zero);
//~ 
  //~ Omega_ =   Axis_*omega_
           //~ + mesh_.time().deltaT().value() * ( inv(I_) & T_ );
//~ 
  //~ //Split Omega_ into Axis_ and omega_
  //~ omega_ = mag(Omega_);
//~ 
  //~ if(omega_ < SMALL)
  //~ {
     //~ Axis_ = vector::one;
  //~ }
  //~ else
  //~ {
     //~ vector oldAxis = Axis_;
     //~ Axis_ =  Omega_/(omega_+SMALL);
     //~ if ((oldAxis & Axis_) < 0) Axis_ *= (-1.0);
     //~ forAll(Axis_,axElI)
     //~ {
         //~ if (mag(Axis_[axElI]) < 1.0e-08) Axis_[axElI] = 0.0;
     //~ }
     //~ Axis_ /= mag(Axis_);
  //~ }
  //~ // Note (MI): added +SMALL for robustness, effect on result not tested
  //~ // Note (MI): ensured that axis will always have the same orientation
  //~ // Note (MI): I am cutting of small elements of Axis_ (robustness)

}
//---------------------------------------------------------------------------//
//update movement variables of the body
void immersedBody::updateMovement()
{
    const uniformDimensionedVectorField g =
        mesh_.lookupObject<uniformDimensionedVectorField>("g");
    
    const dimensionedScalar rhof =
        dimensionedScalar
        (
            transportProperties_.lookup("rho")
        );
      
    dimensionedScalar rho_ =
        dimensionedScalar(immersedDict_.lookup("rho"));
        
    //Update body linear velocity
    Vel_ += mesh_.time().deltaT().value()
            * (
                  //~ F / M_
                  F_ / (M_+SMALL)
                + (1.0-rhof.value()/rho_.value())*g.value()
              );
    // Note (MI): added +SMALL for robustness, effect on result not tested
    
    //Update body angular velocity
    vector Omega_(vector::zero);
    
    Omega_ =   Axis_*omega_
           + mesh_.time().deltaT().value() * ( inv(I_) & T_ );
    
    //Split Omega_ into Axis_ and omega_
    omega_ = mag(Omega_);
    
    if(omega_ < SMALL)
    {
        Axis_ = vector::one;
    }
    else
    {
        vector oldAxis = Axis_;
        Axis_ =  Omega_/(omega_+SMALL);
        if ((oldAxis & Axis_) < 0) Axis_ *= (-1.0);
        forAll(Axis_,axElI)
        {
            if (mag(Axis_[axElI]) < 1.0e-08) Axis_[axElI] = 0.0;
        }
        Axis_ /= mag(Axis_);
    }
    // Note (MI): added +SMALL for robustness, effect on result not tested
    // Note (MI): ensured that axis will always have the same orientation
    // Note (MI): I am cutting of small elements of Axis_ (robustness)
}

//---------------------------------------------------------------------------//
//Create interpolation points
void immersedBody::calculateInterpolationPoints
(
    volScalarField& body
)
{
    scalar sqrtThree_ = sqrt(3.0);
    meshSearch search_(mesh_);
    
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
    scalarField intDists = sqrtThree_*Foam::pow(mesh_.V(),1.0/3.0);
    // Note (MI): doesn't this work only for hexahedral cells? and even
    //            on such not exactly?
    
    // scale the surfNorm (prepare it for the future computations)
    surfNorm = intDists*surfNorm;
    
    // clear the old interpolation data
    interpolationPoints_.clear();
    interpolationCells_.clear();
    // Note (MI): OK, in this cycle, I go through all the surface cells
    //            in each cell, I
    //            1. approximate a position of surface point
    //            2. 
    forAll(surfCells_,cell)
    {        
        //get surface cell label
        label scell = surfCells_[cell];
        
        //create vector for points and cells and add to main vectors
        DynamicPointList intPoints;
        DynamicLabelList intCells;
        
        // Note (MI): these things get extended for each surfCell
        // Note (MI): what happens during mesh update? does the number of
        //            interface cells change? what should I do then?
        // -> destroy old lists and create another ones? (probably)
        // -> I might need a method to call after each mesh update 
        // Note (MI): I should set the size of an element of the global variable
        //            before        
        
        //Approximate distance using body - (MI): what is this comment?
        point surfPoint = mesh_.C()[scell] + surfNorm[scell]*(0.5-body[scell]);
        
        //Add to list
        intPoints.append(surfPoint);
        
        //Add other interpolation points
        for(int order=0;order<ORDER;order++)
        {
            surfPoint = surfPoint + surfNorm[scell];
            intPoints.append(surfPoint);
        }
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
        for(int order=0;order<ORDER;order++)
        {
            //Use findNearestCell()...it is faster
            //Check if outside the domain (if not, then set to -1)
            label cellI = search_.findCell(intPoints[order+1],scell,false);
            intCells.append(cellI);
        }
        // Note (MI): the code should be able to handle outside of domain pts
                
        // assign to global variables
        interpolationPoints_.append(intPoints);
        interpolationCells_.append(intCells);
    
    }

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
    
    forAll(labelCellLst,cell)
    {
        //~ label cellI = intCells_[cell];
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
    tmpCom = addToMAndI(body,surfCells_,tmpCom);
    tmpCom = addToMAndI(body,intCells_,tmpCom);
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
    
    //~ CoM_ = tmpCom / M_;
    CoM_ = tmpCom / (M_+SMALL);
    // Note (MI): added +SMALL for robustness, effect on result not tested
}
//---------------------------------------------------------------------------//
//Move immersed body according to body operation
void immersedBody::moveImmersedBody()
{


    //Rotation angle
    scalar angle = omega_*mesh_.time().deltaT().value();
    vector transIncr = Vel_*mesh_.time().deltaT().value();
    
    //~ pointField bodyPoints (bodySurfMesh_->points());
    pointField bodyPoints (bodySurfMesh_.points());
    
    //~ Info << mag(Axis_) << endl;
    
    //~ tensor rotMatrix(Foam::cos(angle)*tensor::one);
    //~ rotMatrix += Foam::sin(angle)*tensor(0,-Axis_.z(),Axis_.y(),Axis_.z(),0,-Axis_.x(),-Axis_.y(),Axis_.x(),0);
    //~ rotMatrix += (1.0-Foam::cos(angle))*(Axis_ * Axis_);
    
    Info << Axis_ << endl;
    
    //Move points
    forAll(bodyPoints,p)
    {
        vector normV    = (bodyPoints[p]-CoM_)^Axis_;
        scalar axisMod  = ((bodyPoints[p]-CoM_)&Axis_);
        scalar magDist  = mag((bodyPoints[p]-CoM_) - Axis_*axisMod/(mag(Axis_)+SMALL));
        
        //Move in tangential direction
        //~ bodyPoints[p] = bodyPoints[p] + angle*normV;
        
        //Rotate the point
        //~ Info << bodyPoints[p] << endl;
        //~ bodyPoints[p] = rotMatrix & bodyPoints[p];
        bodyPoints[p] = Axis_*(Axis_ & bodyPoints[p]) + Foam::cos(angle)*((Axis_ ^ bodyPoints[p]) ^ Axis_) + Foam::sin(angle)*(Axis_ ^ bodyPoints[p]);
        //~ Info << bodyPoints[p] << endl;
        //~ Info << "----" << endl;
        
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
    
    //move mesh
    bodySurfMesh_.movePoints(bodyPoints);
}
//---------------------------------------------------------------------------//
//Update imposed vector field
//~ void immersedBody::updateVectorField(volVectorField& VS, word VName)
void immersedBody::updateVectorField(volVectorField& VS, word VName,volScalarField& body)
{
  //Check dictionary for parameters (only noSlip allowed)
  word BC = immersedDict_.subDict(VName).lookup("BC");

  if(BC=="noSlip")
  {
   //If STATICBODY set to zero
   if( bodyOperation_==STATICBODY)
   {
        forAll(surfCells_,cell)
        {
            //~ VS[surfCells_[cell]] = Vel_;
            label cellI = surfCells_[cell];
            VS[cellI]   = body[cellI]*Vel_+(scalar(1.0)-body[cellI])*VS[cellI];
        }
        forAll(intCells_,cell)
        {
            label cellI = intCells_[cell];
            VS[cellI] = Vel_;
        }
        // Note (MI): this should be faster than the single loop with if switch
        // Note (MI): if I do noSlip BC in surfCells, should't I make
        //            a convex combination of current velField and 
        //            prescribed velField?
   }
   else
   {
        label cellI;
        //Apply
        forAll(surfCells_,cell)
        {
            cellI=surfCells_[cell];
            
            vector planarVec       =  mesh_.C()[cellI] - CoM_
                                 - Axis_*(
                                      (mesh_.C()[cellI]-CoM_)&Axis_
                                     );
            
            vector VSvalue = (planarVec^Axis_)*omega_ + Vel_;
            VS[cellI] = VSvalue;
        }
        forAll(intCells_,cell)
        {
            cellI=intCells_[cell];
            
            vector planarVec       =  mesh_.C()[cellI] - CoM_
                                 - Axis_*(
                                      (mesh_.C()[cellI]-CoM_)&Axis_
                                     );
            
            vector VSvalue = (planarVec^Axis_)*omega_ + Vel_;
            VS[cellI] = VSvalue;
        }
        // Note (MI): this is the same stupid eror as before, I loop over
        //            all the internal and surface cells (I broke this)
        // =>         code clean up and computation on both internal and
        //            surface cells - compare with v8 of the code
        // Note (MI): I could move the computation into a separate function
        //            but as it is just a few lines...
        // Note (MI): ^ ... cross product, & ... dot product
   }

  }

}
//---------------------------------------------------------------------------//
//Reset body field for this immersed object
void immersedBody::resetBody(volScalarField& body)
{
    forAll(intCells_,cellI) // modified version
    {
        body[intCells_[cellI]] = 0;
    }
    forAll(surfCells_,cellI)
    {
        body[surfCells_[cellI]] = 0;
    }
    
    interpolationPoints_.clear();
    interpolationCells_.clear();
    // Note (MI): I think this is cleaner (no loop and takes into account
    //            the boundary patches as well)
    // Note (MI): this works well for a single body but it will not work
    //            for multiple bodies (this resets all the bodies in the
    //            domain!!!) (20190206)
    // Note (MI): the new loops are faster than the original loops BUT it
    //            does not take into account boundary contitions
    
    //~ isInPrtContact_ = false;
    //~ isInWallContact_= false;
    prtContactCells_.clear();
    wallContactCells_.clear();
}
//---------------------------------------------------------------------------//
//Refine body field for this immersed object using MC-like algorithm
//Cells are assumed to be hexahedral at the particle surface
//(but can have different edge length)
// Note (MI):
// - the refinement seems to be specific to hexahedral meshes
// - would it be ok to just remove this for the sake of making the code
//   usable on arbitrary meshes? this way, I would be able to refine
//   the mesh locally in the vicinity of the IB boundary and to have
//   a nicer shape and refined body
// Note (MI):
// - this approach to body refinement is usable on hexahedral meshes
//   (small and topologically simple domains)
// - as the fields are not well representable in paraview, the program
//   resuls are bound to be "ugly" (big cells and so on) even though the
//   results are "smooth" and accurate
//~ void immersedBody::refineBody(volScalarField& body,
                              //~ triSurfaceSearch * ibTriSurfSearch,
                              //~ const pointField * pp
                           //~ )
//~ {
    //~ if(!immersedDict_.found("refineMC")) return;
//~ 
//~ 
    //~ scalar nPointsEdge = readScalar(immersedDict_.lookup("refineMC"));
//~ 
    //~ //loop over all the surface cells
    //~ forAll(surfCells_,cell)
    //~ {
//~ 
     //~ label cellI = surfCells_[cell];
//~ 
     //~ scalar deltaV = 1.0/(nPointsEdge*nPointsEdge*nPointsEdge);
//~ 
     //~ //Get cell center
     //~ point centerC = mesh_.C()[cellI];
//~ 
     //~ //Get one node
     //~ //Check if partially or completely inside
     //~ const labelList& vertexLabels = mesh_.cellPoints()[cellI];
     //~ const pointField vertexPoints(*pp,vertexLabels);
     //~ point baseNode = vertexPoints[0];
//~ 
     //~ //create vector representing 3d diagonal of the cell
     //~ vector edgesC = 2*(centerC - baseNode);
//~ 
     //~ //create list of points
     //~ List<point> pointsMC;
//~ 
     //~ //create deltas
     //~ scalar delta_i = edgesC[0]/nPointsEdge;
     //~ scalar delta_j = edgesC[1]/nPointsEdge;
     //~ scalar delta_k = edgesC[2]/nPointsEdge;
//~ 
     //~ //add points to list
     //~ for(int i=0;i<nPointsEdge;i++)
     //~ {
      //~ //point i-coordinate
      //~ scalar icoord = baseNode[0] + delta_i*(i+0.5);
//~ 
      //~ for(int j=0;j<nPointsEdge;j++)
      //~ {
       //~ //point j-coordinate
       //~ scalar jcoord = baseNode[1] + delta_j*(j+0.5);
//~ 
       //~ for(int k=0;k<nPointsEdge;k++)
       //~ {
        //~ //point k-coordinate
        //~ scalar kcoord = baseNode[2] + delta_k*(k+0.5);
//~ 
        //~ //create point
        //~ point p(icoord,jcoord,kcoord);
//~ 
        //~ //add to list
        //~ pointsMC.append(p);
//~ 
       //~ }
      //~ }
     //~ }
//~ 
     //~ //Check who is inside
     //~ pointField pField(pointsMC);
     //~ boolList pInside = ibTriSurfSearch->calcInside( pField );
//~ 
     //~ //Calculate new body
     //~ scalar newbody = 0.0;
     //~ forAll(pInside,p)
     //~ {
      //~ if(pInside[p])
       //~ newbody+=deltaV;
//~ 
//~ 
     //~ }
//~ 
//~ 
     //~ body[cellI] = newbody;
//~ 
    //~ }
//~ 
//~ 
//~ }
//---------------------------------------------------------------------------//
void immersedBody::resetIB(volScalarField& body)
{
    resetBody(body);
    isFirstUpdate_ = true;
    // Note (MI): I reset isFirstUpdate_ variable -> on the next call to
    //            updateBodyField the immersedBody will be re-created
}
//---------------------------------------------------------------------------//
// function to compute local particle radii w.r.t. CoM_
DynamicList<scalar> immersedBody::getLocPartRad(DynamicLabelList& cellsOfInt)
{
    DynamicScalarList aLst;
    
    forAll(cellsOfInt,cellI)
    {
        label cCell(cellsOfInt[cellI]);
        aLst.append(mag(mesh_.C()[cCell] - CoM_));
    }
    return aLst;
}
//---------------------------------------------------------------------------//
// function to move the body after the contact
void immersedBody::postContactBodyUpdate(volScalarField& body)
{
    //~ Info << "Made it to postContactBodyUpdate" << endl;
    
    // include the effects of contact
    Info << "pre-contact  linear  velocity: " <<Vel_ << endl;
    Info << "pre-contact  angular velocity: " <<omega_ << endl;
    updateMovement();
    Info << "post-contact linear  velocity: " <<Vel_ << endl;
    Info << "post-contact angular velocity: " <<omega_ << endl;
    
    // move the body due to the contact
    resetBody(body);
    moveImmersedBody();
    createImmersedBody(body);
    
    // reset the contact flags
    isInPrtContact_ = false;
    isInWallContact_= false;
}

