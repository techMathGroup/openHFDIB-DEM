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
    Martin Isoz (2019-*), Martin Šourek (2019-*), 
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "addModelRepeatRandomPosition.H"
#include "meshSearch.H"

using namespace Foam;

//---------------------------------------------------------------------------//
addModelRepeatRandomPosition::addModelRepeatRandomPosition
(
    const dictionary& addModelDict,
    const word        stlName,
    const Foam::dynamicFvMesh& mesh
)
:
addModelDict_(addModelDict),
addMode_(word(addModelDict_.lookup("addModel"))),
stlName_(stlName),
bodyAdded_(false),
mesh_(mesh),

coeffsDict_(addModelDict_.subDict(addMode_+"Coeffs")),

addDomain_(word(coeffsDict_.lookup("addDomain"))),
scalingMode_(word(coeffsDict_.lookup("scalingMode"))),
rotationMode_(word(coeffsDict_.lookup("rotationMode"))),
addModeI_(word(coeffsDict_.lookup("addMode"))),

addDomainCoeffs_(coeffsDict_.subDict(addDomain_ + "Coeffs")),
scalingModeCoeffs_(coeffsDict_.subDict(scalingMode_ + "Coeffs")),
rotationModeCoeffs_(coeffsDict_.subDict(rotationMode_ + "Coeffs")),
addModeICoeffs_(coeffsDict_.subDict(addModeI_ + "Coeffs")),

useNTimes_(0),
timeBetweenUsage_(0),
partPerAdd_(0),
fieldValue_(0),
addedOnTimeLevel_(0),
partPerAddTemp_(0),

zoneName_(),
minBound_(vector::zero),
maxBound_(vector::zero),

scaleParticles_(false),
minScale_(0),
maxScale_(0),
minScaleFit_(0),
scaleStep_(0),
nTriesBeforeScaling_(0),

rotateParticles_(false),
randomAxis_(false),
axisOfRot_(vector::zero),

bodyAdditionAttemptCounter_(0),
scaleCorrectionCounter_(0),

scaleApplication_(false),
scaleRandomApplication_(false),
rescaleRequirement_(false),
succesfulladition_(false),
scalingFactor_(0),
restartPartCountTemp_(false),
reapeatedAddition_(false),
firstTimeRunning_(true),
cellZoneActive_(false),
boundBoxActive_(false),
octreeField_(mesh_.nCells(), 0),
timeBased_(false),
fieldBased_(false),
fieldCurrentValue_(0),
allActiveCellsInMesh_(true),
randGen_(clock::getTime())
{
	init();
}
    
addModelRepeatRandomPosition::~addModelRepeatRandomPosition()
{
}

//---------------------------------------------------------------------------//
// Note (MI): Initialization of this model is so complex that I moved it to a
//            specific function
void addModelRepeatRandomPosition::init()
{
    // Set sizes to necessary datatypes
    cellsInBoundBox_.setSize(Pstream::nProcs());
    cellZonePoints_.setSize(Pstream::nProcs());
    

	if (addModeI_ == "timeBased")
	{
        useNTimes_ = (readLabel(addModeICoeffs_.lookup("useNTimes")));
		timeBetweenUsage_ = (readScalar(addModeICoeffs_.lookup("timeBetweenUsage")));
		partPerAdd_ = (readLabel(addModeICoeffs_.lookup("partPerAdd")));
        timeBased_ = true;
        Info << "-- addModelMessage-- " << "addModel will control simulation time" << endl;
        Info << "-- addModelMessage-- " << "STL will be re-used " << useNTimes_ << " times" << endl;
        Info << "-- addModelMessage-- " << "STL will be added each " << timeBetweenUsage_ << " [T]" << endl;
        Info << "-- addModelMessage-- " << "upon each addition, " << partPerAdd_ << " bodies will be generated from the given STL" << endl;
	}
	else if (addModeI_ == "fieldBased")
	{
		fieldValue_ = (readScalar(addModeICoeffs_.lookup("fieldValue")));
        fieldBased_ = true;
        Info << "-- addModelMessage-- " << "addModel will control particles volume fraction" << endl;
		Info << "-- addModelMessage-- " << "preset volume fraction: " << fieldValue_ << endl;
	}
    else
    {
        Info << "-- addModelMessage-- " << "notImplemented, will crash" << endl;
    }
	
	if (addDomain_ == "cellZone")
	{
		zoneName_ = (word(addDomainCoeffs_.lookup("zoneName")));
		cellZoneActive_ = true;
        initializeCellZone();
        Info << "-- addModelMessage-- " << "cellZone based addition zone" << endl;
	}
	else if (addDomain_ == "boundBox")
	{
		minBound_       = (addDomainCoeffs_.lookup("minBound"));
		maxBound_       = (addDomainCoeffs_.lookup("maxBound"));
		boundBoxActive_ = true;
        initializeBoundBox();
        Info << "-- addModelMessage-- " << "boundBox based addition zone" << endl;
	}
	else if (addDomain_ == "domain")
	{
		Info << "-- addModelMessage-- " << "notImplemented, will crash" << endl;
	}
    else
    {
		Info << "-- addModelMessage-- " << "notImplemented, will crash" << endl;
	}
    
    // check, if the whole zone is in the mesh
    scalarList procZoneVols(Pstream::nProcs());
    procZoneVols[Pstream::myProcNo()] = 0;
    forAll (cellsInBoundBox_[Pstream::myProcNo()],cellI)
    {
        procZoneVols[Pstream::myProcNo()]+=mesh_.V()[cellsInBoundBox_[Pstream::myProcNo()][cellI]];
    }
    scalar zoneVol(gSum(procZoneVols));
    scalar zoneBBoxVol(cellZoneBounds_.volume());
    if (zoneVol - zoneBBoxVol > 1e-5*zoneBBoxVol)
    {
        allActiveCellsInMesh_ = false;
        Info << "-- addModelMessage-- " 
             << "addition zone NOT completely immersed in mesh "
             << "this computation will be EXPENSIVE" << endl;
        Info << zoneVol << " " << zoneBBoxVol << endl;
    }
    else
    {
        Info << "-- addModelMessage-- " 
             << "addition zone completely immersed in mesh -> OK" << endl;
    }
    // Note (MI): the coding should be done in such a way that all the
    //            variables should be present irrespective of addDomain_
    //            (check initializeBoundBox and initializeCellZone)
	
	if (scalingMode_ == "noScaling")
	{
		scaleParticles_ = false;
        Info << "-- addModelMessage-- " << "all particles will have the same scale" << endl;
	}
	else if (scalingMode_ == "randomScaling")
	{
		minScale_               = (readScalar(scalingModeCoeffs_.lookup("minScale")));
		maxScale_               = (readScalar(scalingModeCoeffs_.lookup("maxScale")));
		scaleParticles_         = false;
		scaleRandomApplication_ = true;
        Info << "-- addModelMessage-- " << "particles will be randomly scaled" << endl;
	}
	else if (scalingMode_ == "scaleToFit")
	{
		minScaleFit_        = (readScalar(scalingModeCoeffs_.lookup("minScale")));
		scaleStep_          = (readScalar(scalingModeCoeffs_.lookup("scaleStep")));
		nTriesBeforeScaling_= (readScalar(scalingModeCoeffs_.lookup("nTriesBeforeScaling")));
        Info << "-- addModelMessage-- " << "particles will be downscaled to better fill the domain" << endl;
		Info << "-- addModelMessage-- " << "nTriesBeforeDownScaling: " << nTriesBeforeScaling_ << endl;
		scaleParticles_ = true;
	}
    else
    {
		Info << "-- addModelMessage-- " << "notImplemented, will crash" << endl;
	}
	
	if (rotationMode_ == "noRotation")
	{
		rotateParticles_ = false;
		randomAxis_      = false;
        Info << "-- addModelMessage-- " << "source STL will not be rotated upon addition" << endl;
	}
	else if (rotationMode_ == "randomRotation")
	{
		rotateParticles_ = true;
		randomAxis_      = true;
        Info << "-- addModelMessage-- " << "source STL will be randomly rotated upon addition" << endl;
	}
	else if (rotationMode_ == "fixedAxisRandomRotation")
	{
		axisOfRot_       = (rotationModeCoeffs_.lookup("axis"));
        Info << "-- addModelMessage-- " << "source STL will be rotated by a random angle around a fixed axis upon addition" << endl;
		Info << "-- addModelMessage-- " << "set rotation axis: " << axisOfRot_ << endl;
		rotateParticles_ = true;
		randomAxis_      = false;
	}
    else
    {
		Info << "-- addModelMessage-- " << "notImplemented, will crash" << endl;
	}
	
	partPerAddTemp_ = partPerAdd_;
    
}

//---------------------------------------------------------------------------//
bool addModelRepeatRandomPosition::shouldAddBody(const volScalarField& body)
{
	
    if (timeBased_)
    {
        scalar timeVal(mesh_.time().value());
        scalar deltaTime(mesh_.time().deltaT().value());
        scalar tmFrac(timeVal/timeBetweenUsage_);
        tmFrac -=  floor(tmFrac+deltaTime);
        
        Info << "-- addModelMessage-- " << "Time/(Time beween usage) - floor(Time/Time beween usage): " 
             << tmFrac << endl;
             
        Info << "-- addModelMessage-- " << "Number of bodies added on this time level: " << addedOnTimeLevel_ << endl;
             
        bool tmLevelOk(tmFrac < deltaTime);
        
        if (not tmLevelOk)
        {
            addedOnTimeLevel_ = 0;
            return false;
        }
        
        if (partPerAdd_ <= addedOnTimeLevel_) {return false;}
                                 
        return (tmLevelOk and useNTimes_ > 0);
    }
    
    if (fieldBased_)
    {
        scalar currentLambdaFrac(checkLambdaFraction(body));
        if (currentLambdaFrac < fieldValue_ )
        {
            Info << "-- addModelMessage-- " << "Current lambda fraction = " << currentLambdaFrac << " < then preset lambda fraction = " << fieldValue_ << endl;
            return true;
        }
    }
    
    return false;
    
}
//---------------------------------------------------------------------------//
triSurface addModelRepeatRandomPosition::addBody
(
    const   volScalarField& body
)
{
    bodyAdditionAttemptCounter_++;
    // load the STL file
    triSurfaceMesh bodySurfMesh
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
    );
    
    
    
    // get the working bounding box center
    point bBoxCenter = cellZoneBounds_.midpoint();
    pointField bodyPoints(bodySurfMesh.points());
    
    // get its center of mass
    vector CoM = gSum(bodyPoints/bodyPoints.size());
    Info << "-- addModelMessage-- " << "original STL-based CoM: " << CoM << endl;
    
    // move the body to the center of the bounding box
    bodyPoints += bBoxCenter-CoM;
    CoM         = bBoxCenter;
    bodySurfMesh.movePoints(bodyPoints);
    
    // rotate
    if (rotateParticles_)
    {
        scalar rotAngle = returnRandomAngle();
        if (randomAxis_)
        {
            axisOfRot_ = returnRandomRotationAxis();
        }
        Info << "-- addModelMessage-- " << "Will rotate by " << rotAngle << " PiRad around axis " << axisOfRot_ << endl;
        
        tensor rotMatrix(Foam::cos(rotAngle)*tensor::I);
		
		rotMatrix += Foam::sin(rotAngle)*tensor(
             0.0,      -axisOfRot_.z(),  axisOfRot_.y(),
             axisOfRot_.z(), 0.0,       -axisOfRot_.x(),
            -axisOfRot_.y(), axisOfRot_.x(),  0.0
		);
        
        rotMatrix += (1.0-Foam::cos(rotAngle))*(axisOfRot_ * axisOfRot_);
        
        bodyPoints -= CoM;
        bodyPoints = rotMatrix & bodyPoints;
        bodyPoints += CoM;
        
        bodySurfMesh.movePoints(bodyPoints);
        //~ CoM = gSum(bodySurfMesh.coordinates())/bodySurfMesh.size();
    }
    
    // scale
    if (scaleApplication_ or scaleRandomApplication_)
    {
        if (scaleRandomApplication_){scaleStep_ = returnRandomScale();}
        bodySurfMesh.scalePoints(scaleStep_);
        //~ CoM = gSum(bodySurfMesh.coordinates())/bodySurfMesh.size();
    }
    // Note (MI): there should be no change in CoM after rotation and
    //            scaling BUT CoM is approximate...
    
    // translate
    bodyPoints = bodySurfMesh.points();
    if (mesh_.nGeometricD() < 3)
    {
        const vector validDirs = (mesh_.geometricD() + Vector<label>::one)/2;
        CoM -= cmptMultiply((vector::one - validDirs),CoM);
        CoM += cmptMultiply((vector::one - validDirs),0.5*(mesh_.bounds().max() + mesh_.bounds().min()));
    }//project CoM onto current solution plane (if needed)
    
    vector randomTrans = returnRandomPosition(CoM,bodySurfMesh);
    bodyPoints += randomTrans;
    
    bodySurfMesh.movePoints(bodyPoints);
    
    // check if the body can be added
    bool canAddBodyI(canAddBody(body,bodySurfMesh));
    reduce(canAddBodyI, andOp<bool>());    
    bodyAdded_ = (canAddBodyI);
    
    if(!bodyAdded_)
	{
		scaleCorrectionCounter_++; 
	}
	
	if(bodyAdded_)
	{			
		if(timeBased_)
		{
			Info << "-- addModelMessage-- " << "addedOnTimeLevel:  " << addedOnTimeLevel_<< endl;
			addedOnTimeLevel_++;
			Info << "-- addModelMessage-- " << "bodyAdded: " << bodyAdded_ << " addedOnTimeLevel:  " << addedOnTimeLevel_<<" useNTimes: " << useNTimes_<<  endl;
			if(addedOnTimeLevel_ == partPerAdd_)
			{
				useNTimes_--;
				Info << "-- addModelMessage-- " <<" useNTimes: " << useNTimes_<<  endl;
				reapeatedAddition_ = false;
			}
		}
		
		scaleCorrectionCounter_ = 0;
		
	}
	
	Info << "-- addModelMessage-- " << "bodyAdditionAttemptNr  : " << bodyAdditionAttemptCounter_<< endl;
	Info << "-- addModelMessage-- " << "sameScaleAttempts      : " << scaleCorrectionCounter_<< endl;
	
	if(scaleCorrectionCounter_ > nTriesBeforeScaling_ && scaleParticles_)
	{
		scaleApplication_ = true;
		rescaleRequirement_ = true;
		scalingFactor_++;
		scaleStep_ = pow(scaleStep_,scalingFactor_);
		if(scaleStep_<minScaleFit_)
		{
			scaleStep_=minScaleFit_;
		}
		scaleCorrectionCounter_ = 0;
	}
    
    triSurface triToRet(bodySurfMesh);
    
    return triToRet;
}
//---------------------------------------------------------------------------//
bool addModelRepeatRandomPosition::canAddBody
(
    const volScalarField& body,
    const triSurfaceMesh& bodySurfMesh
)
{
    #include "canAddBodySource.H"
}
// MODEL SPECIFIC FUNCTIONS==================================================//
//---------------------------------------------------------------------------//
void addModelRepeatRandomPosition::initializeCellZone()
{
	
	label zoneID = mesh_.cellZones().findZoneID(zoneName_);
	Info << "-- addModelMessage-- " << "label of the cellZone " << zoneID << endl;
	
	const labelList& cellZoneCells = mesh_.cellZones()[zoneID];
    cellsInBoundBox_[Pstream::myProcNo()] = cellZoneCells;
	
	const pointField& cp = mesh_.C();
	const pointField fCp(cp,cellsInBoundBox_[Pstream::myProcNo()]);
	cellZonePoints_[Pstream::myProcNo()] = fCp;
	
	updateCellZoneBoundBox();
}
//---------------------------------------------------------------------------//
void addModelRepeatRandomPosition::updateCellZoneBoundBox()
{
		boundBox cellZoneBounds(cellZonePoints_[Pstream::myProcNo()]);
		
        reduce(cellZoneBounds.min(), minOp<vector>());
        reduce(cellZoneBounds.max(), maxOp<vector>());
        
        if (Pstream::myProcNo() == 0)
        {
            minBound_ = cellZoneBounds_.min();
            maxBound_ = cellZoneBounds_.max();
            cellZoneBounds_ = boundBox(minBound_,maxBound_);
        }
}
//---------------------------------------------------------------------------//
void addModelRepeatRandomPosition::initializeBoundBox()
{    
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
                    minBound_,maxBound_,bBoxCells
                )
            );
        }
        nextToCheck = auxToCheck;
    }    
    
    cellsInBoundBox_[Pstream::myProcNo()] = bBoxCells[Pstream::myProcNo()];
    
    cellZoneBounds_ = boundBox(minBound_,maxBound_);
}
//---------------------------------------------------------------------------//
labelList addModelRepeatRandomPosition::getBBoxCellsByOctTree
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
scalar addModelRepeatRandomPosition::checkLambdaFraction(const volScalarField& body)
{
	scalarList lambdaIntegrate(Pstream::nProcs());
    scalarList volumeIntegrate(Pstream::nProcs());
	scalar lambdaFraction(0);
    forAll (lambdaIntegrate,k)
    {
        lambdaIntegrate[k] = 0;
        volumeIntegrate[k] = 0;
    }
	forAll (cellsInBoundBox_[Pstream::myProcNo()],k)
	{
		label cell = cellsInBoundBox_[Pstream::myProcNo()][k];
		lambdaIntegrate[Pstream::myProcNo()] += mesh_.V()[cell]*body[cell];
		volumeIntegrate[Pstream::myProcNo()] += mesh_.V()[cell];
	}
	lambdaFraction = gSum(lambdaIntegrate)/gSum(volumeIntegrate);
	Info << "-- addModelMessage-- " << "lambda fraction in controlled region: " << lambdaFraction<< endl;
	return lambdaFraction;
}
//---------------------------------------------------------------------------//
scalar addModelRepeatRandomPosition::returnRandomAngle()
{
    scalar ranNum = 2.0*randGen_.scalar01() - 1.0;
    scalar angle  = ranNum*Foam::constant::mathematical::pi;
	return angle;
}
//---------------------------------------------------------------------------//
scalar addModelRepeatRandomPosition::returnRandomScale()
{
	scalar ranNum       = randGen_.scalar01();
	scalar scaleDiff    = maxScale_ - minScale_; 
    scalar scaleFactor  = minScale_ + ranNum*scaleDiff;
	Info << "-- addModelMessage-- " <<"random scaleFactor " << scaleFactor <<endl;
	return scaleFactor;
}
//---------------------------------------------------------------------------//
vector addModelRepeatRandomPosition::returnRandomRotationAxis()
{
	vector  axisOfRotation(vector::zero);
	scalar ranNum = 0;
    
	for (int i=0;i<3;i++)
	{
		ranNum = randGen_.scalar01();
		axisOfRotation[i] = ranNum;
	}

	axisOfRotation /=mag(axisOfRotation);
	return axisOfRotation;
}
//---------------------------------------------------------------------------//
vector addModelRepeatRandomPosition::returnRandomPosition
(
    const point           CoM,
    const triSurfaceMesh& bodySurfMesh
)
{
    // Note (MI): this function will always return acceptable random
    //            position IF
    //            -> body boundBox is completely inside active boundBox
    //            AND
    //            -> active boundBox is completely contained in the mesh
    //
    // Note (MI): the check if body boundBox is inside active boundBox
    //            is simple and probably unecessary
    // Note (MI): an efficient check if all the active boundBox is inside
    //            mesh is an open issue at the moment
    vector ranVec(vector::zero);
    for (label randomVecAddCounter=0;randomVecAddCounter < 100;randomVecAddCounter++)
    {        
        Info << "-- addModelMessage-- " << "position generation attempt: " << randomVecAddCounter << endl;
        
        // get the bodySurfMesh points
        pointField bSMeshPts = bodySurfMesh.points();
        
        // get the bodySurfMesh boundBox
        const boundBox& bodySurfBounds(bSMeshPts);
        
        // compute the max scales to stay in active bounding box
        vector maxScales(cellZoneBounds_.max() - bodySurfBounds.max());
        maxScales -= cellZoneBounds_.min() - bodySurfBounds.min();
        maxScales *= 0.5*0.9;//0.Y is there just to be sure 
        
        Info << "-- addModelMessage-- " << "acceptable movements: " << maxScales << endl;
        
        scalar ranNum = 0;
        for (int i=0;i<3;i++)
        {
            ranNum = 2.0*maxScales[i]*randGen_.scalar01() - 1.0*maxScales[i];
            ranVec[i] = ranNum;
        }
        
        vector validDirs((mesh_.geometricD() + Vector<label>::one)/2);
        ranVec = cmptMultiply(validDirs,ranVec);//translate only with respect to valid directions
        
        // ok, now I need to check if the generated vector is OK
        // (this is the fun part of the code...)
        // potential fail 1:
        // -> the bodySurfMesh bounding box is NOT fully contained in
        //    active boundBox
        // potential fail 2:
        // -> the active boundBox is NOT fully contained in the mesh
        // if either of these conditions is met, I need to check all the
        // TRANSLATED bodySurfMesh POINTS if they ARE INSIDE the mesh
        
        // 1. check for potential fails 1
        bool    bodySurfBoundsContained(false);
        label   containedDirs(0);
        for (label i=0;i<3;i++)
        {
            if (validDirs[i] < SMALL){containedDirs++;}
            else
            {
                bool maxCheck(cellZoneBounds_.max()[i] - bodySurfBounds.max()[i] > SMALL);
                bool minCheck(cellZoneBounds_.min()[i] - bodySurfBounds.min()[i] < SMALL);
                if (maxCheck and minCheck) {containedDirs++;} 
            }
        }
        if (containedDirs == 3){bodySurfBoundsContained = true;}
        // Note (MI): potential fail 2 is prechecked during addModel
        //            initialization -> bool allActiveCellsInMesh_
        
        if (not bodySurfBoundsContained or not allActiveCellsInMesh_)
        {
            bSMeshPts+=ranVec;
            if (mesh_.nGeometricD() < 3)
            {
                bSMeshPts = cmptMultiply(validDirs,bSMeshPts);
                bSMeshPts += cmptMultiply((vector::one - validDirs),0.5*(mesh_.bounds().max() + mesh_.bounds().min()));
            }
            
            if (not bodySurfBoundsContained)
            {
                Info << "-- addModelMessage-- " << "body boundingBox is NOT fully contained in active boundingBox" << endl;
                Info << "-- addModelMessage-- " << "active boundingBox  : " << cellZoneBounds_ << endl;
                Info << "-- addModelMessage-- " << "surfMesh boundingBox: " << bodySurfBounds << endl;
            }
            else
            {
                Info << "-- addModelMessage-- " << "need to check if all the bodyPoints are inside mesh" << endl;
            }
            meshSearch searchEng(mesh_);
            forAll (bSMeshPts,pointI)
            {
                bool isSurfPtInsideMesh = searchEng.isInside(bSMeshPts[pointI]);
                if (not isSurfPtInsideMesh)
                {
                    Info << "-- addModelMessage-- " << "move resulted in invalid position" << endl;
                    Info << "-- addModelMessage-- " << "discarding vector" << endl;
                    ranVec *= 0.0;
                    break;
                }
            }
        }
        
        if (mag(ranVec) > SMALL) {return ranVec;}
        
    }
        
    return ranVec;
}
