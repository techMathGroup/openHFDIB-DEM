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
#include "addModelDistribution.H"
#include "meshSearch.H"

using namespace Foam;

//---------------------------------------------------------------------------//
addModelDistribution::addModelDistribution
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

distributionDict_
(
    IOobject
    (
        "distributionDict",
        "constant",
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
distribution_(scalarList(distributionDict_.lookup("distribution"))),
particleSize_(scalarList(distributionDict_.lookup("particleSize"))),
convertToMeters_(readScalar(distributionDict_.lookup("convertToMeters"))),
volumeOfAddedBodies_(0),

coeffsDict_(addModelDict_.subDict(addMode_+"Coeffs")),
stlBaseSize_(readScalar(coeffsDict_.lookup("stlBaseSize"))),

addDomain_(word(coeffsDict_.lookup("addDomain"))),
addModeI_(word(coeffsDict_.lookup("addMode"))),

addDomainCoeffs_(coeffsDict_.subDict(addDomain_ + "Coeffs")),
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

bodyAdditionAttemptCounter_(0),

succesfulladition_(false),
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
nGeometricD_(0),
geometricD_(Vector<label>::one),
randGen_(clock::getTime())
{
	init();
}
    
addModelDistribution::~addModelDistribution()
{
}

//---------------------------------------------------------------------------//
// Note (MI): Initialization of this model is so complex that I moved it to a
//            specific function
void addModelDistribution::init()
{
    // Set sizes to necessary datatypes
    cellsInBoundBox_.setSize(Pstream::nProcs());
    cellZonePoints_.setSize(Pstream::nProcs());
    addedParticlesSize_.setSize(particleSize_.size(), 0);
    

	if (addModeI_ == "timeBased")
	{
        Info << "-- addModelMessage-- " << "notImplemented, will crash" << endl;
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
                if (addDomainCoeffs_.found("nGeometricD"))
        {
            nGeometricD_ = readLabel(addDomainCoeffs_.lookup("nGeometricD"));
        }
        else
        {
            nGeometricD_ = mesh_.nGeometricD();
        }
        if (addDomainCoeffs_.found("geometricD"))
        {
            geometricD_ = addDomainCoeffs_.lookup("geometricD");
        }
        else
        {
            geometricD_ = mesh_.geometricD();
        }
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
    
    Pstream::gatherList(procZoneVols, 0);
    Pstream::scatter(procZoneVols, 0);
    
    scalar zoneVol(0);
    forAll (procZoneVols, procI)
    {
        zoneVol += procZoneVols[procI];
    }
    
//     scalar zoneVol(gSum(procZoneVols));
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
	
	partPerAddTemp_ = partPerAdd_;
    
}

//---------------------------------------------------------------------------//
bool addModelDistribution::shouldAddBody(const volScalarField& body)
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
triSurface addModelDistribution::addBody
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
    
    Tuple2<label, scalar> scaleFactor = returnScaleFactor();
    Info << "-- addModelMessage-- " << "scaled STL size: " << stlBaseSize_ * scaleFactor.second() << endl;
    bodySurfMesh.scalePoints(scaleFactor.second());
    scalar partVolume(1.0/6.0*3.14*pow(stlBaseSize_ * scaleFactor.second(),3));
    
    // get the working bounding box center
    point bBoxCenter = cellZoneBounds_.midpoint();    
    pointField bodyPoints(bodySurfMesh.points());
    
    // get its center of mass
    vector CoM(vector::zero);
    forAll(bodyPoints,point)
    {
        CoM += bodyPoints[point];
    }    
    CoM/= bodyPoints.size();
    Info << "-- addModelMessage-- " << "scaled STL CoM: " << CoM << endl;
    
    scalar rotAngle = returnRandomAngle();

    vector axisOfRot = returnRandomRotationAxis();
    
    tensor rotMatrix(Foam::cos(rotAngle)*tensor::I);
    
    rotMatrix += Foam::sin(rotAngle)*tensor(
            0.0,      -axisOfRot.z(),  axisOfRot.y(),
            axisOfRot.z(), 0.0,       -axisOfRot.x(),
        -axisOfRot.y(), axisOfRot.x(),  0.0
    );
    
    rotMatrix += (1.0-Foam::cos(rotAngle))*(axisOfRot * axisOfRot);
    
    bodyPoints -= CoM;
    bodyPoints = rotMatrix & bodyPoints;
    bodyPoints += CoM;
    
    bodySurfMesh.movePoints(bodyPoints);
    
    CoM = vector::zero;
    forAll(bodyPoints,point)
    {
        CoM += bodyPoints[point];
    }    
    CoM/= bodyPoints.size();
    
    // move the body to the center of the bounding box
    bodyPoints += bBoxCenter-CoM;
//     CoM         = bBoxCenter;
    bodySurfMesh.movePoints(bodyPoints);
    CoM = vector::zero;
    forAll(bodyPoints,point)
    {
        CoM += bodyPoints[point];
    }    
    CoM/= bodyPoints.size();
    Info << "-- addModelMessage-- " << "moved STL CoM: " << CoM << endl;
    
    // translate
    bodyPoints = bodySurfMesh.points();
    if (nGeometricD_ < 3)
    {
        const vector validDirs = (geometricD_ + Vector<label>::one)/2;
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
		
		volumeOfAddedBodies_ += partVolume;
        addedParticlesSize_[scaleFactor.first()] += partVolume;
	}
	
	Info << "-- addModelMessage-- " << "bodyAdditionAttemptNr  : " << bodyAdditionAttemptCounter_<< endl;
    
    triSurface triToRet(bodySurfMesh);
    
    return triToRet;
}
//---------------------------------------------------------------------------//
bool addModelDistribution::canAddBody
(
    const volScalarField& body,
    const triSurfaceMesh& bodySurfMesh
)
{
    #include "canAddBodySource.H"
}
// MODEL SPECIFIC FUNCTIONS==================================================//
//---------------------------------------------------------------------------//
void addModelDistribution::initializeCellZone()
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
void addModelDistribution::updateCellZoneBoundBox()
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
void addModelDistribution::initializeBoundBox()
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
    
    Info << "-- addModelMessage-- " << "initiliazed boundBox size " << cellsInBoundBox_[Pstream::myProcNo()].size() << endl;
    
    cellZoneBounds_ = boundBox(minBound_,maxBound_);
}
//---------------------------------------------------------------------------//
void addModelDistribution::recreateBoundBox()
{    
    octreeField_ = Field<label>(mesh_.nCells(), 0);
    cellsInBoundBox_[Pstream::myProcNo()].clear();
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
    
    Info << "-- addModelMessage-- " << "recreated boundBox size " << cellsInBoundBox_[Pstream::myProcNo()].size() << endl;
    
    cellZoneBounds_ = boundBox(minBound_,maxBound_);
}
//---------------------------------------------------------------------------//
labelList addModelDistribution::getBBoxCellsByOctTree
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
scalar addModelDistribution::checkLambdaFraction(const volScalarField& body)
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
scalar addModelDistribution::returnRandomAngle()
{
    scalar ranNum = 2.0*randGen_.scalar01() - 1.0;
    scalar angle  = ranNum*Foam::constant::mathematical::pi;
	return angle;
}
//---------------------------------------------------------------------------//
vector addModelDistribution::returnRandomRotationAxis()
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
Tuple2<label, scalar> addModelDistribution::returnScaleFactor()
{
    DynamicScalarList  distributionDiff;
    forAll (addedParticlesSize_,size)
    {
        distributionDiff.append(distribution_[size] - 100*addedParticlesSize_[size]/(volumeOfAddedBodies_+SMALL));
    }
    
    label highestDiff(0);
    forAll (distributionDiff,size)
    {
        if(distributionDiff[size] > distributionDiff[highestDiff])
        {
            highestDiff = size;
        }
    }
    
    scalar factor(particleSize_[highestDiff - 1] + (particleSize_[highestDiff] - particleSize_[highestDiff - 1]) * randGen_.scalar01());
    factor *= convertToMeters_/stlBaseSize_;
    
    Tuple2<label, scalar> returnValue(highestDiff, factor);
    
    return returnValue;
}
//---------------------------------------------------------------------------//
vector addModelDistribution::returnRandomPosition
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
    
    meshSearch searchEng(mesh_);
    pointField bSMeshPts = bodySurfMesh.points();
    for (label randomVecAddCounter=0;randomVecAddCounter < 100;randomVecAddCounter++)
    {        
        Info << "-- addModelMessage-- " << "position generation attempt: " << randomVecAddCounter << endl;
        // get the bodySurfMesh points
        bSMeshPts = bodySurfMesh.points();
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
        
        vector validDirs((geometricD_ + Vector<label>::one)/2);
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
            if (nGeometricD_ < 3)
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
            
            forAll (bSMeshPts,pointI)
            {
                bool isSurfPtInsideMesh = searchEng.isInside(bSMeshPts[pointI]);
                reduce(isSurfPtInsideMesh, orOp<bool>());
                if (not isSurfPtInsideMesh)
                {                    
                    Info << "-- addModelMessage-- " << "move resulted in invalid position. Point not in mesh: " << bSMeshPts[pointI] << endl;
                    Info << "-- addModelMessage-- " << "discarding vector: " << ranVec << endl;
                    ranVec *= 0.0;
                    break;
                }
            }
        }
        
        if (mag(ranVec) > SMALL) {return ranVec;}
        
    }
        
    return ranVec;
}
