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
stlName_(stlName),
bodyAdded_(false),
mesh_(mesh),

addMode_(word(addModelDict_.lookup("addModel"))),
coeffsDict_(addModelDict_.subDict(addMode_+"Coeffs")),

addDomain_(word(coeffsDict_.lookup("addDomain"))),
scalingMode_(word(coeffsDict_.lookup("scalingMode"))),
rotationMode_(word(coeffsDict_.lookup("rotationMode"))),
addModeI_(word(coeffsDict_.lookup("addMode"))),

addDomainCoeffs_(coeffsDict_.subDict(addDomain_ + "Coeffs")),
scalingModeCoeffs_(coeffsDict_.subDict(scalingMode_ + "Coeffs")),
rotationModeCoeffs_(coeffsDict_.subDict(rotationMode_ + "Coeffs")),
addModeICoeffs_(coeffsDict_.subDict(addModeI_ + "Coeffs")),

//~ useNTimes_(readLabel(timeBasedCoeffs_.lookup("useNTimes"))),
//~ timeBetweenUsage_(readScalar(timeBasedCoeffs_.lookup("timeBetweenUsage"))),
//~ partPerAdd_(readLabel(timeBasedCoeffs_.lookup("partPerAdd"))),
useNTimes_(0),
timeBetweenUsage_(0),
partPerAdd_(0),
partPerAddTemp_(0),
addedOnTimeLevel_(0),
zoneName_(),
minBound_(vector::zero),
maxBound_(vector::zero),
minScale_(0),
maxScale_(0),
minScaleFit_(0),
scaleStep_(0),
nTriesBeforeScaling_(0),
axisOfRot_(vector::zero),
scaleParticles_(false),
randomAxis_(false),
rotateParticles_(false),
iteration_(0),
correctorOfScale_(0),
scaleApplication_(false),
scaleRandomApplication_(false),
rescaleRequirement_(false),
succesfulladition_(false),
scalingFactor_(0),
restartPartCountTemp_(false),
ReapeatedAddition_(false),
firstTimeRunning_(true),
cellZoneActive_(false),
boundBoxActive_(false),
octreeField_(mesh_.nCells(), 0),
timeBased_(false),
fieldBased_(false),
fieldCurrentValue_(0),
fieldValue_(0)
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
    Info << "Debug #11 " << endl;
	if (addModeI_ == "timeBased")
	{
		useNTimes_ = (readLabel(addModeICoeffs_.lookup("useNTimes")));
		timeBetweenUsage_ = (readScalar(addModeICoeffs_.lookup("timeBetweenUsage")));
		partPerAdd_ = (readLabel(addModeICoeffs_.lookup("partPerAdd")));
		timeBased_ = true;
	}
	else if (addModeI_ == "fieldBased")
	{
		fieldValue_ = (readScalar(addModeICoeffs_.lookup("fieldValue")));
		Info << "fieldValue_ " << fieldValue_ << endl;
		fieldBased_ = true;
	}
    else
    {
        Info << "notImplemented, will crash" << endl;
    }
	
	if (addDomain_ == "cellZone")
	{
		zoneName_ = (word(addDomainCoeffs_.lookup("zoneName")));
		cellZoneActive_ = true;
        initializeCellZone();
	}
	else if (addDomain_ == "boundBox")
	{
		minBound_       = (addDomainCoeffs_.lookup("minBound"));
		maxBound_       = (addDomainCoeffs_.lookup("maxBound"));
		boundBoxActive_ = true;
		Info << "maxBound_ "<< maxBound_ << " minBound_ "<< minBound_ << endl;
        
        List<DynamicLabelList> cellsInBoundZone(Pstream::nProcs());
		octreeField_ *= 0;
		bool insideBB = false;
		Info << "Label size is 0 " << cellsInBoundZone.size() << endl;
		
		label cellToCheck = mesh_.findCell(point(minBound_.x(), minBound_.y(), minBound_.z()));
		Info << "Label is" << cellToCheck << endl;
		getBBoxCellsByOctTree(cellToCheck,insideBB,minBound_,maxBound_,cellsInBoundZone);
		cellsInBoundBox_ = cellsInBoundZone[0];
		Info << "Label size is " << cellsInBoundBox_.size() << endl;
	}
	else if (addDomain_ == "domain")
	{
		Info << "notImplemented, will crash" << endl;
	}
    else
    {
		Info << "notImplemented, will crash" << endl;
	}
	
	if (scalingMode_ == "noScaling")
	{
		scaleParticles_ = false;
	}
	else if (scalingMode_ == "randomScaling")
	{
		minScale_               = (readScalar(scalingModeCoeffs_.lookup("minScale")));
		maxScale_               = (readScalar(scalingModeCoeffs_.lookup("maxScale")));
		scaleParticles_         = false;
		scaleRandomApplication_ = true;
	}
	else if (scalingMode_ == "scaleToFit")
	{
		minScaleFit_        = (readScalar(scalingModeCoeffs_.lookup("minScale")));
		scaleStep_          = (readScalar(scalingModeCoeffs_.lookup("scaleStep")));
		nTriesBeforeScaling_= (readScalar(scalingModeCoeffs_.lookup("nTriesBeforeScaling")));
		Info << "nTriesBeforeScaling " << nTriesBeforeScaling_ << endl;
		scaleParticles_ = true;
	}
    else
    {
		Info << "notImplemented, will crash" << endl;
	}
	
	if (rotationMode_ == "noRotation")
	{
		rotateParticles_ = false;
		randomAxis_      = false;
	}
	else if (rotationMode_ == "randomRotation")
	{
		rotateParticles_ = true;
		randomAxis_      = true;
	}
	else if (rotationMode_ == "fixedAxisRandomRotation")
	{
		axisOfRot_       = (rotationModeCoeffs_.lookup("axis"));
		Info << "axisOfRot_ " << axisOfRot_ << endl;
		rotateParticles_ = true;
		randomAxis_      = false;
	}
    else
    {
		Info << "notImplemented, will crash" << endl;
	}
	
	partPerAddTemp_ = partPerAdd_;
	Info << "Debug #12 " << endl;
    
}

//---------------------------------------------------------------------------//
bool addModelRepeatRandomPosition::shouldAddBody(const volScalarField& body)
{
	if (cellZoneActive_ && firstTimeRunning_ )
	{
		//~ initializeCellZone();
		firstTimeRunning_ = false;
		cellZoneActive_ = false;
		addedOnTimeLevel_ = partPerAdd_;
	}
	
	if (boundBoxActive_ && firstTimeRunning_ && fieldBased_)
	{
		//~ List<DynamicLabelList> cellsInBoundZone(Pstream::nProcs());
		//~ octreeField_ *= 0;
		//~ bool insideBB = false;
		//~ Info << "Label size is 0 " << cellsInBoundZone.size() << endl;
		
		//~ label cellToCheck = mesh_.findCell(point(minBound_.x(), minBound_.y(), minBound_.z()));
		//~ Info << "Label is" << cellToCheck << endl;
		//~ getBBoxCellsByOctTree(cellToCheck,insideBB,minBound_,maxBound_,cellsInBoundZone);
		//~ cellsInBoundBox_ = cellsInBoundZone[0];
		//~ Info << "Label size is " << cellsInBoundBox_.size() << endl;
		firstTimeRunning_ = false;
		boundBoxActive_ = false;
	}
	
	if (fieldBased_)
	{
	
		fieldCurrentValue_ = lambdaFraction(body);
		scalar timeVal(mesh_.time().value());
		
		//~ if (timeVal == 0) 
		//~ {
			//~ return false;
		//~ }
		
		Info << "Current Time value is     : " << timeVal << endl;
		Info << "Current lambdaFraction is : " << fieldCurrentValue_ << endl;
		Info << "Set lambdaFraction is : " << fieldValue_ << endl;
		
		return (fieldCurrentValue_ < fieldValue_);
	}
	
	if (timeBased_)
	{
		scalar timeVal(mesh_.time().value());
		scalar deltaTime(mesh_.time().deltaT().value());
		scalar tmFrac(timeVal/timeBetweenUsage_);
		tmFrac -=  floor(tmFrac+deltaTime);
		
		//~ if (timeVal == 0) 
		//~ {
			//~ return false;
		//~ }
		
		Info << "Time/(Time beween usage) - floor(Time/Time beween usage): " 
			 << tmFrac << endl;
			 
		bool tmLevelOk(tmFrac < deltaTime); 
		
		if (not tmLevelOk){addedOnTimeLevel_ = partPerAdd_;}
		
		if (tmLevelOk && addedOnTimeLevel_ > 0 && useNTimes_ > 0)
		{
			Info << " Condition #1 activated " << endl;
			ReapeatedAddition_ = true;
			return true;
		}
		else if (addedOnTimeLevel_> 0 && ReapeatedAddition_)
		{
			return true;
		}
		else
		{
			return (tmLevelOk and useNTimes_ > 0 and addedOnTimeLevel_ == 0 and partPerAddTemp_> 0);
		}
	}
    
}
//---------------------------------------------------------------------------//
triSurface addModelRepeatRandomPosition::addBody
(
    const   volScalarField& body
)
{
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
//====================================================================================================//
	iteration_ = iteration_ +1;
    
	const Vector<label> validDirs = (mesh_.geometricD() + Vector<label>::one)/2;		// Geometric directions
//====================================================================================================//	

	
	if (scaleApplication_ || scaleRandomApplication_)
    {
		if (scaleRandomApplication_)
		{
			scaleStep_ = returnRandomScale();
		}
		Info << "Scaling had happened "<< " scaleStep "<< scaleStep_ <<endl;
		bodySurfMesh.triSurface::scalePoints(scaleStep_);
	}
	
	scalar rotAngel = 0;
	Info << "State of rotateParticlesII is : " << rotateParticles_ << endl;
	if (rotateParticles_)
	{
		rotAngel  = returnRandomAngle();
		
		if (randomAxis_)
		{
			axisOfRot_ = returnRandomRotationAxis();	
		}
		Info << " This is random axis " << axisOfRot_ << endl;
		Info << " This is random angel " << rotAngel << endl;
		Info << " Cosinus of randomAngle " << Foam::cos(rotAngel) << endl;
	    pointField bodyPoints (bodySurfMesh.points());
	    
		tensor rotMatrix(Foam::cos(rotAngel)*tensor::I);
		
		rotMatrix += Foam::sin(rotAngel)*tensor(
         0.0,      -axisOfRot_.z(),  axisOfRot_.y(),
         axisOfRot_.z(), 0.0,       -axisOfRot_.x(),
        -axisOfRot_.y(), axisOfRot_.x(),  0.0
		);
		
		rotMatrix += (1.0-Foam::cos(rotAngel))*(axisOfRot_ * axisOfRot_);
		
		vector  centerOfMass(vector::zero);

		centerOfMass = correctCoM(bodySurfMesh);
		
		forAll (bodyPoints,p)
		{

			if (rotateParticles_)
			{
				bodyPoints[p] -= centerOfMass;
				bodyPoints[p] = rotMatrix & bodyPoints[p];
				bodyPoints[p] += centerOfMass;
			}
		}

		bodySurfMesh.movePoints(bodyPoints);
		Info << " Rotation happend - Debug "<< endl;
	}
	
	vector transIncr = returnRandomPosition(maxBound_,minBound_,validDirs,bodySurfMesh);

	pointField bodyPoints (bodySurfMesh.points());

    forAll (bodyPoints,p)
    {

        bodyPoints[p] += transIncr;
    
    }

    bodySurfMesh.movePoints(bodyPoints);
    
    boundBox bound(bodyPoints);
    //~ minBoundPoint = bound.min();
    //~ maxBoundPoint = bound.max();
    Info << "Bound box of newly moved object " << bound<< endl;
    Info << "Center of newly moved object " << bound.midpoint()<< endl;
//====================================================================================================//
    bool canAddBodyI(canAddBody(body,bodySurfMesh));
    //~ bool boundBoxCondition(isBodyInBoundBox(minBound_,maxBound_,bodySurfMesh));
    reduce(canAddBodyI, andOp<bool>());
    if (!isBodyInBoundBox(maxBound_,minBound_,bodySurfMesh))
    {
		canAddBodyI = false;
	}
    
    bodyAdded_ = (canAddBodyI);
//====================================================================================================//
	//~ if (bodyAdded_ && rescaleRequirement)
	//~ {
		//~ scaleImmersedBodyRandom(1/scaleStep_, bodySurfMesh);
		//~ Info <<" Body has been rescaled to the original size "<<endl;
		//~ rescaleRequirement = false;
	//~ }
	if (!bodyAdded_)
	{
		correctorOfScale_ = correctorOfScale_ +1; 
	}
	
	if (bodyAdded_)
	{			
		if (timeBased_)
		{
			Info << "addedOnTimeLevel_  " << addedOnTimeLevel_<< endl;
			addedOnTimeLevel_--;
			Info<< "bodyAdded_ " << bodyAdded_ << " addedOnTimeLevel_  " << addedOnTimeLevel_<<" useNTimes_ " << useNTimes_<<  endl;
			if (addedOnTimeLevel_ == 0)
			{
				useNTimes_--;
				Info <<" useNTimes_ " << useNTimes_<<  endl;
				ReapeatedAddition_ = false;
			}
		}
		
		correctorOfScale_ = 0;
		
	}
	
	Info << "iteration_  " << iteration_<< endl;
	Info << "correctorOfScale_  " << correctorOfScale_<< endl;
	
	if (correctorOfScale_ > nTriesBeforeScaling_ && scaleParticles_)
	{
		scaleApplication_ = true;
		rescaleRequirement_ = true;
		scalingFactor_++;
		scaleStep_ = pow(scaleStep_,scalingFactor_);
		if (scaleStep_<minScaleFit_)
		{
			scaleStep_=minScaleFit_;
		}
		correctorOfScale_ = 0;
	}
	
    
    Info << "will try to use the body " << useNTimes_ << " more times" << endl;
//====================================================================================================//
    triSurface triToRet(bodySurfMesh);
    
    return triToRet;
}
//====================================================================================================//
vector addModelRepeatRandomPosition::returnRandomPosition
(	
    vector maxBound_,
	vector minBound_,
	vector validDirs,
	const triSurfaceMesh& bodySurfMesh
)
{
	Info << "maxBound_ "<< maxBound_ << " minBound_ "<< minBound_ << endl;
	Random randObj(clock::getTime());
	scalar ranNum = 0;
	for (int k=(0) ;k<6*iteration_; k++) //~ To get always new random number
	{
		ranNum = randObj.scalar01();
	}
	vector  position(vector::zero);
	vector  center(vector::zero);
	vector  range(vector::zero);

	//==================================
	pointField bodyPointsTranslaste (bodySurfMesh.points());
	boundBox littleBounadry (bodyPointsTranslaste);
	vector 	CoM = gSum(bodySurfMesh.coordinates())/bodySurfMesh.size();
	vector  particleDiameter(vector::zero);
	vector  centerOfMass(vector::zero);
	vector  transIncr(vector::zero);

	for (int i = 0; i < 3 ; i = i + 1)
	{
		particleDiameter[i] = mag(littleBounadry.max()[i]-littleBounadry.min()[i]);
		
		if (validDirs[i] == 1)
		{
			centerOfMass[i] = CoM[i];
		}
		else if (validDirs[i] == 0)
		{
			point tempPoint = mesh_.C()[i];
			centerOfMass[i] = tempPoint[i]; 
		}	
		range[i]=(maxBound_[i]-minBound_[i]);
		center[i]=((range[i]/2)+minBound_[i]);
		if (range[i] <=  particleDiameter[i])
		{
			range[i]=0.0;
		}
		ranNum = randObj.scalar01();
		scalar ranNumUp=2*ranNum -1;
		position[i]=(center[i]+ranNumUp*range[i]/2);

		if (validDirs[i] == 1)
		{
				
			if ((particleDiameter[i]/range[i])<0.2)
			{
			 if (range[i] > 0 || range[i] < 0)
			 {
				if (position[i] < (center[i]-(1-(particleDiameter[i]/range[i]))*range[i]/2) && range[i]>particleDiameter[i])
				{
					position[i]=position[i]+0.55*particleDiameter[i];
				}	
				if (position[i] > (center[i]+(1-(particleDiameter[i]/range[i]))*range[i]/2) && range[i]>particleDiameter[i])
				{
					position[i]=position[i]-0.55*particleDiameter[i];
				}
			 }
			}
			if ((particleDiameter[i]/range[i])>0.2)
			{
			 position[i]=(center[i]+ranNumUp*(range[i]/2-0.5*particleDiameter[i]/2));
			 if (range[i] > 0 || range[i] < 0)
			 {
				if (position[i] < (center[i]-(1-(particleDiameter[i]/(range[i])))*range[i]/2) && range[i]>particleDiameter[i])
				{
					position[i]=position[i]+0.25*particleDiameter[i];
				}
					
				if (position[i] > (center[i]+(1-(particleDiameter[i]/(range[i])))*range[i]/2) && range[i]>particleDiameter[i])
				{
					position[i]=position[i]-0.25*particleDiameter[i];
				}
			 }
			}
		}
		transIncr[i] = position[i] - centerOfMass[i];
	}
	Info <<" random vector " << transIncr << endl;
	Info <<" newposition shall be " <<centerOfMass + transIncr << endl;	
	return transIncr;
}
//---------------------------------------------------------------------------//
bool addModelRepeatRandomPosition::isBodyInBoundBox	
(
    vector maxBound_,
    vector minBound_,
    const triSurfaceMesh& bodySurfMesh
)
{
	const Vector<label> validDirs = (mesh_.geometricD() + Vector<label>::one)/2;
	pointField bodyPointsScale (bodySurfMesh.points());								
	boundBox scalingBounadry (bodyPointsScale);
	int outOfTheBox = 0;
	Info << "===============================================" << endl;
	Info << "Particle boundbox" <<scalingBounadry<< endl;
	Info << "Zone boundbox" <<minBound_<<" "<<maxBound_ << endl;
	Info << "===============================================" << endl;
	for (int i = 0;i< 3; i++)
	{
		if (validDirs[i] == 1)
		{
			
			if (scalingBounadry.min()[i]<= minBound_[i])
			{
				outOfTheBox = outOfTheBox + 1;
			}
			if (scalingBounadry.max()[i] >= maxBound_[i])
			{
				outOfTheBox = outOfTheBox + 1;
			}
		}
	}
	if (outOfTheBox == 0)
	{
		return true;
	}
	if (outOfTheBox > 0)
	{
		return false;
	}
}
//---------------------------------------------------------------------------//
vector addModelRepeatRandomPosition::correctCoM(const triSurfaceMesh& bodySurfMesh)
{
	pointField bodyPoints (bodySurfMesh.points());
	const Vector<label> validDirs = (mesh_.geometricD() + Vector<label>::one)/2;
	boundBox littleBounadry (bodyPoints);
	vector  centerOfMass(vector::zero);
	centerOfMass = littleBounadry.midpoint();
	vector CoM = gSum(bodySurfMesh.coordinates())/bodySurfMesh.size();
	//=====================================================================
	Info << " Center of Mass base on BoundBox " << centerOfMass << endl;
	for (int i=0; i<3;i++)
	{	
		if (validDirs[i] == 1)
		{
			centerOfMass[i] = CoM[i];
		}
		else if (validDirs[i] == 0)
		{
			point tempPoint = mesh_.C()[0];
			centerOfMass[i] = tempPoint[i]; 
		}	
	}
	//=====================================================================
	return centerOfMass;
}
//---------------------------------------------------------------------------//
vector addModelRepeatRandomPosition::returnRandomRotationAxis()
{
	vector  axisOfRotation(vector::zero);
	Random randObj(clock::getTime());
	scalar ranNum = 0;
	for (int k=0 ;k<3*iteration_; k++) //~ To get always new random number
	{
		 ranNum = randObj.scalar01();
	}
	for (int i=0;i<3;i++)
	{
		ranNum = randObj.scalar01();
		axisOfRotation[i] = ranNum;
	}

	axisOfRotation /=mag(axisOfRotation);
	return axisOfRotation;
}
//---------------------------------------------------------------------------//
scalar addModelRepeatRandomPosition::returnRandomAngle()
{
	Random randObj(clock::getTime());
	Info << "Script had run " << endl;
	scalar angle = 0;
	scalar ranNum = 0;
	for (int k=0 ;k<iteration_; k++)
	{
		ranNum = randObj.scalar01();
		angle = ranNum*6.28318530718;
	}
	Info <<"random angle " << angle <<endl;
	return angle;
}
//---------------------------------------------------------------------------//
scalar addModelRepeatRandomPosition::returnRandomScale()
{
	Random randObj(clock::getTime());
	Info << "Script had run " << endl;
	scalar scaleFactor = 0;
	scalar ranNum = 0;
	scalar scaleDiff = maxScale_ - minScale_; 
	for (int k=0 ;k<iteration_; k++)
	{
		ranNum = randObj.scalar01();
		scaleFactor = minScale_ + ranNum*scaleDiff;
	}
	Info <<"random scaleFactor " << scaleFactor <<endl;
	return scaleFactor;
}
//---------------------------------------------------------------------------//
void addModelRepeatRandomPosition::initializeCellZone()
{
	
	label zoneID = mesh_.cellZones().findZoneID(zoneName_);
	Info << "label of the cellZone " << zoneID << endl;
	
	const labelList& cellZoneCells = mesh_.cellZones()[zoneID];
		
	Info << "Cells in cellZoneCells size  " << cellZoneCells.size() << ":" << endl;
	Info << "Cells in cellZone  " << zoneName_ << ":" << endl;
	
	cellsInBoundBox_ = cellZoneCells;
	const pointField& cp = mesh_.C();
	const pointField fCp(cp,cellsInBoundBox_);
	cellZonePoints_ = fCp;
	
	updateCellZoneBoundBox();
}
//---------------------------------------------------------------------------//
void addModelRepeatRandomPosition::updateCellZoneBoundBox()
{
		//~ Info << "=================================================================="<< endl;
		boundBox cellZoneBounds (cellZonePoints_);
		const Vector<label> validDirs = (mesh_.geometricD() + Vector<label>::one)/2;
        
        cellZoneBounds_ = cellZoneBounds;
        minBound_ = cellZoneBounds_.min();
		maxBound_ = cellZoneBounds_.max();
        
        // TO BE CORRECTED LATER - IT DOES NOT WORK ANY WAY
		//~ vector correctBbox(vector::zero);
		//~ for (int i=0; i<3;i++)
		//~ {
			//~ if (validDirs[i] == 1)
			//~ {
				//~ correctBbox[i] = 0.001;
			//~ }
			//~ else if (validDirs[i] == 0)
			//~ {
				//~ correctBbox[i] = 0; 
			//~ }	
		//~ }
		//~ Info << "cellZoneBounds size  "<< cellZoneBounds << " Corrector Vector   " << correctBbox <<endl;
		//~ Info << mesh_.findCell(point(cellZoneBounds.min().x(), cellZoneBounds.min().y(), cellZoneBounds.min().z())) << " <-----   Return of cell boundBox minimum"<< endl;
		//~ Info << mesh_.findCell(point(cellZoneBounds.max().x(), cellZoneBounds.max().y(),cellZoneBounds.max().z())) << " <-----   Return of cell boundBox maximum"<< endl;
		//~ int iterator = 0;
		//~ int iteratorII = 0;
		//~ int iteratorIII = 0;
		//~ int iteratorIV = 0;
		//~ while ((mesh_.findCell(point(cellZoneBounds.min().x(), cellZoneBounds.min().y(), cellZoneBounds.min().z())) == -1))
		//~ {
			//~ iterator = iterator + 1; 
			//~ if (mesh_.findCell(point(cellZoneBounds.min().x(), cellZoneBounds.min().y(), cellZoneBounds.min().z())) == -1)
			//~ {
				//~ cellZoneBounds.min().x() = cellZoneBounds.min().x() + correctBbox.x();
			//~ }
			//~ if (mesh_.findCell(point(cellZoneBounds.min().x(), cellZoneBounds.min().y(), cellZoneBounds.min().z())) == -1)
			//~ {
				//~ cellZoneBounds.min().y() = cellZoneBounds.min().y() + correctBbox.y();
			//~ }
			//~ if (mesh_.findCell(point(cellZoneBounds.min().x(), cellZoneBounds.min().y(), cellZoneBounds.min().z())) == -1)
			//~ {
				//~ cellZoneBounds.min().z() = cellZoneBounds.min().z() + correctBbox.z();
			//~ }
		//~ }
		//~ while ((mesh_.findCell(point(cellZoneBounds.max().x(), cellZoneBounds.max().y(), cellZoneBounds.max().z())) == -1))
		//~ {
			//~ iteratorII = iteratorII + 1;
			//~ if (mesh_.findCell(point(cellZoneBounds.max().x(), cellZoneBounds.max().y(), cellZoneBounds.max().z())) == -1)
			//~ {
				//~ cellZoneBounds.max().x() = cellZoneBounds.max().x() - correctBbox.x();
			//~ }
			//~ if (mesh_.findCell(point(cellZoneBounds.max().x(), cellZoneBounds.max().y(), cellZoneBounds.max().z())) == -1)
			//~ {
				//~ cellZoneBounds.max().y() = cellZoneBounds.max().y() - correctBbox.y();
			//~ }
			//~ if (mesh_.findCell(point(cellZoneBounds.max().x(), cellZoneBounds.max().y(), cellZoneBounds.max().z())) == -1)
			//~ {
				//~ cellZoneBounds.max().z() = cellZoneBounds.max().z() - correctBbox.z();
			//~ }
		//~ }
		//~ while ((mesh_.findCell(point(cellZoneBounds.max().x(), cellZoneBounds.min().y(), cellZoneBounds.min().z())) == -1))
		//~ {
			//~ iteratorIII = iteratorIII + 1;
			//~ if (mesh_.findCell(point(cellZoneBounds.max().x(), cellZoneBounds.min().y(), cellZoneBounds.min().z())) == -1)
			//~ {
				//~ cellZoneBounds.max().x() = cellZoneBounds.max().x() - correctBbox.x();
			//~ }
			//~ if (mesh_.findCell(point(cellZoneBounds.max().x(), cellZoneBounds.min().y(), cellZoneBounds.min().z())) == -1)
			//~ {
				//~ cellZoneBounds.min().y() = cellZoneBounds.min().y() + correctBbox.y();
			//~ }
			//~ if (mesh_.findCell(point(cellZoneBounds.max().x(), cellZoneBounds.min().y(), cellZoneBounds.min().z())) == -1)
			//~ {
				//~ cellZoneBounds.min().z() = cellZoneBounds.min().z() - correctBbox.z();
			//~ }
		//~ }
		//~ while ((mesh_.findCell(point(cellZoneBounds.min().x(), cellZoneBounds.max().y(), cellZoneBounds.max().z())) == -1))
		//~ {
			//~ iteratorIV = iteratorIV + 1;
			//~ if (mesh_.findCell(point(cellZoneBounds.min().x(), cellZoneBounds.max().y(), cellZoneBounds.max().z())) == -1)
			//~ {
				//~ cellZoneBounds.min().x() = cellZoneBounds.min().x() + correctBbox.x();
			//~ }
			//~ if (mesh_.findCell(point(cellZoneBounds.min().x(), cellZoneBounds.max().y(), cellZoneBounds.max().z())) == -1)
			//~ {
				//~ cellZoneBounds.max().y() = cellZoneBounds.max().y() - correctBbox.y();
			//~ }
			//~ if (mesh_.findCell(point(cellZoneBounds.min().x(), cellZoneBounds.max().y(), cellZoneBounds.max().z())) == -1)
			//~ {
				//~ cellZoneBounds.min().z() = cellZoneBounds.min().z() - correctBbox.z();
			//~ }
		//~ }
		//~ cellZoneBounds_ = cellZoneBounds;
		//~ Info << "cellZoneBounds size  "<< cellZoneBounds_ << " Corrector Vector   " << correctBbox << " iterator : " << iterator << " iteratorII : "<< iteratorII << " iteratorIII : " << iteratorIII << " iteratorIV : "<< iteratorIV << endl;
		//~ Info << mesh_.findCell(point(cellZoneBounds_.min().x(), cellZoneBounds_.min().y(), cellZoneBounds_.min().z())) << " <-----   Return of cell boundBox minimum"<< endl;
		//~ Info << mesh_.findCell(point(cellZoneBounds_.max().x(), cellZoneBounds_.max().y(),cellZoneBounds_.max().z())) << " <-----   Return of cell boundBox maximum"<< endl;
		//~ Info << "PointField size  "<< cellZonePoints_.size()<<endl;
		//~ Info << "cellZoneBounds size  "<< cellZoneBounds_ <<endl;
		//~ minBound_ = cellZoneBounds_.min();
		//~ maxBound_ = cellZoneBounds_.max();
		//~ Info << "minBound_ size  "<< minBound_<<endl;
		//~ Info << "maxBound_ size  "<< maxBound_ <<endl;
		//~ Info << "=================================================================="<< endl;
}
//---------------------------------------------------------------------------//
void addModelRepeatRandomPosition::getBBoxCellsByOctTree
	(
		label cellToCheck,
		bool& insideBB,
		vector& bBoxMin,
		vector& bBoxMax,
		List<DynamicLabelList>& bBoxCells
	)
{
    //~ Info << "Descended into octree search" << endl;
    //~ Info << bBoxCells << endl;
    if (octreeField_[cellToCheck] ==0)
    {
        octreeField_[cellToCheck] = 1;
        vector cCenter = mesh_.C()[cellToCheck];
        label   partCheck(0);
        forAll (bBoxMin,vecI)
        {
            if (cCenter[vecI] >= bBoxMin[vecI] && cCenter[vecI] <= bBoxMax[vecI])
            {
                partCheck += 1;
            }
        }
        bool cellInside = (partCheck == 3) ? true : false;
        if (cellInside)
        {
            bBoxCells[Pstream::myProcNo()].append(cellToCheck);
            insideBB = true;
        }
        if (!insideBB || cellInside)
        {
            labelList cellNb(mesh_.cellCells()[cellToCheck]);
            forAll (cellNb,nbCellI)
            {
                getBBoxCellsByOctTree(cellNb[nbCellI],insideBB,bBoxMin,bBoxMax,bBoxCells);
            }
        }
    }
}
//---------------------------------------------------------------------------//
scalar addModelRepeatRandomPosition::lambdaFraction(const volScalarField& body)
{
	//~ scalar activeCells = 0;
	//~ scalar zoneVolume = 0;
	//~ scalar activeFraction = 0;
	//~ forAll (cellsInBoundBox_,k)
	//~ {
		//~ label cell = cellsInBoundBox_[k];
		//~ activeCells += mesh_.V()[cell]*body[cell];
		//~ zoneVolume += mesh_.V()[cell];
	//~ }
    Info << cellsInBoundBox_ << endl;
    scalarField fBody(body,cellsInBoundBox_);
    scalarField fVols(mesh_.V(),cellsInBoundBox_);
    scalar activeFraction = gSum(fBody*fVols)/gSum(fVols);
	//~ activeFraction = activeCells/zoneVolume;
    //~ activeFraction = gSum(mesh_.V()[cellsInBoundBox_]*body[cellsInBoundBox_])/gSum(mesh_.V()[cellsInBoundBox_]);
	Info << "Active Fraction " << activeFraction<< endl;
	return activeFraction;
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
