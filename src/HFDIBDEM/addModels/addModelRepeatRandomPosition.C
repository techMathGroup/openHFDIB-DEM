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
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "addModelRepeatRandomPosition.H"
#include "meshSearch.H"

using namespace Foam;

//---------------------------------------------------------------------------//
addModelRepeatRandomPosition::addModelRepeatRandomPosition
(
    const dictionary& addModelDict,
    const Foam::fvMesh& mesh,
    std::unique_ptr<geomModel> bodyGeomModel,
    List<labelList>& cellPoints
)
:
addModel(mesh, std::move(bodyGeomModel), cellPoints),
addModelDict_(addModelDict),
addMode_(word(addModelDict_.lookup("addModel"))),
bodyAdded_(false),
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
void addModelRepeatRandomPosition::init()
{
    // set sizes to necessary datatypes
    cellsInBoundBox_.setSize(Pstream::nProcs());
    cellZonePoints_.setSize(Pstream::nProcs());


	if (addModeI_ == "timeBased")
	{
        useNTimes_ = (readLabel(addModeICoeffs_.lookup("useNTimes")));
		timeBetweenUsage_ = (readScalar(addModeICoeffs_.lookup("timeBetweenUsage")));
		partPerAdd_ = (readLabel(addModeICoeffs_.lookup("partPerAdd")));
        timeBased_ = true;
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "addModel will control simulation time" << endl;
        InfoH << "-- addModelMessage-- " << "STL will be re-used "
            << useNTimes_ << " times" << endl;
        InfoH << "-- addModelMessage-- " << "STL will be added each "
            << timeBetweenUsage_ << " [T]" << endl;
        InfoH << "-- addModelMessage-- " << "upon each addition, "
            << partPerAdd_ << " bodies will be generated from the given STL"
            << endl;
	}
	else if (addModeI_ == "fieldBased")
	{
		fieldValue_ = (readScalar(addModeICoeffs_.lookup("fieldValue")));
        fieldBased_ = true;
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "addModel will control particles volume fraction" << endl;
		InfoH << "-- addModelMessage-- "
            << "preset volume fraction: " << fieldValue_ << endl;
	}
    else
    {
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "notImplemented, will crash" << endl;
    }

	if (addDomain_ == "cellZone")
	{
		zoneName_ = (word(addDomainCoeffs_.lookup("zoneName")));
		cellZoneActive_ = true;
        initializeCellZone();
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "cellZone based addition zone" << endl;
	}
	else if (addDomain_ == "boundBox")
	{
		minBound_       = vector(addDomainCoeffs_.lookup("minBound"));
		maxBound_       = vector(addDomainCoeffs_.lookup("maxBound"));
		boundBoxActive_ = true;
        initializeBoundBox();
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "boundBox based addition zone" << endl;
	}
	else if (addDomain_ == "domain")
	{
		InfoH << addModel_Info << "-- addModelMessage-- "
            << "notImplemented, will crash" << endl;
	}
    else
    {
		InfoH << addModel_Info << "-- addModelMessage-- "
            << "notImplemented, will crash" << endl;
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
        InfoH << addModel_Info << "-- addModelMessage-- "
             << "addition zone NOT completely immersed in mesh "
             << "this computation will be EXPENSIVE" << endl;
        InfoH << zoneVol << " " << zoneBBoxVol << endl;
    }
    else
    {
        InfoH << addModel_Info << "-- addModelMessage-- "
             << "addition zone completely immersed in mesh -> OK" << endl;
    }

	if (scalingMode_ == "noScaling")
	{
		scaleParticles_ = false;
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "all particles will have the same scale" << endl;
	}
	else if (scalingMode_ == "randomScaling")
	{
		minScale_               = (readScalar(scalingModeCoeffs_.lookup("minScale")));
		maxScale_               = (readScalar(scalingModeCoeffs_.lookup("maxScale")));
		scaleParticles_         = false;
		scaleRandomApplication_ = true;
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "particles will be randomly scaled" << endl;
	}
	else if (scalingMode_ == "scaleToFit")
	{
		minScaleFit_        = (readScalar(scalingModeCoeffs_.lookup("minScale")));
		scaleStep_          = (readScalar(scalingModeCoeffs_.lookup("scaleStep")));
		nTriesBeforeScaling_= (readScalar(scalingModeCoeffs_.lookup("nTriesBeforeScaling")));
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "particles will be downscaled to better fill the domain" << endl;
		InfoH << "-- addModelMessage-- "
            << "nTriesBeforeDownScaling: " << nTriesBeforeScaling_ << endl;
		scaleParticles_ = true;
	}
    else
    {
		InfoH << addModel_Info << "-- addModelMessage-- "
            << "notImplemented, will crash" << endl;
	}

	if (rotationMode_ == "noRotation")
	{
		rotateParticles_ = false;
		randomAxis_      = false;
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "source STL will not be rotated upon addition" << endl;
	}
	else if (rotationMode_ == "randomRotation")
	{
		rotateParticles_ = true;
		randomAxis_      = true;
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "source STL will be randomly rotated upon addition" << endl;
	}
	else if (rotationMode_ == "fixedAxisRandomRotation")
	{
		axisOfRot_       = vector(rotationModeCoeffs_.lookup("axis"));
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "source STL will be rotated by a random angle around a fixed axis upon addition" << endl;
		InfoH << "-- addModelMessage-- " << "set rotation axis: "
            << axisOfRot_ << endl;
		rotateParticles_ = true;
		randomAxis_      = false;
	}
    else
    {
		InfoH << addModel_Info << "-- addModelMessage-- "
            << "notImplemented, will crash" << endl;
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

        InfoH << addModel_Info << "-- addModelMessage-- "
            << "Time/(Time beween usage) - floor(Time/Time beween usage): "
            << tmFrac << endl;

        InfoH << addModel_Info << "-- addModelMessage-- "
            << "Number of bodies added on this time level: "
            << addedOnTimeLevel_ << endl;

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
            InfoH << addModel_Info << "-- addModelMessage-- "
                << "Current lambda fraction = " << currentLambdaFrac
                << " < then preset lambda fraction = " << fieldValue_ << endl;
            return true;
        }
    }

    return false;

}
//---------------------------------------------------------------------------//
std::shared_ptr<geomModel> addModelRepeatRandomPosition::addBody
(
    const volScalarField& body,
    PtrList<immersedBody>& immersedBodies
)
{
    geomModel_->resetBody();

    bodyAdditionAttemptCounter_++;

    // rotate
    if (rotateParticles_)
    {
        scalar rotAngle = returnRandomAngle();
        if (randomAxis_)
        {
            axisOfRot_ = returnRandomRotationAxis();
        }
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "Will rotate by " << rotAngle
            << " PiRad around axis " << axisOfRot_ << endl;

        geomModel_->bodyRotatePoints(rotAngle,axisOfRot_);
    }

    // scale
    if (scaleApplication_ or scaleRandomApplication_)
    {
        if (scaleRandomApplication_){scaleStep_ = returnRandomScale();}
        geomModel_->bodyScalePoints(scaleStep_);
    }

    geomModel_->bodyScalePoints(1.02);

    vector CoM(geomModel_->getCoM());
    point bBoxCenter = cellZoneBounds_.midpoint();
    geomModel_->bodyMovePoints(bBoxCenter - CoM);

    // translate

    vector randomTrans = geomModel_->addModelReturnRandomPosition(allActiveCellsInMesh_,cellZoneBounds_,randGen_);
    geomModel_->bodyMovePoints(randomTrans);

    // check if the body can be added
    volScalarField helpBodyField_ = body;
    geomModel_->createImmersedBody(
        helpBodyField_,
        octreeField_,
        cellPoints_
    );

    bool canAddBodyI = !isBodyInContact(immersedBodies);

    geomModel_->bodyScalePoints(1.0/1.02);

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
			InfoH << addModel_Info << "-- addModelMessage-- "
                << "addedOnTimeLevel:  " << addedOnTimeLevel_<< endl;
			addedOnTimeLevel_++;
			InfoH << "-- addModelMessage-- " << "bodyAdded: "
                << bodyAdded_ << " addedOnTimeLevel:  " << addedOnTimeLevel_
                << " useNTimes: " << useNTimes_<<  endl;
			if(addedOnTimeLevel_ == partPerAdd_)
			{
				useNTimes_--;
				InfoH << "-- addModelMessage-- "
                    <<" useNTimes: " << useNTimes_<<  endl;
				reapeatedAddition_ = false;
			}
		}

		scaleCorrectionCounter_ = 0;

	}

	InfoH << addModel_Info << "-- addModelMessage-- "
        << "bodyAdditionAttemptNr  : " << bodyAdditionAttemptCounter_<< endl;
	InfoH << "-- addModelMessage-- " << "sameScaleAttempts      : "
        << scaleCorrectionCounter_<< endl;

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

    return geomModel_->getCopy();
}
// MODEL SPECIFIC FUNCTIONS==================================================//
//---------------------------------------------------------------------------//
void addModelRepeatRandomPosition::initializeCellZone()
{

	label zoneID = mesh_.cellZones().findZoneID(zoneName_);
	InfoH << addModel_Info << "-- addModelMessage-- "
        << "label of the cellZone " << zoneID << endl;

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
void addModelRepeatRandomPosition::recreateBoundBox()
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

    InfoH << addModel_Info << "-- addModelMessage-- "
        << "recreated boundBox size "
        << cellsInBoundBox_[Pstream::myProcNo()].size() << endl;

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
    InfoH << addModel_Info << "-- addModelMessage-- "
          << "lambda fraction in controlled region: " << lambdaFraction<< endl;
    return lambdaFraction;
}
//---------------------------------------------------------------------------//
scalar addModelRepeatRandomPosition::returnRandomAngle()
{
    scalar ranNum = 2.0*randGen_.sample01<scalar>() - 1.0;
    scalar angle  = ranNum*Foam::constant::mathematical::pi;
    return angle;
}
//---------------------------------------------------------------------------//
scalar addModelRepeatRandomPosition::returnRandomScale()
{
    scalar ranNum       = randGen_.sample01<scalar>();
    scalar scaleDiff    = maxScale_ - minScale_;
    scalar scaleFactor  = minScale_ + ranNum*scaleDiff;
    InfoH << addModel_Info << "-- addModelMessage-- "
          <<"random scaleFactor " << scaleFactor <<endl;
    return scaleFactor;
}
//---------------------------------------------------------------------------//
vector addModelRepeatRandomPosition::returnRandomRotationAxis()
{
    vector  axisOfRotation(vector::zero);
    scalar ranNum = 0;

    for (int i=0;i<3;i++)
    {
        ranNum = randGen_.sample01<scalar>();
        axisOfRotation[i] = ranNum;
    }

    axisOfRotation /=mag(axisOfRotation);
    return axisOfRotation;
}
//---------------------------------------------------------------------------//
