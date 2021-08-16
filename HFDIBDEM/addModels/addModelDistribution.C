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
    const Foam::dynamicFvMesh& mesh,
    const Vector<label> geomDir,
    geomModel* bodyGeomModel
)
:
addModel(mesh),
addModelDict_(addModelDict),
addMode_(word(addModelDict_.lookup("addModel"))),
bodyAdded_(false),
geomModel_(bodyGeomModel),
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
geometricD_(geomDir),
randGen_(clock::getTime())
{
	init();
}

addModelDistribution::~addModelDistribution()
{
}

//---------------------------------------------------------------------------//
void addModelDistribution::init()
{
    // set sizes to necessary datatypes
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

	partPerAddTemp_ = partPerAdd_;

	forAll (geometricD_, direction)
    {
        if (geometricD_[direction] == 1)
        {
            nGeometricD_++;
        }
    }
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
geomModel* addModelDistribution::addBody
(
    const   volScalarField& body
)
{
    bodyAdditionAttemptCounter_++;

    geomModel_->resetBody();

    Tuple2<label, scalar> scaleFactor = returnScaleFactor();
    Info << "-- addModelMessage-- " << "scaled STL size: " << stlBaseSize_ * scaleFactor.second() << endl;
    geomModel_->bodyScalePoints(scaleFactor.second());
    scalar partVolume(1.0/6.0*3.14*pow(stlBaseSize_ * scaleFactor.second(),3));

    scalar rotAngle = returnRandomAngle();

    vector axisOfRot = returnRandomRotationAxis();

    geomModel_->bodyRotatePoints(rotAngle,axisOfRot);

    vector CoM(geomModel_->getCoM());
    point bBoxCenter = cellZoneBounds_.midpoint();
    geomModel_->bodyMovePoints(bBoxCenter - CoM);

    vector randomTrans = geomModel_->addModelReturnRandomPosition(allActiveCellsInMesh_,cellZoneBounds_,randGen_);
    geomModel_->bodyMovePoints(randomTrans);
    // check if the body can be added
    bool canAddBodyI(geomModel_->canAddBody(body));
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

    return geomModel_->getGeomModel();;
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
