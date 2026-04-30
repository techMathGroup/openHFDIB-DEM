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
#include "addModelDistribution.H"
#include "meshSearch.H"

using namespace Foam;

//---------------------------------------------------------------------------//
addModelDistribution::addModelDistribution
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
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "notImplemented, will crash" << endl;
	}
	else if (addModeI_ == "fieldBased")
	{
		fieldValue_ = (readScalar(addModeICoeffs_.lookup("fieldValue")));
        fieldBased_ = true;
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "addModel will control particles volume fraction" << endl;
		InfoH << "-- addModelMessage-- " << "preset volume fraction: "
            << fieldValue_ << endl;
	}
    else
    {
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "notImplemented, will crash" << endl;
    }

	if (addDomain_ == "cellZone")
	{
		zoneName_ = word(addDomainCoeffs_.lookup("zoneName"));
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

    Pstream::gatherList(procZoneVols, 0);
    Pstream::broadcast(procZoneVols, 0);                                //OF.com compatibility: scatter -> broadcast


    scalar zoneVol(0);
    forAll (procZoneVols, procI)
    {
        zoneVol += procZoneVols[procI];
    }

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

        InfoH << addModel_Info << "-- addModelMessage-- "
            << "Time/(Time beween usage) - floor(Time/Time beween usage): "
            << tmFrac << endl;

        InfoH << "-- addModelMessage-- "
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
std::shared_ptr<geomModel> addModelDistribution::addBody
(
    const volScalarField& body,
    PtrList<immersedBody>& immersedBodies
)
{
    bodyAdditionAttemptCounter_++;

    geomModel_->resetBody();

    Tuple2<label, scalar> scaleFactor = returnScaleFactor();
    InfoH << addModel_Info << "-- addModelMessage-- "
        << "scaled STL size: " << stlBaseSize_ * scaleFactor.second() << endl;
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
    volScalarField helpBodyField_ = body;
    geomModel_->createImmersedBody(
        helpBodyField_,
        octreeField_,
        cellPoints_
    );

    bool canAddBodyI = !isBodyInContact(immersedBodies);

    reduce(canAddBodyI, andOp<bool>());
    bodyAdded_ = (canAddBodyI);

    if(bodyAdded_)
    {
        if(timeBased_)
		{
			InfoH << addModel_Info << "-- addModelMessage-- "
                << "addedOnTimeLevel:  " << addedOnTimeLevel_<< endl;
			addedOnTimeLevel_++;
			InfoH << addModel_Info << "-- addModelMessage-- " << "bodyAdded: "
                << bodyAdded_ << " addedOnTimeLevel:  " << addedOnTimeLevel_
                << " useNTimes: " << useNTimes_<<  endl;
			if(addedOnTimeLevel_ == partPerAdd_)
			{
				useNTimes_--;
				InfoH << addModel_Info << "-- addModelMessage-- "
                    << " useNTimes: " << useNTimes_<<  endl;
				reapeatedAddition_ = false;
			}
		}

		volumeOfAddedBodies_ += partVolume;
        addedParticlesSize_[scaleFactor.first()] += partVolume;
	}

	InfoH << addModel_Info << "-- addModelMessage-- "
        << "bodyAdditionAttemptNr  : " << bodyAdditionAttemptCounter_<< endl;

    return geomModel_->getCopy();;
}
// MODEL SPECIFIC FUNCTIONS==================================================//
//---------------------------------------------------------------------------//
void addModelDistribution::initializeCellZone()
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
	InfoH << addModel_Info << "-- addModelMessage-- "
        << "lambda fraction in controlled region: " << lambdaFraction<< endl;
	return lambdaFraction;
}
//---------------------------------------------------------------------------//
scalar addModelDistribution::returnRandomAngle()
{
    scalar ranNum = 2.0*randGen_.sample01<scalar>() - 1.0;
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
        ranNum = randGen_.sample01<scalar>();
        axisOfRotation[i] = ranNum;
    }

    axisOfRotation /=mag(axisOfRotation);
    return axisOfRotation;
}
//---------------------------------------------------------------------------//
Tuple2<label, scalar> addModelDistribution::returnScaleFactor()
{
    DynamicLabelList  missingParticles;
    scalar maxSizeVolume = 1.0/6.0*3.14*pow(stlBaseSize_
                           *particleSize_[particleSize_.size()-1]
                           *convertToMeters_/stlBaseSize_,3);
    forAll (addedParticlesSize_,size)
    {
        scalar distribDiff = distribution_[size] - 100*addedParticlesSize_[size]/(volumeOfAddedBodies_+SMALL);
        if(distribDiff > 0)
        {
            scalar meanFactor = particleSize_[size]*convertToMeters_/stlBaseSize_;
            scalar meanVolume = 1.0/6.0*3.14*pow(stlBaseSize_ * meanFactor,3);
            missingParticles.append(floor((distribDiff*maxSizeVolume)/meanVolume));
        }
        else
        {
            missingParticles.append(0);
        }
    }

    label totalMissParts(0);
    forAll (missingParticles,size)
    {
        totalMissParts += missingParticles[size];
    }
    label randomMissPart = floor(randGen_.sample01<scalar>()*totalMissParts);
    label missingPart(0);
    forAll (missingParticles,size)
    {
        missingPart = size;
        randomMissPart -= missingParticles[size];
        if(randomMissPart < 0)
        {
            break;
        }
    }
    scalar factor(particleSize_[missingPart - 1] + (particleSize_[missingPart] - particleSize_[missingPart - 1]) * randGen_.sample01<scalar>());
    factor *= convertToMeters_/stlBaseSize_;

    Tuple2<label, scalar> returnValue(missingPart, factor);

    return returnValue;
}
