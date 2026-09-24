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

    openHFDIB-DEM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (Version 3) as published
    by the Free Software Foundation.

    openHFDIB-DEM is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with openHFDIB-DEM. If not, see <http://www.gnu.org/licenses/>.

InNamespace
    Foam

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/
#include "openHFDIBDEM.H"
#include "polyMesh.H"
#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"

#include "interpolationCellPoint.H"
#include "interpolationCell.H"

#include "scalarMatrices.H"
#include "OFstream.H"
#include <iostream>
#include "defineExternVars.H"
#include "parameters.H"
#include "demTimeStepInfo.H"

#define ORDER 2

using namespace Foam;
using namespace contactModel;

//---------------------------------------------------------------------------//
openHFDIBDEM::openHFDIBDEM(const Foam::fvMesh& mesh)
:
mesh_(mesh),
HFDIBDEMDict_
(
    IOobject
    (
        "HFDIBDEMDict",
        "constant",
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
transportProperties_
(
    IOobject
    (
        "transportProperties",
        "constant",
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
bodyNames_(HFDIBDEMDict_.lookup("bodyNames")),
prtcInfoTable_(0),
stepDEM_(readScalar(HFDIBDEMDict_.lookup("stepDEM"))),
adaptiveDEM_(HFDIBDEMDict_.lookupOrDefault<bool>("adaptiveStepDEM", false)),
sweepSafetyTrans_(1.5),
sweepSafetyRot_(3.0),
useDEMdtEstimator_(true),
dtEstimatorVelocityAware_(true),
dtAreaCoeff_(1.0),
dtTangentialFactor_(1.0),
dtRotationalFactor_(1.0),
dtMinDEMdt_(0.0),
dtMaxDEMSubCycles_(0),
dtReportInterval_(1),
deltaTDEM_(GREAT),
dtDiagHertz_(GREAT),
dtDiagRayleigh_(GREAT),
recordSimulation_(readBool(HFDIBDEMDict_.lookup("recordSimulation")))
{
    materialProperties::matProps_insert(
        "None",
        materialInfo("None", 1, 1, 1, 1, 1)
    );

    if(HFDIBDEMDict_.found("recordFirstTimeStep"))
    {
        recordFirstTimeStep_ = readBool(HFDIBDEMDict_.lookup("recordFirstTimeStep"));
    }

    if(HFDIBDEMDict_.found("nSolidsInDomain"))
    {
        solverInfo::setNSolidsThreshold(readLabel(HFDIBDEMDict_.lookup("nSolidsInDomain")));
    }

    //--- adaptiveStepDEM settings. top-level sweepSafety* in an old
    //    case is honored with a deprecation warning
    if (HFDIBDEMDict_.found("adaptiveStepDEMDict"))
    {
        dictionary asdDic(HFDIBDEMDict_.subDict("adaptiveStepDEMDict"));

        useDEMdtEstimator_ =
            asdDic.lookupOrDefault<bool>("useDEMdtEstimator", true);

        word dtEstimator(asdDic.lookupOrDefault<word>("DEMdtEstimator", "velocityAware"));
        if (dtEstimator == "velocityAware")
        {
            dtEstimatorVelocityAware_ = true;
        }
        else if (dtEstimator == "conservative")
        {
            dtEstimatorVelocityAware_ = false;
        }
        else
        {
            FatalIOErrorInFunction(HFDIBDEMDict_)
                << "unknown DEMdtEstimator " << dtEstimator
                << ", expected conservative or velocityAware"
                << exit(FatalIOError);
        }

        sweepSafetyTrans_ =
            asdDic.lookupOrDefault<scalar>("sweepSafetyTrans", 1.5);
        sweepSafetyRot_ =
            asdDic.lookupOrDefault<scalar>("sweepSafetyRot", 3.0);

        dtAreaCoeff_ = asdDic.lookupOrDefault<scalar>("areaCoeff", 1.0);
        dtTangentialFactor_ =
            asdDic.lookupOrDefault<scalar>("tangentialFactor", 1.0);
        dtRotationalFactor_ =
            asdDic.lookupOrDefault<scalar>("rotationalFactor", 1.0);
        dtMinDEMdt_ = asdDic.lookupOrDefault<scalar>("minDEMdt", 0.0);
        dtMaxDEMSubCycles_ = asdDic.lookupOrDefault<label>("maxDEMSubCycles", 0);
        dtReportInterval_ = asdDic.lookupOrDefault<label>("dtReportInterval", 1);
    }
    else
    {
        if (HFDIBDEMDict_.found("sweepSafetyTrans"))
        {
            sweepSafetyTrans_ =
                readScalar(HFDIBDEMDict_.lookup("sweepSafetyTrans"));

            WarningInFunction
                << "top-level sweepSafetyTrans is deprecated, move it "
                << "into the adaptiveStepDEMDict sub-dictionary"
                << endl;
        }

        if (HFDIBDEMDict_.found("sweepSafetyRot"))
        {
            sweepSafetyRot_ =
                readScalar(HFDIBDEMDict_.lookup("sweepSafetyRot"));

            WarningInFunction
                << "top-level sweepSafetyRot is deprecated, move it "
                << "into the adaptiveStepDEMDict sub-dictionary"
                << endl;
        }
    }

    dictionary demDic = HFDIBDEMDict_.subDict("DEM");
    dictionary materialsDic = demDic.subDict("materials");
    List<word> materialsNames = materialsDic.toc();
    if(demDic.found("increasedDamping"))
    {
        contactModelInfo::setIncreasedDamping(readBool(demDic.lookup("increasedDamping")));
    }
    else
    {
        contactModelInfo::setIncreasedDamping(false);
    }

    forAll(materialsNames, matI)
    {
        dictionary matIDic = materialsDic.subDict(materialsNames[matI]);
        scalar eps = readScalar(matIDic.lookup("eps"));
        if(!contactModelInfo::getIncreasedDamping())
        {
            eps = 0.906463027*eps + 0.093538298;                        // LinearRegression based on LIGGGHTS data testing
            eps = min(eps, 1.0);
        }

        materialProperties::matProps_insert(
            materialsNames[matI],
            materialInfo(
                materialsNames[matI],
                readScalar(matIDic.lookup("Y")),
                readScalar(matIDic.lookup("nu")),
                readScalar(matIDic.lookup("mu")),
                readScalar(matIDic.lookup("adhN")),
                // readScalar(matIDic.lookup("eps"))
                eps
            )
        );
    }

    if(demDic.found("interfaceAdh"))
    {
        dictionary interfAdhDic = demDic.subDict("interfaceAdh");
        List<word> interNames = interfAdhDic.toc();
        forAll(interNames, interI)
        {
            dictionary interDicI = interfAdhDic.subDict(interNames[interI]);
            wordList interMat = interDicI.lookup("materials");
            string interKey;
            if(interMat[0] < interMat[1])
            {
                interKey += interMat[0];
                interKey += "-";
                interKey += interMat[1];
            }
            else
            {
                interKey += interMat[1];
                interKey += "-";
                interKey += interMat[0];
            }

            interAdhesion::interAdhesion_insert(
                interKey,
                readScalar(interDicI.lookup("value"))
            );
        }
    }

    if(demDic.found("LcCoeff"))
    {
        contactModelInfo::setLcCoeff(readScalar(demDic.lookup("LcCoeff")));
    }
    else
    {
        contactModelInfo::setLcCoeff(4.0);
    }

    if(demDic.found("rotationModel"))
    {

        word rotModel = word(demDic.lookup("rotationModel"));
        if(rotModel == "chen2012")
        {
            contactModelInfo::setRotationModel(0);
        }
        else if(rotModel == "mindlin1953")
        {
            contactModelInfo::setRotationModel(1);
        }
        else
        {
            Info << "Rotation Model not recognized, setting to default mindlin1953" << endl;
            contactModelInfo::setRotationModel(1);
        }
    }
    else
    {
        Info << "Rotation Model not recognized, setting to default mindlin1953" << endl;
        contactModelInfo::setRotationModel(1);
    }


    Info <<" -- Coefficient for characteristic Lenght Lc is set to : "<< contactModelInfo::getLcCoeff() << endl;

    dictionary patchDic = demDic.subDict("collisionPatches");
    List<word> patchNames = patchDic.toc();
    forAll(patchNames, patchI)
    {
        word patchMaterial = word(patchDic.subDict(patchNames[patchI]).lookup("material"));
        vector patchNVec   = vector(patchDic.subDict(patchNames[patchI]).lookup("nVec"));
        vector planePoint  = vector(patchDic.subDict(patchNames[patchI]).lookup("planePoint"));

        wallPlaneInfo::wallPlaneInfo_insert(
            patchNames[patchI],
            patchNVec,
            planePoint
        );

        wallMatInfo::wallMatInfo_insert(
            patchNames[patchI],
            materialProperties::getMatProps()[patchMaterial]
        );
    }

    if(demDic.found("cyclicPatches"))
    {
        Info << "CyclicPatches Found " << endl;
        dictionary cyclicPatchDic = demDic.subDict("cyclicPatches");
        List<word> cyclicPatchNames = cyclicPatchDic.toc();
        forAll(cyclicPatchNames, patchI)
        {
            vector patchNVec    = vector(cyclicPatchDic.subDict(cyclicPatchNames[patchI]).lookup("nVec"));
            vector planePoint   = vector(cyclicPatchDic.subDict(cyclicPatchNames[patchI]).lookup("planePoint"));
            word neighbourPatch = word(cyclicPatchDic.subDict(cyclicPatchNames[patchI]).lookup("neighbourPatch"));

            cyclicPlaneInfo::insert(
                cyclicPatchNames[patchI],
                patchNVec,
                planePoint,
                neighbourPatch
            );
        }
        Info << "CyclicPatches  " <<  cyclicPatchNames <<endl;
    }

    if (HFDIBDEMDict_.found("geometricD"))
    {
        geometricD = vector(HFDIBDEMDict_.lookup("geometricD"));
    }
    else
    {
        geometricD = mesh_.geometricD();
    }

    forAll (geometricD, direction)
    {
        if (geometricD[direction] == -1)
        {
            case3D = false;
            emptyDir[direction] = 1;
            emptyDim = direction;
            break;
        }
    }

    if (HFDIBDEMDict_.isDict("virtualMesh"))
    {
        dictionary vMDic = HFDIBDEMDict_.subDict("virtualMesh");
        virtualMeshLevel::setVirtualMeshLevel(readScalar(vMDic.lookup("level")),readScalar(vMDic.lookup("charCellSize")));
        virtualMeshLevel::setMaxSubVolumes(vMDic.lookupOrDefault<scalar>("maxSubVolumes",virtualMeshLevel::getMaxSubVolumes()));
        Info <<" -- VirtMesh Decomposition Level is set to        : "<< virtualMeshLevel::getVirtualMeshLevel() << endl;
        Info <<" -- VirtMesh charCellSize for boundary is set to  : "<< virtualMeshLevel::getCharCellSize() << endl;
        Info <<" -- VirtMesh max number of subVolumes is set to   : "<< virtualMeshLevel::getMaxSubVolumes() << endl;

    }
    else
    {
        virtualMeshLevel::setVirtualMeshLevel(1,1);
        Info <<" -- VirtMesh Decomposition Level is set to        : "<< virtualMeshLevel::getVirtualMeshLevel() << endl;
        Info <<" -- VirtMesh charCellSize for boundary is set to  : "<< virtualMeshLevel::getCharCellSize() << endl;
        Info <<" -- VirtMesh max number of subVolumes is set to   : "<< virtualMeshLevel::getMaxSubVolumes() << endl;

    }

    recordOutDir_ = mesh_.time().rootPath() + "/" + mesh_.time().globalCaseName() + "/bodiesInfo";
}
//---------------------------------------------------------------------------//
openHFDIBDEM::~openHFDIBDEM()
{}
//---------------------------------------------------------------------------//
void openHFDIBDEM::initialize
(
    volScalarField& body,
    volVectorField& U,
    volScalarField& refineF,
    label recomputeM0,
    word runTime
)
{
    if(HFDIBDEMDict_.found("outputSetup"))
    {
        dictionary outputDic = HFDIBDEMDict_.subDict("outputSetup");
        bool basicOutput = readBool(outputDic.lookup("basic"));
        bool iBoutput = readBool(outputDic.lookup("iB"));
        bool DEMoutput = readBool(outputDic.lookup("DEM"));
        bool addModelOutput = readBool(outputDic.lookup("addModel"));
        bool parallelDEMOutput = readBool(outputDic.lookup("parallelDEM"));
        InfoH.setOutput(
            basicOutput,
            iBoutput,
            DEMoutput,
            addModelOutput,
            parallelDEMOutput
        );
    }

    preCalculateCellPoints();

    if(HFDIBDEMDict_.found("interpolationSchemes"))
    {
        HFDIBinterpDict_ = HFDIBDEMDict_.subDict("interpolationSchemes");

        if(HFDIBinterpDict_.found("method"))
        {
            word intMethod = word(HFDIBinterpDict_.lookup("method"));

            if(intMethod == "leastSquares")
            {
                dictionary lsCoeffsDict
                    = HFDIBinterpDict_.subDict("leastSquaresCoeffs");

                ibInterp_.reset(new leastSquaresInt(
                    mesh_,
                    readScalar(lsCoeffsDict.lookup("distFactor")),
                    readScalar(lsCoeffsDict.lookup("radiusFactor")),
                    readScalar(lsCoeffsDict.lookup("angleFactor")),
                    readScalar(lsCoeffsDict.lookup("maxCCRows"))
                ));
            }
            else if(intMethod == "line")
            {
                ibInterp_.reset(new lineInt(HFDIBinterpDict_));
            }
        }

        else
        {
            InfoH << basic_Info << "WARN: No interpolation"
                << " method specified, using line as default" << endl;
            ibInterp_.reset(new lineInt(HFDIBinterpDict_));
        }
    }

    else
    {
        InfoH << basic_Info << "WARN: Dictionary"
            << " interpolationSchemes not found, interpolation"
            << " at IB not initialized" << endl;
    }

    bool startTime0(runTime == "0");

    // initialize addModels
    addModels_.setSize(bodyNames_.size());
    immersedBodies_.setSize(0);                                         //on the fly creation

    refineF *= 0.0;
    recomputeM0_ = recomputeM0;

    if(!startTime0)
    {
        if(!isDir(recordOutDir_))
            mkDir(recordOutDir_);
        else
        {
            fileNameList entries(readDir(recordOutDir_,fileName::DIRECTORY));
            scalar runTimeS(stod(runTime));
            forAll(entries,entry)
            {
                scalar dirTime(stod(entries[entry].name()));
                if(dirTime > runTimeS)
                {
                    word pathI(recordOutDir_ + "/" + entries[entry]);
                    rmDir(pathI);
                }
            }
        }

        restartSimulation(body, refineF, runTime);
    }
    else
    {
        if(!isDir(recordOutDir_))
            mkDir(recordOutDir_);
        else
        {
            rmDir(recordOutDir_);
            mkDir(recordOutDir_);
        }
    }

    #include "initializeAddModels.H"

    forAll (addModels_,modelI)
    {
        word bodyName(bodyNames_[modelI]);
        InfoH << basic_Info << "Creating immersed body based on: " << bodyName << endl;

        label maxAdditions(1000);
        label cAddition(0);

        // Note (MI): add solid body only if the number of immersed 
        //            bodies is below the threshold (if set)
        // Note (MI): the second condition is required to avoid skipping
        //            addition of all the bodies on first time step
        while (addModels_[modelI].shouldAddBody(body) and cAddition < maxAdditions and (solverInfo::getNSolidsThreshold() < 0 or immersedBodies_.size() < solverInfo::getNSolidsThreshold()))
        {
            InfoH << addModel_Info << "addModel invoked action, trying to add new body" << endl;
            std::shared_ptr<geomModel> bodyGeomModel(addModels_[modelI].addBody(body, immersedBodies_));
            cAddition++;

            // initialize the immersed bodies
            if (addModels_[modelI].getBodyAdded())
            {
                label newIBSize(immersedBodies_.size()+1);
                label addIBPos(newIBSize - 1);
                immersedBodies_.setSize(newIBSize);

                InfoH << addModel_Info << "Trying to set immersedBodies" << endl;
                immersedBodies_.set
                (
                    addIBPos,
                    new immersedBody
                    (
                        bodyName,
                        mesh_,
                        HFDIBDEMDict_,
                        transportProperties_,
                        addIBPos,
                        recomputeM0_,
                        bodyGeomModel,
                        ibInterp_,
                        cellPoints_
                    )
                );
                immersedBodies_[addIBPos].createImmersedBody(body,refineF);
                immersedBodies_[addIBPos].computeBodyCharPars();
                if (immersedBodies_[addIBPos].getStartSynced())
                {
                    immersedBodies_[addIBPos].initSyncWithFlow(U);
                }
                verletList_.addBodyToVList(immersedBodies_[addIBPos]);
                InfoH << addModel_Info << "Body based on: " << bodyName << " successfully added" << endl;
                InfoH << addModel_Info << "Current count of solids within the domain : " << immersedBodies_.size() << endl;
                cAddition = 0;
            }
            else
            {
                InfoH << addModel_Info << "Body based on: "
                    << bodyName << " should have been added but was not "
                    << "(probably overlap with an already existing body)"
                    << endl;
            }
        }
    }

    verletList_.initialSorting();

    body.correctBoundaryConditions();
    refineF.correctBoundaryConditions();
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::createBodies(volScalarField& body,volScalarField& refineF)
{
    // add/remove bookkeeping is accumulated over a single time step,
    // createBodies runs first, so this is where it gets reset
    nBodiesRemovedLastStep_ = 0;

    // static-body re-creation skip
    // Note (MI): if required, make this a HFDIBDEMDict entry
    const bool skipStaticRecreation(true);

    // per-body re-creation decision, shared by the reset and the create
    boolList recreateBody(immersedBodies_.size(), true);
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            if (skipStaticRecreation
                && immersedBodies_[bodyId].getbodyOperation() == 0)
            {
                recreateBody[bodyId] = !immersedBodies_[bodyId]
                    .getGeomModel().bodyFieldValid();
            }

            reduce(recreateBody[bodyId], andOp<bool>());
        }
    }

    // reset all bodies before recreation
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive() && recreateBody[bodyId])
        {
            immersedBodies_[bodyId].resetBody(body);
        }
    }

    // recreate all bodies after contact update
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            if (recreateBody[bodyId])
            {
                immersedBodies_[bodyId].createImmersedBody
                (
                    body,
                    refineF,
                    false                                               //create without "synchronization"
                );
            }
            else
            {
                // skipped static body: refresh intpInfo from
                // the (unchanged) cached surfCells here
                immersedBodies_[bodyId].refreshIntpInfo();
            }
        }
    }

    createBodiesComputeDynamicsVars(body, refineF);

    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].checkIfInDomain(body);
            if (!immersedBodies_[bodyId].getIsActive())
            {
                nBodiesRemovedLastStep_++;
                removeBodyContacts(bodyId);
            }
            immersedBodies_[bodyId].updateOldMovementVars();
            immersedBodies_[bodyId].checkBodyOp();
        }
    }

    // make the counter rank-uniform
    reduce(nBodiesRemovedLastStep_, maxOp<label>());
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::recreateBodies
(
    volScalarField& body,
    volScalarField& refineF
)
{
    refineF *= 0.0;
    preCalculateCellPoints();
    forAll (addModels_,modelI)
    {
        addModels_[modelI].recreateBoundBox();
    }
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].recreateBodyField(body,refineF);
        }
    }

    createBodiesComputeDynamicsVars(body, refineF);

    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].checkIfInDomain(body);
            if (!immersedBodies_[bodyId].getIsActive())
            {
                nBodiesRemovedLastStep_++;
                removeBodyContacts(bodyId);
            }
            if(immersedBodies_[bodyId].getRecomputeM0() > 0)
            {
                immersedBodies_[bodyId].computeBodyCharPars();
                immersedBodies_[bodyId].recomputedM0();
            }
            InfoH << iB_Info << "-- body "
                << immersedBodies_[bodyId].getBodyIdStr() << " re-created" << endl;
        }
    }

    // make the counter rank-uniform
    reduce(nBodiesRemovedLastStep_, maxOp<label>());
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::createBodiesComputeDynamicsVars(
    volScalarField& body,
    volScalarField& refineF
)
{
    DynamicList<scalar> particleMasses;
    DynamicList<label> particleCells;
    DynamicList<symmTensor> particleInertiaTensors;

    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].syncImmersedBodyGeometry(body,refineF);
            if (immersedBodies_[bodyId].getGeomModel().isCluster())
            {
                clusterBody& cBody = dynamic_cast<clusterBody&>(immersedBodies_[bodyId].getGeomModel());
                std::vector<std::shared_ptr<geomModel>>& cBodies = cBody.getClusterBodies();
                for (auto& cB : cBodies)
                {
                    particleMasses.append(cB->getM());
                    particleCells.append(cB->getNCells());
                    particleInertiaTensors.append(cB->getI());
                }
            }
            else
            {
                particleMasses.append(immersedBodies_[bodyId].getGeomModel().getM());
                particleCells.append(immersedBodies_[bodyId].getGeomModel().getNCells());
                particleInertiaTensors.append(immersedBodies_[bodyId].getGeomModel().getI());
            }
        }
    }
    reduce(particleMasses,sumOp<List<scalar>>());
    reduce(particleCells,sumOp<List<label>>());
    reduce(particleInertiaTensors,sumOp<List<symmTensor>>());

    label bodyIndex(0);
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            if (immersedBodies_[bodyId].getGeomModel().isCluster())
            {
                clusterBody& cBody = dynamic_cast<clusterBody&>(immersedBodies_[bodyId].getGeomModel());
                std::vector<std::shared_ptr<geomModel>>& cBodies = cBody.getClusterBodies();
                for (auto& cB : cBodies)
                {
                    cB->setM(particleMasses[bodyIndex]);
                    cB->setNCells(particleCells[bodyIndex]);
                    cB->setI(particleInertiaTensors[bodyIndex]);
                    bodyIndex++;
                }
                cBody.setMassAndInertia();
            }
            else
            {
                immersedBodies_[bodyId].getGeomModel().setM(particleMasses[bodyIndex]);
                immersedBodies_[bodyId].getGeomModel().setNCells(particleCells[bodyIndex]);
                immersedBodies_[bodyId].getGeomModel().setI(particleInertiaTensors[bodyIndex]);
                bodyIndex++;
            }

            immersedBodies_[bodyId].syncImmersedBodyRefinement(body,refineF);
        }
    }

    body.correctBoundaryConditions();
    refineF.correctBoundaryConditions();
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::updateSurface
(
    const scalar surfaceThreshold,
    const volScalarField& body,
    volScalarField& surface
)
{
    forAll(surface, sI)
    {
        if (body[sI] > surfaceThreshold)
            surface[sI] = 1;
        else
            surface[sI] = 0;
    }
    surface.correctBoundaryConditions();
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::updateSurface
(
    const scalar surfaceThreshold,
    const volScalarField& body,
    volScalarField& surface,
    surfaceScalarField& surfaceF
)
{
    updateSurface(surfaceThreshold, body, surface);
    surfaceF = fvc::interpolate(surface);
    surfaceF.correctBoundaryConditions();
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::updateGlobalFluidDensity(
    const volScalarField& body,
    volScalarField& rho
)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].updateLocalFluidDensity(body,rho);
        }
    }
    rho.correctBoundaryConditions();
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::updateGlobalFluidDensity(
    const volScalarField& body,
    volScalarField& rho,
    volScalarField& rhoS
)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].updateLocalFluidDensity(body,rho,rhoS);
        }
    }
    rho.correctBoundaryConditions();
    rhoS.correctBoundaryConditions();
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::preUpdateBodies
(
    volScalarField& body
)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            // create body or compute body-fluid coupling and estimate
            // potential contacts with walls
            immersedBodies_[bodyId].inContactWithStatic(false);

            immersedBodies_[bodyId].updateOldMovementVars();
            immersedBodies_[bodyId].printStats();
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::postUpdateBodies
(
    const volScalarField& body,
    const volVectorField& f
)
{
    postUpdateBodies(body, f, false, false);
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::postUpdateBodies
(
    const volScalarField& body,
    const volVectorField& f,
    const bool kinematicForce,
    const bool applyAddedMass
)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].postPimpleUpdateImmersedBody
            (
                body,
                f,
                kinematicForce,
                applyAddedMass
            );
            immersedBodies_[bodyId].clearIntpInfo();
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::postUpdateBodies
(
    const volScalarField& body,
    const volVectorField& f,
    const volScalarField& rho,
    const bool kinematicForce,
    const bool applyAddedMass
)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].postPimpleUpdateImmersedBody
            (
                body,
                f,
                rho,
                kinematicForce,
                applyAddedMass
            );
            immersedBodies_[bodyId].clearIntpInfo();
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::interpolateIB
( 
    volVectorField & V,
    volVectorField & Vs,
    volScalarField & body
)
{
    if(ibInterp_.valid())
    {
        ibInterp_->resetInterpolator(V);
    }
    // reset imposed field
    Vs = V;

    // loop over all the immersed bodies
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            // update imposed field according to body
            immersedBodies_[bodyId].updateVectorField(Vs, V.name(),body);

            if(ibInterp_.valid())
            {
                ibInterp_->ibInterpolate
                (
                    immersedBodies_[bodyId].getIntpInfo(),
                    Vs,
                    immersedBodies_[bodyId].getUatIbPoints(),
                    mesh_
                );
            }
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::writeBodiesInfo()
{
    if(!recordSimulation_)
        return;

    word curOutDir(recordOutDir_ + "/" + mesh_.time().timeName());


    mkDir(curOutDir);
    mkDir(curOutDir +"/stlFiles");
    DynamicLabelList activeIB;
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            activeIB.append(bodyId);
        }
    }
    wordList bodyNames;
    scalar listZize(activeIB.size());
    label bodiesPerProc = ceil(listZize/Pstream::nProcs());
    InfoH << basic_Info << "Active IB listZize      : " << listZize<< endl;
    InfoH << basic_Info << "bodiesPerProc : " << bodiesPerProc<< endl;
    // Pout << "Processor "<< Pstream::myProcNo() << endl;

    for(int assignProc = Pstream::myProcNo()*bodiesPerProc; assignProc < min((Pstream::myProcNo()+1)*bodiesPerProc,activeIB.size()); assignProc++)
    {
        const label bodyId(activeIB[assignProc]);
        // Pout <<"Processor "<< Pstream::myProcNo() << " writes Body " << bodyId << endl;
        word path(curOutDir + "/body" + std::to_string(immersedBodies_[bodyId].getBodyId()) +".info");
        OFstream ofStream(path);
        IOobject outClass
            (
                path,
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            );
        IOdictionary outDict(outClass);

        outDict.writeHeader(ofStream);
        immersedBodies_[bodyId].recordBodyInfo(outDict,curOutDir);
        outDict.writeData(ofStream);
    }

}
//---------------------------------------------------------------------------//
void openHFDIBDEM::updateDEM(volScalarField& body,volScalarField& refineF)
{
    if (cyclicPlaneInfo::getCyclicPlaneInfo().size() > 0)
    {
        forAll (immersedBodies_,bodyId)
        {
            if (!immersedBodies_[bodyId].getGeomModel().isCluster())
            {
                vector transVec = vector::zero;

                if (detectCyclicContact(
                    immersedBodies_[bodyId].getWallCntInfo(),
                    transVec
                ))
                {
                    verletList_.removeBodyFromVList(immersedBodies_[bodyId]);

                    std::shared_ptr<periodicBody> newPeriodicBody
                        = std::make_shared<periodicBody>(mesh_);

                    newPeriodicBody->setRhoS(immersedBodies_[bodyId].getGeomModel().getRhoS());
                    std::shared_ptr<geomModel> iBcopy(immersedBodies_[bodyId].getGeomModel().getCopy());
                    iBcopy->bodyMovePoints(transVec);
                    newPeriodicBody->addBodyToCluster(immersedBodies_[bodyId].getGeomModelPtr());
                    newPeriodicBody->addBodyToCluster(iBcopy);
                    immersedBodies_[bodyId].getGeomModelPtr() = newPeriodicBody;

                    verletList_.addBodyToVList(immersedBodies_[bodyId]);
                    Info << "Periodic body created for body " << bodyId << endl;
                }
            }
            else
            {
                periodicBody& cBody = dynamic_cast<periodicBody&>(immersedBodies_[bodyId].getGeomModel());

                if(cBody.shouldBeUnclustered())
                {
                    verletList_.removeBodyFromVList(immersedBodies_[bodyId]);

                    immersedBodies_[bodyId].getGeomModelPtr() = cBody.getRemGeomModel();

                    verletList_.addBodyToVList(immersedBodies_[bodyId]);
                    Info << "Periodic body unclustered for body " << bodyId << endl;
                }
            }
        }
    }

    scalar deltaTime(mesh_.time().deltaT().value());
    scalar pos(0.0);

    // adaptive DEM stepping (adaptiveStepDEM): bodies are sorted into the
    // sub-cycled set A (contact potential during this CFD time step) and
    // the single-stepped set B (contact-free).
    HashSet<label,Hash<label>> subCycleSet;
    HashSet<label,Hash<label>> singleStepSet;

    // classification internals reused by the dt estimator
    HashSet<label,Hash<label>> hardFlags;
    HashTable<scalar,label,Hash<label>> sweepDist;

    if (adaptiveDEM_)
    {
        computeAdaptiveSets(subCycleSet, singleStepSet, hardFlags, sweepDist);
    }
    // DEM step selection: with the estimator on, deltaTDEM comes from
    // the pair stability steps (floor, caps and the stepDEM safety
    // factor applied in computeDEMdtEstimate); otherwise the legacy
    // stepDEM fraction of the CFD step (1.0 for an empty set)
    scalar step
    (
        adaptiveDEM_ && subCycleSet.empty()
        ? 1.0
        : stepDEM_
    );

    if (adaptiveDEM_)
    {
        if (useDEMdtEstimator_ && !subCycleSet.empty())
        {
            computeDEMdtEstimate(subCycleSet, hardFlags, sweepDist);

            step = min(deltaTDEM_/deltaTime, 1.0);

            // a zero estimate means stale contact state (a
            // switched-off body feeding the estimator); integrating
            // would need 2^31 subcycles, so fall back to the
            // conservative fixed step instead
            if (step < SMALL)
            {
                WarningInFunction
                    << "zero DEM stability estimate - stale contact"
                    << " state suspected (body switch-off?),"
                    << " falling back to a single DEM step"
                    << endl;
                step = 1.0;
            }

            // report the subcycle count and apply the cap
            label nSub(ceil(1.0/max(step, SMALL)));

            if (dtMaxDEMSubCycles_ > 0 && nSub > dtMaxDEMSubCycles_)
            {
                nSub = dtMaxDEMSubCycles_;
                step = 1.0/scalar(nSub);

                WarningInFunction
                    << "DEM subcycles capped at maxDEMSubCycles "
                    << dtMaxDEMSubCycles_
                    << ": the estimated stability bound is being "
                    << "exceeded"
                    << endl;
            }

            InfoH << DEM_Info << " adaptiveStepDEM: "
                << subCycleSet.size() << " sub-cycled / "
                << singleStepSet.size() << " single-stepped bodies"
                << ", DEM subcycles: " << nSub
                << ", DEM step: " << step << endl;
        }
        else
        {
            InfoH << DEM_Info << " adaptiveStepDEM: "
                << subCycleSet.size() << " sub-cycled / "
                << singleStepSet.size() << " single-stepped bodies"
                << ", DEM step: " << step << endl;
        }
    }
    // scalar timeStep(step*deltaTime);
    List<DynamicList<pointField>> bodiesPositionList(Pstream::nProcs());
    // Infos <<bodiesPositionList.size() << endl;
    HashTable <label,Tuple2<label, label>,Hash<Tuple2<label, label>>> syncOutForceKeyTable;
    HashTable <label,Tuple2<label, label>,Hash<Tuple2<label, label>>> contactResolvedKeyTable;
    HashTable <label,label,Hash<label>> wallContactIBTable;
    while( pos < 1)
    {
        bodiesPositionList[Pstream::myProcNo()].clear();

        InfoH << DEM_Info << " Start DEM pos: " << pos
            << " DEM step: " << step << endl;

        InfoH << basic_Info << " DEM - CFD Time: "
            << mesh_.time().value() + deltaTime*pos << endl;

        // single-stepped bodies move only in the first loop
        // iteration (pos == 0, stepDEM = 1)
        const bool firstIter(pos < SMALL);

        auto ibSingleStep = [&](label ib) -> bool
        {
            return adaptiveDEM_ && singleStepSet.found(ib);
        };

        forAll (immersedBodies_,ib)
        {
            if (ibSingleStep(ib) && !firstIter) continue;

            const scalar ibStep(ibSingleStep(ib) ? 1.0 : step);

            immersedBodies_[ib].updateMovement(deltaTime*ibStep*0.5);

            // rotate the cached inertia on all ranks with the same
            // sub-step rotation moveImmersedBody applies on rank 0 -
            // keeps I_ rank-uniform for the next updateMovement
            // half-step without any communication
            immersedBodies_[ib].rotateCachedInertia(deltaTime*ibStep);

            if(Pstream::myProcNo() == 0 )
            {
                immersedBodies_[ib].moveImmersedBody(deltaTime*ibStep);
                if(immersedBodies_[ib].getGeomModel().getcType() != cluster)
                {
                    bodiesPositionList[Pstream::myProcNo()].append(immersedBodies_[ib].getGeomModel().getBodyPoints());
                }
                else
                {
                    clusterBody& cBody = dynamic_cast<clusterBody&>(immersedBodies_[ib].getGeomModel());
                    std::vector<std::shared_ptr<geomModel>>& cBodies = cBody.getClusterBodies();
                    for (auto& cB : cBodies)
                    {
                        bodiesPositionList[Pstream::myProcNo()].append(cB->getBodyPoints());
                    }
                }
            }
        }

        Pstream::gatherList(bodiesPositionList,0);
        Pstream::scatterList(bodiesPositionList,0);

        label bodyIndex(0);
        forAll (immersedBodies_,ib)
        {
            if (ibSingleStep(ib) && !firstIter) continue;

            if(immersedBodies_[ib].getGeomModel().getcType() != cluster)
            {
                immersedBodies_[ib].getGeomModel().setBodyPosition(bodiesPositionList[0][bodyIndex++]);
            }
            else
            {
                clusterBody& cBody = dynamic_cast<clusterBody&>(immersedBodies_[ib].getGeomModel());
                std::vector<std::shared_ptr<geomModel>>& cBodies = cBody.getClusterBodies();
                for (auto& cB : cBodies)
                {
                    cB->setBodyPosition(bodiesPositionList[0][bodyIndex++]);
                }
            }
        }

        bodiesPositionList[Pstream::myProcNo()].clear();

        verletList_.update(immersedBodies_);

        DynamicLabelList wallContactIB;
        wallContactIBTable.clear();
        forAll (immersedBodies_,bodyId)
        {
            immersedBody& cIb(immersedBodies_[bodyId]);
            if (cIb.getIsActive())
            {
                // set F_ and T_ to zero.
                cIb.resetContactForces();

                if(cIb.getbodyOperation() != 0)
                {
                    // detect wall contact
                    if(
                        detectWallContact
                        (
                            mesh_,
                            cIb.getibContactClass(),
                            cIb.getWallCntInfo()
                        )
                    )
                    {
                        cIb.getibContactClass().setWallContact(true);
                        cIb.getibContactClass().inContactWithStatic(true);
                        wallContactIB.append(bodyId);
                        wallContactIBTable.insert(bodyId,wallContactIB.size()-1);
                        // cIb.getWallCntInfo().registerSubContactList(wallContactList);
                    }
                }
            }
        }
        // possibleWallContacts = wallContactIB.size();
        List<bool> wallContactResolvedList(wallContactIB.size(),false);

        if(wallContactIB.size() > 0)
        {
            label wallContactPerProc(ceil(double(wallContactIB.size())/Pstream::nProcs()));
            // Info <<" wallContactPerProc : "<< wallContactPerProc << endl;
            if(wallContactIB.size() <= Pstream::nProcs())
            {
                wallContactPerProc = 1;
            }
            for
            (
                int assignProc = Pstream::myProcNo()*wallContactPerProc; 
                assignProc < min((Pstream::myProcNo()+1)*wallContactPerProc,wallContactIB.size()); 
                assignProc++
            )
            {
                immersedBody& cIb(immersedBodies_[wallContactIB[assignProc]]);
                if(cIb.getGeomModel().getcType() != sphere && cIb.getGeomModel().getcType() != cluster)
                {
                    cIb.getWallCntInfo().findContactAreas();
                }

                DynamicList<wallSubContactInfo*> wallContactList;
                cIb.getWallCntInfo().registerSubContactList(wallContactList);
                List<bool> wallcRList(wallContactList.size(),false);

                forAll(wallContactList,sC)
                {
                    wallSubContactInfo* sCW = wallContactList[sC];
                    bool resolved(solveWallContact(
                        mesh_,
                        cIb.getWallCntInfo(),
                        deltaTime*step,
                        *sCW
                        ));
                    sCW->setResolvedContact(resolved);
                    wallContactResolvedList[assignProc] += resolved;
                    wallcRList[sC] = resolved;
                }
            }

            reduce(wallContactResolvedList,sumOp<List<bool>>());

            List<vector> iBodyOutForceList(wallContactIB.size(),vector::zero);
            List<vector> iBodyOutTorqueList(wallContactIB.size(),vector::zero);

            forAll (wallContactIB,iB)
            {
                immersedBody& cIb(immersedBodies_[wallContactIB[iB]]);
                if(wallContactIBTable.found(cIb.getBodyId()))
                {
                    label cKey(wallContactIBTable[cIb.getBodyId()]);
                    if(wallContactResolvedList[cKey])
                    {
                        std::vector<std::shared_ptr<wallSubContactInfo>>& subCList
                            = cIb.getWallCntInfo().getWallSCList();

                        for(auto sCW : subCList)
                        {
                            iBodyOutForceList[cKey] += sCW->getOutForce().F;
                            iBodyOutTorqueList[cKey] += sCW->getOutForce().T;
                        }
                    }
                }
            }
            reduce(iBodyOutForceList,sumOp<List<vector>>());
            reduce(iBodyOutTorqueList,sumOp<List<vector>>());

            forAll (wallContactIB,iB)
            {
                immersedBody& cIb(immersedBodies_[wallContactIB[iB]]);
                forces cF;
                cF.F = iBodyOutForceList[iB];
                cF.T = iBodyOutTorqueList[iB];

                cIb.updateContactForces(cF);
                cIb.getWallCntInfo().clearOldContact();
                
                cIb.resetCouplingHistory();
            }
        }

        wallContactIB.clear();

        DynamicList<prtSubContactInfo*> contactList;
        // check only pairs whose bounding boxes are intersected for the contact
        label vListSize(0);
        for (auto it = verletList_.begin(); it != verletList_.end(); ++it)
        {
            const Tuple2<label, label> cPair = Tuple2<label, label>(it->first, it->second);

            label cInd(cPair.first());
            bool cStatic(immersedBodies_[cInd].getbodyOperation() == 0);

            label tInd(cPair.second());
            bool tStatic(immersedBodies_[tInd].getbodyOperation() == 0);

            if
            (
                (immersedBodies_[cInd].getIsActive() && immersedBodies_[tInd].getIsActive())
                &&
                !(cStatic && tStatic)
            )
            {
                if(cStatic)
                    immersedBodies_[tInd].inContactWithStatic(true);

                if(tStatic)
                    immersedBodies_[cInd].inContactWithStatic(true);

                prtContactInfo& prtcInfo(getPrtcInfo(
                    cPair)
                );

                prtcInfo.clearData();
                getContacts(
                    mesh_,
                    prtcInfo
                );

                // prtcInfo.syncContactList();

                prtcInfo.registerContactList(contactList);
            }
            vListSize++;
        }

        List<bool> contactResolved(contactList.size(),false);
        List<label> contactResolvedcKey(contactList.size(),0);
        List<label> contactResolvedtKey(contactList.size(),0);

        if(contactList.size() > 0 )
        {
            label contactPerProc(ceil(double(contactList.size())/Pstream::nProcs()));
            if(contactList.size() <= Pstream::nProcs())
            {
                contactPerProc = 1;
            }

            for
            (
                int assignProc = Pstream::myProcNo()*contactPerProc; 
                assignProc < min((Pstream::myProcNo()+1)*contactPerProc,contactList.size()); 
                assignProc++
            )
            {
                prtSubContactInfo* sCI = contactList[assignProc];
                const Tuple2<label, label>& cPair = sCI->getCPair();

                contactResolvedcKey[assignProc] = cPair.first();
                contactResolvedtKey[assignProc] = cPair.second();

                ibContactClass& cClass(immersedBodies_[cPair.first()].getibContactClass());
                ibContactClass& tClass(immersedBodies_[cPair.second()].getibContactClass());

                if(detectPrtPrtContact(mesh_,cClass,tClass,*sCI))
                {
                    // resolvedPrtContacts++;
                    prtContactInfo& prtcInfo(getPrtcInfo(cPair));

                    bool resolved(solvePrtContact(mesh_, prtcInfo, *sCI, deltaTime*step));
                    sCI->setResolvedContact(resolved);

                    contactResolved[assignProc] += resolved;
                }
            }
        }

        reduce(contactResolved,sumOp<List<bool>>());
        reduce(contactResolvedcKey,sumOp<List<label>>());
        reduce(contactResolvedtKey,sumOp<List<label>>());
        contactResolvedKeyTable.clear();
        forAll(contactResolvedcKey,cKey)
        {
            contactResolvedKeyTable.insert(Tuple2<label, label>(contactResolvedcKey[cKey],contactResolvedtKey[cKey]),cKey);
        }
        List<vector> cBodyOutForceList(vListSize,vector::zero);
        List<vector> cBodyOutTorqueList(vListSize,vector::zero);
        List<vector> tBodyOutForceList(vListSize,vector::zero);
        List<vector> tBodyOutTorqueList(vListSize,vector::zero);
        syncOutForceKeyTable.clear();

        label nIter(0);
        for (auto it = verletList_.begin(); it != verletList_.end(); ++it)
        {
            const Tuple2<label, label> cPair = Tuple2<label, label>(it->first, it->second);

            prtContactInfo& prtcInfo(getPrtcInfo(cPair));

            if(contactResolvedKeyTable.found(cPair))
            {
                label nSubContact(0);
                std::vector<std::shared_ptr<prtSubContactInfo>>& subCList
                    = prtcInfo.getPrtSCList();
                for(auto sC : subCList)
                {
                    nSubContact++;
                    cBodyOutForceList[nIter] += sC->getOutForce().first().F;
                    cBodyOutTorqueList[nIter] += sC->getOutForce().first().T;
                    tBodyOutForceList[nIter] += sC->getOutForce().second().F;
                    tBodyOutTorqueList[nIter] += sC->getOutForce().second().T;
                }
            }
            syncOutForceKeyTable.insert(cPair,nIter);
            nIter++;
        }

        reduce(cBodyOutForceList,sumOp<List<vector>>());
        reduce(cBodyOutTorqueList,sumOp<List<vector>>());
        reduce(tBodyOutForceList,sumOp<List<vector>>());
        reduce(tBodyOutTorqueList,sumOp<List<vector>>());

        label nvListIter(0);

        for (auto it = verletList_.begin(); it != verletList_.end(); ++it)
        {
            const Tuple2<label, label> cPair = Tuple2<label, label>(it->first, it->second);
            label cInd(cPair.first());
            label tInd(cPair.second());

            if(!contactResolvedKeyTable.found(cPair))
            {
                if(prtcInfoTable_.found(cPair))
                {
                    prtcInfoTable_.erase(cPair);
                    continue;
                }
            }
            else if(contactResolved[contactResolvedKeyTable[cPair]])
            {
                if(!syncOutForceKeyTable.found(cPair))
                {
                    Pout <<" -- cPair  "<<cInd << " - "<<tInd << " not found in syncOutForceKeyTable" << endl;
                    continue;
                }

                nvListIter = syncOutForceKeyTable[cPair];
                if(nvListIter > cBodyOutForceList.size())
                {
                    Pout <<" -- cPair  "<<cInd << " - "<<tInd << " nvListIter > bodiesOutForceList[Pstream::myProcNo()].size()" << endl;
                    continue;
                }

                vector F1 = vector::zero;
                vector T1 = vector::zero;
                vector F2 = vector::zero;
                vector T2 = vector::zero;

                F1 += cBodyOutForceList[nvListIter];
                T1 += cBodyOutTorqueList[nvListIter];
                F2 += tBodyOutForceList[nvListIter];
                T2 += tBodyOutTorqueList[nvListIter];

                forces cF;
                cF.F = F1;
                cF.T = T1;
                forces tF;
                tF.F = F2;
                tF.T = T2;

                immersedBodies_[cInd].updateContactForces(cF);
                immersedBodies_[tInd].updateContactForces(tF);
                
                immersedBodies_[cInd].resetCouplingHistory();
                immersedBodies_[tInd].resetCouplingHistory();
            }
            else
            {
                if(prtcInfoTable_.found(cPair))
                {
                    prtcInfoTable_.erase(cPair);
                    continue;
                }
            }
        }

        forAll (immersedBodies_,ib)
        {
            if (ibSingleStep(ib) && !firstIter) continue;

            const scalar ibStep(ibSingleStep(ib) ? 1.0 : step);

            immersedBodies_[ib].updateMovement(deltaTime*ibStep*0.5);
            immersedBodies_[ib].printBodyInfo();
            // immersedBodies_[ib].computeBodyCoNumber();
            // if (maxCoNum < immersedBodies_[ib].getCoNum())
            // {
                // maxCoNum = immersedBodies_[ib].getCoNum();
                // bodyId = ib;
            // }
        }
        // InfoH << basic_Info << "Max CoNum = " << maxCoNum << " at body " << bodyId << endl;

        pos += step;

        if (pos + step + SMALL >= 1)
            step = 1 - pos;
//OS Time effitiency Testing
        // demItegrationTime_ = DEMIntergrationRun.timeIncrement();
//OS Time effitiency Testing
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::computeAdaptiveSets
(
    HashSet<label,Hash<label>>& subCycle,
    HashSet<label,Hash<label>>& singleStep,
    HashSet<label,Hash<label>>& hardFlags,
    HashTable<scalar,label,Hash<label>>& sweepDist
)
{
    subCycle.clear();
    singleStep.clear();
    hardFlags.clear();
    sweepDist.clear();

    const scalar deltaTime(mesh_.time().deltaT().value());

    const HashTable<List<vector>,string,Hash<string>>& cyclicPlanes(
        cyclicPlaneInfo::getCyclicPlaneInfo());
    const HashTable<List<vector>,string,Hash<string>>& wallPlanes(
        wallPlaneInfo::getWallPlaneInfo());

    forAll (immersedBodies_, ib)
    {
        immersedBody& cIb(immersedBodies_[ib]);

        if (!cIb.getIsActive()) continue;

        // static bodies are neither sub-cycled nor single-stepped - they
        // do not integrate at all, they only act as contact obstacles
        if (cIb.getbodyOperation() == 0) continue;

        bool hardFlag(false);

        // residual contact force from the previous step: acceleration
        // evolves with the contact and cannot be bounded
        if (mag(cIb.getFContact().F) > SMALL
            || mag(cIb.getFContact().T) > SMALL)
        {
            hardFlag = true;
        }

        // clusters/periodic bodies span cyclic planes
        if (cIb.getGeomModel().isCluster())
        {
            hardFlag = true;
        }

        // body currently in wall contact (keeps resolving the contact)
        if (cIb.checkWallContact())
        {
            hardFlag = true;
        }

        const scalar s
        (
            cIb.computeSweepDistance
            (
                deltaTime,
                sweepSafetyTrans_,
                sweepSafetyRot_
            )
        );
        sweepDist.insert(ib, s);

        // inflating the bbox by s for the plane tests below is only
        // meaningful for soft bodies; hard-flagged bodies are always
        // sub-cycled and need no test
        if (hardFlag)
        {
            hardFlags.insert(ib);
            continue;
        }

        //--- wall contact potential: the body stays inside bbox XOR ball(s)
        //    during the step; the swept region reaches plane (n, p0) if
        //    the most outboard bbox corner is within s of the plane.
        //    This generalizes the corner-plane test of
        //    wallContactInfo::detectWallContact to the whole step.
        if (wallPlanes.size() > 0)
        {
            const stringList wallNames(wallPlanes.toc());

            if (cIb.getGeomModel().getcType() == sphere)
            {
                const vector CoM(cIb.getGeomModel().getCoM());
                const scalar rad(cIb.getGeomModel().getDC()/2.0);

                forAll(wallNames, wI)
                {
                    const List<vector>& planeInfo(wallPlanes[wallNames[wI]]);
                    const vector& n(planeInfo[0]);
                    const vector& p0(planeInfo[1]);

                    // signed distance of the outboard sphere point (the
                    // one closest to crossing the plane; nVec points from
                    // the fluid into the wall, so the fluid side is
                    // negative); potential iff it has come within s of the
                    // plane, i.e. the distance has risen above -s
                    if (((CoM - p0) & n) + rad > -s)
                    {
                        hardFlag = true;
                        break;
                    }
                }
            }
            else
            {
                pointField bbPoints(cIb.getGeomModel().getBounds().points());

                forAll(wallNames, wI)
                {
                    const List<vector>& planeInfo(wallPlanes[wallNames[wI]]);
                    const vector& n(planeInfo[0]);
                    const vector& p0(planeInfo[1]);

                    scalar dOut(-GREAT);
                    forAll(bbPoints, bP)
                    {
                        dOut = max(dOut, ((bbPoints[bP] - p0) & n));
                    }

                    if (dOut > -s)
                    {
                        hardFlag = true;
                        break;
                    }
                }
            }

            if (hardFlag)
            {
                hardFlags.insert(ib);
                continue;
            }
        }

        //--- cyclic proximity: the body may meet its own periodic image,
        //    which the pairwise test cannot see
        if (cyclicPlanes.size() > 0)
        {
            pointField bbPoints(cIb.getGeomModel().getBounds().points());
            const stringList cyclicNames(cyclicPlanes.toc());

            forAll(cyclicNames, cI)
            {
                const List<vector>& planeInfo(cyclicPlanes[cyclicNames[cI]]);
                const vector& n(planeInfo[0]);
                const vector& p0(planeInfo[1]);

                scalar dOut(-GREAT);
                forAll(bbPoints, bP)
                {
                    dOut = max(dOut, ((bbPoints[bP] - p0) & n));
                }

                if (dOut > -s)
                {
                    hardFlag = true;
                    break;
                }
            }

            if (hardFlag)
            {
                hardFlags.insert(ib);
                continue;
            }
        }
    }

    //--- pairwise potential: both members of every sweep-inflated pair
    //    are sub-cycled
    HashSet<label,Hash<label>> potentialBodies;
    verletList_.computePotentialBodies
    (
        immersedBodies_,
        sweepDist,
        hardFlags,
        potentialBodies
    );

    // hard flags go to the sub-cycled set directly
    {
        const labelList hardFlagList(hardFlags.toc());
        forAll(hardFlagList, hI)
        {
            subCycle.insert(hardFlagList[hI]);
        }
    }

    // pair members (static members included in the verlet output are
    // filtered here since they never integrate)
    {
        const labelList potentialList(potentialBodies.toc());
        forAll(potentialList, pI)
        {
            const label ib(potentialList[pI]);

            if (immersedBodies_[ib].getbodyOperation() == 0) continue;

            subCycle.insert(ib);
        }
    }

    // remaining active, non-static bodies: no contact potential - one
    // DEM step per CFD step
    forAll (immersedBodies_, ib)
    {
        immersedBody& cIb(immersedBodies_[ib]);

        if (!cIb.getIsActive()) continue;
        if (cIb.getbodyOperation() == 0) continue;
        if (subCycle.found(ib)) continue;

        singleStep.insert(ib);
    }

    // bodies are replicated across ranks but their force state is only
    // as synchronized as postUpdateBodies makes it
    {
        List<label> setFlag(immersedBodies_.size(), 0);

        forAll(immersedBodies_, ib)
        {
            if (subCycle.found(ib)) setFlag[ib] = 2;
            else if (singleStep.found(ib)) setFlag[ib] = 1;
        }

        reduce(setFlag, maxOp<List<label>>());

        subCycle.clear();
        singleStep.clear();

        forAll(setFlag, ib)
        {
            if (setFlag[ib] == 2) subCycle.insert(ib);
            else if (setFlag[ib] == 1) singleStep.insert(ib);
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::computeDEMdtEstimate
(
    const HashSet<label,Hash<label>>& subCycle,
    const HashSet<label,Hash<label>>& hardFlags,
    const HashTable<scalar,label,Hash<label>>& sweepDist
)
{
    deltaTDEM_ = GREAT;
    dtDiagHertz_ = GREAT;
    dtDiagRayleigh_ = GREAT;

    // governing pair bookkeeping for the log (argmin with the
    // reduce); govIsContact: 0 = contact, 1 = pre-contact, -1 = none;
    // govIsRot: 1 = the rotational bound of the governing pair
    // governs, 0 otherwise
    label govC(-1);
    label govT(-1);
    label govIsContact(-1);
    label govIsRot(0);

    //--- the estimate runs over two pair classes:
    //    contact - pairs with an ongoing contact (prtcInfoTable_):
    //    sized from the cached contact state of the previous substep
    //    (contactVolume_, contactArea_, Lc_, reduced moduli and mass)
    //    pre-contact - sweep-inflated pairs that may first touch
    //    within this CFD step: no contact object exists yet, so the
    //    contact geometry is estimated (exact pair materials, masses
    //    and charCellSize; first-touch area and Lc from h). hard-flagged
    //    members are excluded from the inflated set (their GREAT sweep
    //    would pair them with everything)
    {
        for (auto it = verletList_.begin(); it != verletList_.end(); ++it)
        {
            const Tuple2<label, label> cPair
            (
                Tuple2<label, label>(it->first, it->second)
            );

            // defense in depth: a switched-off body's table entries
            // are purged at switch-off, but never size from them
            // even if one slips through
            if (!immersedBodies_[cPair.first()].getIsActive()
                || !immersedBodies_[cPair.second()].getIsActive())
            {
                continue;
            }

            if (prtcInfoTable_.found(cPair))
            {
                std::vector<std::shared_ptr<prtSubContactInfo>>& subCList
                    = prtcInfoTable_[cPair]->getPrtSCList();

                for (auto sC : subCList)
                {
                    if (mag(sC->getprtCntVars().contactVolume_) > SMALL)
                    {
                        // lambda_min(I) per body; GREAT for static
                        // bodies, so their rotational bound never
                        // governs
                        const scalar iEffC
                        (
                            bodyIeff(immersedBodies_[cPair.first()])
                        );
                        const scalar iEffT
                        (
                            bodyIeff(immersedBodies_[cPair.second()])
                        );

                        // rotational bound of each body (tangential
                        // force on the lever arm vs lambda_min(I))
                        // enters the same min
                        const scalar dtPair
                        (
                            sC->getPairDtCritRot
                            (
                                iEffC,
                                iEffT,
                                dtTangentialFactor_,
                                dtRotationalFactor_
                            )
                        );

                        // the bound governs rotationally when it is
                        // below the translational one
                        const bool isRot
                        (
                            dtPair
                                < sC->getPairDtCrit(dtTangentialFactor_)
                        );

                        if (dtPair < deltaTDEM_)
                        {
                            deltaTDEM_ = dtPair;
                            govC = cPair.first();
                            govT = cPair.second();
                            govIsContact = 0;
                            govIsRot = isRot;
                        }
                    }
                }
            }
        }

        // sweep-inflated potential pairs minus hard-flagged members
        std::unordered_set<std::pair<label,label>, hashFunction> inflatedPairs;
        verletList_.computeInflatedPairs
        (
            immersedBodies_,
            sweepDist,
            hardFlags,
            inflatedPairs
        );

        for (const auto& pair : inflatedPairs)
        {
            const Tuple2<label, label> cPair
            (
                Tuple2<label, label>(pair.first, pair.second)
            );

            if (prtcInfoTable_.found(cPair)) continue;                  // contact: covered above

            immersedBody& cIb(immersedBodies_[pair.first]);
            immersedBody& tIb(immersedBodies_[pair.second]);

            if (!cIb.getIsActive() || !tIb.getIsActive()) continue;

            // exact pair materials
            const materialInfo& cMat(cIb.getibContactClass().getMatInfo());
            const materialInfo& tMat(tIb.getibContactClass().getMatInfo());

            const scalar aY
            (
                1.0/((1.0 - sqr(cMat.getNu()))/cMat.getY()
                    + (1.0 - sqr(tMat.getNu()))/tMat.getY())
            );
            const scalar aG
            (
                1.0/(2.0*(2.0 - cMat.getNu())*(1.0 + cMat.getNu())/cMat.getY()
                    + 2.0*(2.0 - tMat.getNu())*(1.0 + tMat.getNu())/tMat.getY())
            );

            const scalar Mc(cIb.getGeomModel().getM0());
            const scalar Mt(tIb.getGeomModel().getM0());
            const scalar reduceM
            (
                (Mc + Mt) > SMALL ? Mc*Mt/(Mc + Mt) : GREAT
            );

            // mesh scale of the contact zone
            const scalar h
            (
                min
                (
                    cIb.getCharCellSize(),
                    tIb.getCharCellSize()
                )
            );

            // equivalent radii from body masses and densities
            const scalar rC
            (
                demTimeStepInfo::equivRadius
                (
                    Mc/max(cIb.getGeomModel().getRhoS().value(), SMALL)
                )
            );
            const scalar rT
            (
                demTimeStepInfo::equivRadius
                (
                    Mt/max(tIb.getGeomModel().getRhoS().value(), SMALL)
                )
            );
            const scalar R(demTimeStepInfo::harmonicRadius(rC, rT));

            // approach-speed upper bound (Vel + omega*rMax per body)
            const scalar vN
            (
                pairApproachSpeed(cIb) + pairApproachSpeed(tIb)
            );

            const scalar aEst
            (
                demTimeStepInfo::aEst
                (
                    h,
                    R,
                    reduceM,
                    aY,
                    vN,
                    dtAreaCoeff_,
                    dtEstimatorVelocityAware_
                )
            );

            // pre-contact assumptions: contact volume = area over one
            // cell depth, Lc = h (the smallest possible Lc of any
            // contact zone spanning the cell - the most conservative
            // reading of the force laws)
            const scalar dtTrans
            (
                demTimeStepInfo::pairDtCrit
                (
                    demTimeStepInfo::effK
                    (
                        aY,
                        aEst*h,
                        aG,
                        aEst,
                        h,
                        dtTangentialFactor_
                    ),
                    reduceM
                )
            );

            // rotational bound: the same tangential stiffness
            // estimate acting on the lever-arm upper bound of each
            // body (bbox corner distance; 0 for static bodies) vs
            // lambda_min(I) (GREAT for static bodies, so the bound
            // never governs for them)
            const scalar kT
            (
                dtRotationalFactor_*dtTangentialFactor_
                    *demTimeStepInfo::contactKT(aG, aEst, h)
            );
            const scalar dtPair
            (
                min
                (
                    min
                    (
                        dtTrans,
                        demTimeStepInfo::rotationalDtCrit
                        (
                            kT, bodyRMax(cIb), bodyIeff(cIb)
                        )
                    ),
                    demTimeStepInfo::rotationalDtCrit
                    (
                        kT, bodyRMax(tIb), bodyIeff(tIb)
                    )
                )
            );

            if (dtPair < deltaTDEM_)
            {
                deltaTDEM_ = dtPair;
                govC = pair.first;
                govT = pair.second;
                govIsContact = 1;
                govIsRot = dtPair < dtTrans;
            }

            // hertz contact-duration diagnostic for the same pair
            const scalar dtH
            (
                demTimeStepInfo::hertzContactDt(reduceM, aY, R, vN)
            );
            if (dtH < dtDiagHertz_) dtDiagHertz_ = dtH;
        }

        // rayleigh diagnostic: min over sub-cycled bodies
        {
            const labelList subList(subCycle.toc());
            forAll(subList, sI)
            {
                immersedBody& cIb(immersedBodies_[subList[sI]]);

                if (!cIb.getIsActive()) continue;

                const scalar M0(cIb.getGeomModel().getM0());
                const materialInfo& cMat
                (
                    cIb.getibContactClass().getMatInfo()
                );

                // stiffest partner over all case materials: the
                // smallest rayleigh step (conservative diagnostic)
                scalar aGMin(GREAT);
                const HashTable<materialInfo,string,Hash<string>>& matProps
                (
                    materialProperties::getMatProps()
                );
                for
                (
                    auto mIter = matProps.cbegin();
                    mIter != matProps.cend();
                    ++mIter
                )
                {
                    const materialInfo& pMat(mIter());
                    aGMin = min
                    (
                        aGMin,
                        1.0/(2.0*(2.0 - cMat.getNu())*(1.0 + cMat.getNu())/cMat.getY()
                            + 2.0*(2.0 - pMat.getNu())*(1.0 + pMat.getNu())/pMat.getY())
                    );
                }

                const scalar R
                (
                    demTimeStepInfo::equivRadius
                    (
                        M0/max(cIb.getGeomModel().getRhoS().value(), SMALL)
                    )
                );

                // effective density of the equivalent sphere
                const scalar rhoEff
                (
                    M0/max
                    (
                        4.0/3.0*Foam::constant::mathematical::pi*pow3(R),
                        VSMALL
                    )
                );

                const scalar dtR
                (
                    demTimeStepInfo::rayleighDt
                    (
                        rhoEff,
                        aGMin,
                        cMat.getNu(),
                        R
                    )
                );
                if (dtR < dtDiagRayleigh_) dtDiagRayleigh_ = dtR;
            }
        }
    }

    //--- parallel reduce: min dt over ranks; the governing pair ids
    //    are rank-local, so the master picks the row of the minimum
    //    and scatters the winner
    {
        List<scalar> dtPerRank(Pstream::nProcs(), GREAT);
        List<labelList> govPerRank
        (
            Pstream::nProcs(),
            labelList(4, -1)
        );

        dtPerRank[Pstream::myProcNo()] = deltaTDEM_;
        govPerRank[Pstream::myProcNo()][0] = govIsContact;
        govPerRank[Pstream::myProcNo()][1] = govC;
        govPerRank[Pstream::myProcNo()][2] = govT;
        govPerRank[Pstream::myProcNo()][3] = govIsRot;

        Pstream::gatherList(dtPerRank, 0);
        Pstream::scatterList(dtPerRank, 0);
        Pstream::gatherList(govPerRank, 0);
        Pstream::scatterList(govPerRank, 0);

        label winRank(0);
        for (label rI = 1; rI < Pstream::nProcs(); ++rI)
        {
            if (dtPerRank[rI] < dtPerRank[winRank]) winRank = rI;
        }

        deltaTDEM_ = dtPerRank[winRank];
        govIsContact = govPerRank[winRank][0];
        govC = govPerRank[winRank][1];
        govT = govPerRank[winRank][2];
        govIsRot = govPerRank[winRank][3];

        reduce(dtDiagHertz_, minOp<scalar>());
        reduce(dtDiagRayleigh_, minOp<scalar>());
    }

    //--- floors and caps
    if (dtMinDEMdt_ > 0 && deltaTDEM_ < dtMinDEMdt_)
    {
        WarningInFunction
            << "estimated deltaTDEM " << deltaTDEM_
            << " below minDEMdt " << dtMinDEMdt_
            << ": clamping - the estimated stability bound is being "
            << "exceeded"
            << endl;
        deltaTDEM_ = dtMinDEMdt_;
    }

    //--- apply the stepDEM safety factor
    deltaTDEM_ *= stepDEM_;

    if (mesh_.time().timeIndex() % max(dtReportInterval_, 1) == 0)
    {
        if (govIsContact >= 0)
        {
            InfoH << DEM_Info << " deltaTDEM estimate: " << deltaTDEM_
                << " (pair " << govC << "-" << govT
                << ", "
                << (govIsContact == 0 ? "contact" : "pre-contact")
                << (govIsRot == 1 ? ", rotational" : "")
                << ")"
                << "; Hertz-based stepDEM: " << dtDiagHertz_
                << "; Rayleigh-based stepDEM: " << dtDiagRayleigh_
                << endl;
        }
        else
        {
            InfoH << DEM_Info << " deltaTDEM estimate: " << deltaTDEM_
                << "; Hertz-based stepDEM: " << dtDiagHertz_
                << "; Rayleigh-based stepDEM: " << dtDiagRayleigh_
                << endl;
        }
    }
}
//---------------------------------------------------------------------------//
scalar openHFDIBDEM::pairApproachSpeed(immersedBody& ib) const
{
    // upper bound on the contact-point speed magnitude (same terms as
    // computeSweepDistance minus the accel contribution)
    const ibContactVars& cVars(ib.getContactVars());

    return mag(cVars.Vel_) + mag(cVars.omega_)*bodyRMax(ib);
}
//---------------------------------------------------------------------------//
scalar openHFDIBDEM::bodyRMax(immersedBody& ib) const
{
    // static bodies do not rotate, so their contact points have no
    // lever arm
    if (ib.getbodyOperation() == 0) return 0;

    scalar rMax(0);
    boundBox bb(ib.getGeomModel().getBounds());
    pointField bbPoints(bb.points());
    vector CoM(ib.getGeomModel().getCoM());
    forAll(bbPoints, bP)
    {
        rMax = max(rMax, mag(bbPoints[bP] - CoM));
    }

    return rMax;
}
//---------------------------------------------------------------------------//
scalar openHFDIBDEM::bodyIeff(immersedBody& ib) const
{
    // static bodies are not integrated rotationally, so no rotational
    // stability step applies to them
    if (ib.getbodyOperation() == 0) return GREAT;

    // switched-off bodies have a zeroed mass state; their stale
    // contacts must not pull the estimate to zero
    if (!ib.getIsActive()) return GREAT;

    return demTimeStepInfo::minEigenvalue(ib.getGeomModel().getI());
}
//---------------------------------------------------------------------------//
prtContactInfo& openHFDIBDEM::getPrtcInfo(Tuple2<label,label> cPair)
{
    if(!prtcInfoTable_.found(cPair))
    {
        prtcInfoTable_.insert(cPair, autoPtr<prtContactInfo>( new prtContactInfo(
            immersedBodies_[cPair.first()].getibContactClass(),
            immersedBodies_[cPair.first()].getContactVars(),
            immersedBodies_[cPair.second()].getibContactClass(),
            immersedBodies_[cPair.second()].getContactVars()
        )));
    }

    return prtcInfoTable_[cPair]();
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::removeBodyContacts(label bodyId)
{
    // collect first, erase second: the OpenFOAM HashTable erase
    // invalidates the iterator (it returns bool, not the next
    // position like std::map)
    List<Tuple2<label,label>> keys;
    for (auto it = prtcInfoTable_.begin(); it != prtcInfoTable_.end(); ++it)
    {
        if (it.key().first() == bodyId || it.key().second() == bodyId)
        {
            keys.append(it.key());
        }
    }
    for (const auto& key : keys)
    {
        prtcInfoTable_.erase(key);
    }

    verletList_.removeBodyFromVList(immersedBodies_[bodyId]);

    InfoH << iB_Info << "-- body " << bodyId
        << " switched off: purged " << keys.size()
        << " contact pair(s), Verlet entry removed" << endl;
}
//---------------------------------------------------------------------------//
// function to either add or remove bodies from the simulation
void openHFDIBDEM::addRemoveBodies
(
    volScalarField& body,
    volVectorField& U,
    volScalarField& refineF
)
{
    nBodiesAddedLastStep_ = 0;

    forAll (addModels_,modelI)
    {
        word bodyName(bodyNames_[modelI]);

        label maxAdditions(50);
        label cAddition(0);

        // Note (MI): add solid body only if the number of immersed 
        //            bodies is below the threshold (if set)
        while (addModels_[modelI].shouldAddBody(body) and cAddition < maxAdditions and (solverInfo::getNSolidsThreshold() < 0 or immersedBodies_.size() < solverInfo::getNSolidsThreshold()))
        {
            InfoH << addModel_Info << "addModel invoked action, trying to add new body" << endl;
            std::shared_ptr<geomModel> bodyGeomModel(addModels_[modelI].addBody(body, immersedBodies_));

            cAddition++;

            if (addModels_[modelI].getBodyAdded())
            {
                InfoH << addModel_Info << "STL file correctly generated, registering the new body" << endl;

                // prepare pointer list for IBs (increase its size)
                label newIBSize(immersedBodies_.size()+1);
                label addIBPos(newIBSize - 1);
                immersedBodies_.setSize(newIBSize);

                // create the new body
                immersedBodies_.set
                (
                    addIBPos,
                    new immersedBody
                    (
                        bodyName,
                        mesh_,
                        HFDIBDEMDict_,
                        transportProperties_,
                        addIBPos,
                        recomputeM0_,
                        bodyGeomModel,
                        ibInterp_,
                        cellPoints_
                    )
                );

                // get reference for further processing
                immersedBody& nBody(immersedBodies_[addIBPos]);
                nBody.createImmersedBody(body,refineF);
                nBody.computeBodyCharPars();
                if (nBody.getStartSynced())
                {
                    nBody.initSyncWithFlow(U);
                }
                verletList_.addBodyToVList(nBody);

                InfoH << addModel_Info
                    << "new body included into the simulation" << endl;
                nBodiesAddedLastStep_++;
                cAddition = 0;
            }
            else
            {
                InfoH << addModel_Info
                    << "new body should have been added but was not "
                    << "(probably overlap with an existing body)"
                    << endl;
            }
        }
    }

    // make the counter rank-uniform
    reduce(nBodiesAddedLastStep_, maxOp<label>());
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::updateFSCoupling
(
    volScalarField& body,
    volVectorField& f
)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].pimpleUpdate(body,f,false,false);
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::updateFSCoupling
(
    volScalarField& body,
    volVectorField& f,
    const bool kinematicForce,
    const bool applyAddedMass
)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].pimpleUpdate(body,f,kinematicForce,applyAddedMass);
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::restartSimulation
(
    volScalarField& body,
    volScalarField& refineF,
    word runTime
)
{
    word timePath(recordOutDir_+"/"+runTime);
    fileNameList files(readDir(timePath));

    forAll(files,f)
    {
        IOdictionary bodyDict
        (
            IOobject
            (
                timePath + "/" + files[f],
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        word bodyId(std::to_string(readLabel(bodyDict.lookup("bodyId"))));
        word bodyName(bodyDict.lookup("bodyName"));
        vector Vel(bodyDict.lookup("Vel"));
        scalar omega(readScalar(bodyDict.lookup("omega")));
        vector Axis(bodyDict.lookup("Axis"));
        bool isStatic(readBool(bodyDict.lookup("static")));
        label timeStepsInContWStatic(readLabel(bodyDict.lookup("timeStepsInContWStatic")));

        std::shared_ptr<geomModel> bodyGeomModel;
        word bodyGeom;
        // check if the immersedDict_ contains bodyGeom
        if (HFDIBDEMDict_.subDict(bodyName).found("bodyGeom"))
        {
            word input = word(HFDIBDEMDict_.subDict(bodyName).lookup("bodyGeom"));
            bodyGeom = input;
            InfoH << iB_Info << "Found bodyGeom for "
                << bodyName << ", the body is: " << bodyGeom << endl;
        }
        else
        {
            bodyGeom = "convex";
            InfoH << iB_Info << "Did not find bodyGeom for "
                << bodyName << ", using bodyGeom: " << bodyGeom << endl;
        }

        if(bodyGeom == "convex")
        {
            word stlPath(timePath + "/stlFiles/"+bodyId+".stl");
            bodyGeomModel = std::make_shared<convexBody>(mesh_,stlPath);
        }
        else if(bodyGeom == "nonConvex")
        {
            word stlPath(timePath + "/stlFiles/"+bodyId+".stl");
            bodyGeomModel = std::make_shared<nonConvexBody>(mesh_,stlPath);
        }
        else if(bodyGeom == "sphere")
        {
            vector startPosition = vector(bodyDict.subDict("sphere").lookup("position"));
            scalar radius = readScalar(bodyDict.subDict("sphere").lookup("radius"));

            bodyGeomModel = std::make_shared<sphereBody>(mesh_,startPosition,radius);
        }
        else
        {
            word stlPath(timePath + "/stlFiles/"+bodyId+".stl");
            InfoH << iB_Info << "bodyGeom: " << bodyGeom
                << " not supported, using bodyGeom nonConvex" << endl;
            bodyGeom = "nonConvex";
            bodyGeomModel = std::make_shared<nonConvexBody>(mesh_,stlPath);
        }

        // propagate the per-body creation mode to the
        // restart-constructed geomModel (initializeIB.H path is not
        // taken here, so the dict entry would be silently ignored)
        {
            word bodyCreation("connectivity");
            if (HFDIBDEMDict_.subDict(bodyName).found("bodyCreation"))
            {
                bodyCreation = word
                    (HFDIBDEMDict_.subDict(bodyName)
                        .lookup("bodyCreation"));
                if
                (
                    bodyCreation != "connectivity"
                    && bodyCreation != "legacy"
                )
                {
                    InfoH << iB_Info << "Unknown bodyCreation: "
                        << bodyCreation
                        << ", using connectivity" << endl;
                    bodyCreation = "connectivity";
                }
            }
            bodyGeomModel->setBodyCreation(bodyCreation);
        }

        label newIBSize(immersedBodies_.size()+1);
        label addIBPos(newIBSize - 1);
        immersedBodies_.setSize(newIBSize);

        InfoH << iB_Info << "Restarting body: " << bodyId << " as "
            << addIBPos << " bodyName: " << bodyName << endl;
        immersedBodies_.set
        (
            addIBPos,
            new immersedBody
            (
                bodyName,
                mesh_,
                HFDIBDEMDict_,
                transportProperties_,
                addIBPos,
                recomputeM0_,
                bodyGeomModel,
                ibInterp_,
                cellPoints_
            )
        );

        immersedBodies_[addIBPos].createImmersedBody(body,refineF);
        immersedBodies_[addIBPos].computeBodyCharPars();
        immersedBodies_[addIBPos].setRestartSim(Vel,omega,Axis,isStatic,timeStepsInContWStatic);

        // coupling state (optional entries: restart files keep the
        //  constructor defaults)
        immersedBodies_[addIBPos].getCouplingModel().readCouplingInfo(bodyDict);

        verletList_.addBodyToVList(immersedBodies_[addIBPos]);
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::preCalculateCellPoints()
{
    cellPoints_.clear();
    cellPoints_.setSize(mesh_.nCells());
    forAll(mesh_.C(), cellI)
    {
        cellPoints_[cellI] = mesh_.cellPoints()[cellI];
    }

    forAll (immersedBodies_,bodyId)
    {
        immersedBodies_[bodyId].getGeomModel().resetHashTable();
        // mesh changed -> cached seeds and count baselines are invalid
        immersedBodies_[bodyId].getGeomModel().resetSeeds();
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::writeFirtsTimeBodiesInfo()
{
    word curOutDir(recordOutDir_ + "/" + mesh_.time().timeName());
    bool checkExistance(false);
    if(!recordSimulation_ || isDir(curOutDir))
        return;
    reduce(checkExistance,orOp<bool>());
    if(Pstream::myProcNo() == 0)
    {
        mkDir(curOutDir);
        mkDir(curOutDir +"/stlFiles");
    }
    reduce(checkExistance,orOp<bool>());

    DynamicLabelList activeIB;
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            activeIB.append(bodyId);
        }
    }

    wordList bodyNames;
    scalar listZize(activeIB.size());
    label bodiesPerProc = ceil(listZize/Pstream::nProcs());
    InfoH << basic_Info << "Active IB listZize      : " << listZize<< endl;
    InfoH << basic_Info << "bodiesPerProc : " << bodiesPerProc<< endl;

    for(int assignProc = Pstream::myProcNo()*bodiesPerProc; assignProc < min((Pstream::myProcNo()+1)*bodiesPerProc,activeIB.size()); assignProc++)
    {
        const label bodyId(activeIB[assignProc]);
        word path(curOutDir + "/body" + std::to_string(immersedBodies_[bodyId].getBodyId()) +".info");
        OFstream ofStream(path);
        IOobject outClass
            (
                path,
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            );
        IOdictionary outDict(outClass);

        outDict.writeHeader(ofStream);
        immersedBodies_[bodyId].recordBodyInfo(outDict,curOutDir);
        outDict.writeData(ofStream);
    }

}
//---------------------------------------------------------------------------//
void openHFDIBDEM::setSolverInfo()
{
    solverInfo::setOnlyDEM(true);
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::updateBodiesRhoF(const scalar rho)
{
    forAll (immersedBodies_,bodyId)
    {
        immersedBodies_[bodyId].updateRhoF(rho);
    }
}
void openHFDIBDEM::updateBodiesRhoF
(
    const volScalarField& rho,
    const volScalarField& body
)
{
    forAll (immersedBodies_,bodyId)
    {
        immersedBodies_[bodyId].updateRhoF(rho,body);
    }
}
//---------------------------------------------------------------------------//
scalar openHFDIBDEM::computeBodiesLinCourantNo()
{
    scalar maxCoNum(0);
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            maxCoNum = max(maxCoNum, immersedBodies_[bodyId].computeBodyLinCoNumber());
        }
    }
    return maxCoNum;
}
