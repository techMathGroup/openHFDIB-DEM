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
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*)
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

    dictionary demDic = HFDIBDEMDict_.subDict("DEM");
    dictionary materialsDic = demDic.subDict("materials");
    List<word> materialsNames = materialsDic.toc();
    forAll(materialsNames, matI)
    {
        dictionary matIDic = materialsDic.subDict(materialsNames[matI]);
        materialProperties::matProps_insert(
            materialsNames[matI],
            materialInfo(
                materialsNames[matI],
                readScalar(matIDic.lookup("Y")),
                readScalar(matIDic.lookup("nu")),
                readScalar(matIDic.lookup("gamma")),
                readScalar(matIDic.lookup("mu")),
                readScalar(matIDic.lookup("adhN"))
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

    dictionary patchDic = demDic.subDict("collisionPatches");
    List<word> patchNames = patchDic.toc();
    forAll(patchNames, patchI)
    {
        word patchMaterial = patchDic.subDict(patchNames[patchI]).lookup("material");
        vector patchNVec = patchDic.subDict(patchNames[patchI]).lookup("nVec");
        vector planePoint = patchDic.subDict(patchNames[patchI]).lookup("planePoint");

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
        dictionary cyclicPatchDic = demDic.subDict("cyclicPatches");
        List<word> cyclicPatchNames = cyclicPatchDic.toc();
        forAll(cyclicPatchNames, patchI)
        {
            vector patchNVec = cyclicPatchDic.subDict(cyclicPatchNames[patchI]).lookup("nVec");
            vector planePoint = cyclicPatchDic.subDict(cyclicPatchNames[patchI]).lookup("planePoint");
            word neighbourPatch = cyclicPatchDic.subDict(cyclicPatchNames[patchI]).lookup("neighbourPatch");

            cyclicPlaneInfo::insert(
                cyclicPatchNames[patchI],
                patchNVec,
                planePoint,
                neighbourPatch
            );
        }
    }

    if (HFDIBDEMDict_.found("geometricD"))
    {
        geometricD = HFDIBDEMDict_.lookup("geometricD");
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
        Info <<" -- VirtMesh Decomposition Level is set to        : "<< virtualMeshLevel::getVirtualMeshLevel() << endl;
        Info <<" -- VirtMesh charCellSize for boundary is set to  : "<< virtualMeshLevel::getCharCellSize() << endl;

    }
    else
    {
        virtualMeshLevel::setVirtualMeshLevel(1,1);
        Info <<" -- VirtMesh Decomposition Level is set to        : "<< virtualMeshLevel::getVirtualMeshLevel() << endl;
        Info <<" -- VirtMesh charCellSize for boundary is set to  : "<< virtualMeshLevel::getCharCellSize() << endl;

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

    // get data from HFDIBDEMDict
    //HFDIBinterpDict_ = HFDIBDEMDict_.subDict("interpolationSchemes");
    preCalculateCellPoints();

    if(HFDIBDEMDict_.found("interpolationSchemes"))
    {
        HFDIBinterpDict_ = HFDIBDEMDict_.subDict("interpolationSchemes");

        if(HFDIBinterpDict_.found("method"))
        {
            word intMethod = HFDIBinterpDict_.lookup("method");

            if(intMethod == "leastSquares")
            {
                dictionary lsCoeffsDict
                    = HFDIBinterpDict_.subDict("leastSquaresCoeffs");
                ibInterp_.set(new leastSquaresInt(
                    mesh_,
                    readScalar(lsCoeffsDict.lookup("distFactor")),
                    readScalar(lsCoeffsDict.lookup("radiusFactor")),
                    readScalar(lsCoeffsDict.lookup("angleFactor")),
                    readScalar(lsCoeffsDict.lookup("maxCCRows"))
                ));
            }
            else if(intMethod == "line")
            {
                ibInterp_.set(new lineInt(HFDIBinterpDict_));
            }
        }
    }

    bool startTime0(runTime == "0");

    // initialize addModels
    addModels_.setSize(bodyNames_.size());
    immersedBodies_.setSize(0);                                         //on the fly creation
    refineF *= 0;
    recomputeM0_ = recomputeM0;

    if(!startTime0)
    {
        if(!isDir(recordOutDir_))
            mkDir(recordOutDir_);
        else
        {
            fileNameList entries(readDir(recordOutDir_,fileType::directory)); // OF version 8, For version 6 use fileName::DIRECTORY instead of fileType::directory
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

        while (addModels_[modelI].shouldAddBody(body) and cAddition < maxAdditions)
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
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::createBodies(volScalarField& body,volScalarField& refineF)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].postContactUpdateBodyField(body,refineF);
        }
    }

    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].syncCreateImmersedBody(body,refineF);
            immersedBodies_[bodyId].checkIfInDomain(body);
            immersedBodies_[bodyId].updateOldMovementVars();
        }
    }

    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].chceckBodyOp();
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::preUpdateBodies
(
    volScalarField& body,
    volVectorField& f
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
    volScalarField& body,
    volVectorField& f
)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].clearIntpInfo();
            immersedBodies_[bodyId].postPimpleUpdateImmersedBody(body,f);
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::recreateBodies
(
    volScalarField& body,
    volScalarField& refineF
)
{
    refineF *= 0;
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
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].syncCreateImmersedBody(body,refineF);
            immersedBodies_[bodyId].checkIfInDomain(body);
            if(immersedBodies_[bodyId].getrecomputeM0() > 0)
            {
                immersedBodies_[bodyId].computeBodyCharPars();
                immersedBodies_[bodyId].recomputedM0();
            }
            InfoH << iB_Info << "-- body "
                << immersedBodies_[bodyId].getBodyId() << " Re-created" << endl;
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::interpolateIB( volVectorField & V
                              ,volVectorField & Vs
                              ,volScalarField & body)
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

    if(!isDir(curOutDir))
    {
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

                    scalar thrSurf(readScalar(HFDIBDEMDict_.lookup("surfaceThreshold")));
                    std::shared_ptr<periodicBody> newPeriodicBody
                        = std::make_shared<periodicBody>(mesh_, thrSurf);

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
    scalar step(stepDEM_);
    // scalar timeStep(step*deltaTime);

    while( pos < 1)
    {
        InfoH << DEM_Info << " Start DEM pos: " << pos
            << " DEM step: " << step << endl;

        InfoH << basic_Info << " DEM - CFD Time: "
            << mesh_.time().value() + deltaTime*pos << endl;

        forAll (immersedBodies_,ib)
        {
            immersedBodies_[ib].updateMovement(deltaTime*step*0.5);
            immersedBodies_[ib].moveImmersedBody(deltaTime*step);
        }

        verletList_.update(immersedBodies_);

        DynamicLabelList wallContactIB;
        DynamicList<wallSubContactInfo*> wallContactList;
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
                    if(detectWallContact
                    (
                        mesh_,
                        cIb.getibContactClass(),
                        cIb.getWallCntInfo()
                    ))
                    {
                        cIb.getibContactClass().setWallContact(true);
                        cIb.getibContactClass().inContactWithStatic(true);
                        wallContactIB.append(bodyId);
                        cIb.getWallCntInfo().registerSubContactList(wallContactList);
                    }
                }
            }
        }
        if(wallContactIB.size() > 0)
        {
            label wallContactPerProc(ceil(double(wallContactList.size())/Pstream::nProcs()));

            if( wallContactList.size() <= Pstream::nProcs())
            {
                wallContactPerProc = 1;
            }

            for(int assignProc = Pstream::myProcNo()*wallContactPerProc; assignProc < min((Pstream::myProcNo()+1)*wallContactPerProc,wallContactList.size()); assignProc++)
            {
                wallSubContactInfo* sCW = wallContactList[assignProc];
                immersedBody& cIb(immersedBodies_[sCW->getBodyId()]);
                sCW->setResolvedContact(solveWallContact(
                    mesh_,
                    cIb.getWallCntInfo(),
                    deltaTime*step,
                    *sCW
                    ));
            }
            // Pout <<" Survived " << endl;

            forAll (wallContactIB,iB)
            {
                immersedBody& cIb(immersedBodies_[wallContactIB[iB]]);

                std::vector<std::shared_ptr<wallSubContactInfo>>& subCList
                    = cIb.getWallCntInfo().getWallSCList();

                for(auto sC : subCList)
                {
                    sC->syncContactResolve();
                    if(sC->getContactResolved())
                    {
                        sC->syncData();

                        cIb.updateContactForces(
                            sC->getOutForce()
                        );
                    }
                }
                cIb.getWallCntInfo().clearOldContact();
            }
            // Pout <<" Survived 2 " << endl;
        }
        wallContactList.clear();
        wallContactIB.clear();

        DynamicList<prtSubContactInfo*> contactList;
        // check only pairs whose bounding boxes are intersected for the contact
        for (auto it = verletList_.begin(); it != verletList_.end(); ++it)
        {
            const Tuple2<label, label> cPair = Tuple2<label, label>(it->first, it->second);

            label cInd(cPair.first());
            bool cStatic(immersedBodies_[cInd].getbodyOperation() == 0);

            label tInd(cPair.second());
            bool tStatic(immersedBodies_[tInd].getbodyOperation() == 0);

            if((immersedBodies_[cInd].getIsActive() && immersedBodies_[tInd].getIsActive())
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

                prtcInfo.syncContactList();

                prtcInfo.registerContactList(contactList);
            }
        }

        if(contactList.size() > 0 )
        {
            label contactPerProc(ceil(double(contactList.size())/Pstream::nProcs()));
            if( contactList.size() <= Pstream::nProcs())
            {
                contactPerProc = 1;
            }

            for(int assignProc = Pstream::myProcNo()*contactPerProc; assignProc < min((Pstream::myProcNo()+1)*contactPerProc,contactList.size()); assignProc++)
            {
                prtSubContactInfo* sCI = contactList[assignProc];
                const Tuple2<label, label>& cPair = sCI->getCPair();

                ibContactClass& cClass(immersedBodies_[cPair.first()].getibContactClass());
                ibContactClass& tClass(immersedBodies_[cPair.second()].getibContactClass());

                if(detectPrtPrtContact(mesh_,cClass,tClass,*sCI))
                {
                    prtContactInfo& prtcInfo(getPrtcInfo(cPair));
                    sCI->setResolvedContact(solvePrtContact(mesh_, prtcInfo, *sCI, deltaTime*step));
                }
            }
        }

        for (auto it = verletList_.begin(); it != verletList_.end(); ++it)
        {
            const Tuple2<label, label> cPair = Tuple2<label, label>(it->first, it->second);
            label cInd(cPair.first());
            label tInd(cPair.second());

            prtContactInfo& prtcInfo(getPrtcInfo(cPair));
            if(!prtcInfo.contactResolved())
            {
                if(prtcInfoTable_.found(cPair))
                {
                    prtcInfoTable_.erase(cPair);
                    continue;
                }
            }

            std::vector<std::shared_ptr<prtSubContactInfo>>& subCList
                = prtcInfo.getPrtSCList();

            for(auto sC : subCList)
            {
                sC->syncData();

                immersedBodies_[cInd].updateContactForces
                (
                    sC->getOutForce().first()
                );

                immersedBodies_[tInd].updateContactForces
                (
                    sC->getOutForce().second()
                );
            }
        }

        scalar maxCoNum = 0;
        label  bodyId = 0;
        forAll (immersedBodies_,ib)
        {
            immersedBodies_[ib].updateMovement(deltaTime*step*0.5);

            immersedBodies_[ib].computeBodyCoNumber();
            if (maxCoNum < immersedBodies_[ib].getCoNum())
            {
                maxCoNum = immersedBodies_[ib].getCoNum();
                bodyId = ib;
            }
        }
        InfoH << basic_Info << "Max CoNum = " << maxCoNum << " at body " << bodyId << endl;

        pos += step;

        if (pos + step + SMALL >= 1)
            step = 1 - pos;
    }
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
// function to either add or remove bodies from the simulation
void openHFDIBDEM::addRemoveBodies
(
    volScalarField& body,
    volVectorField& U,
    volScalarField& refineF
)
{
    forAll (addModels_,modelI)
    {
        word bodyName(bodyNames_[modelI]);

        label maxAdditions(50);
        label cAddition(0);

        while (addModels_[modelI].shouldAddBody(body) and cAddition < maxAdditions)
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
            immersedBodies_[bodyId].pimpleUpdate(body,f);
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
    scalar thrSurf(readScalar(HFDIBDEMDict_.lookup("surfaceThreshold")));

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
            word input = HFDIBDEMDict_.subDict(bodyName).lookup("bodyGeom");
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
            bodyGeomModel = std::make_shared<convexBody>(mesh_,stlPath,thrSurf);
        }
        else if(bodyGeom == "nonConvex")
        {
            word stlPath(timePath + "/stlFiles/"+bodyId+".stl");
            bodyGeomModel = std::make_shared<nonConvexBody>(mesh_,stlPath,thrSurf);
        }
        else if(bodyGeom == "sphere")
        {
            vector startPosition = bodyDict.subDict("sphere").lookup("position");
            scalar radius = readScalar(bodyDict.subDict("sphere").lookup("radius"));

            bodyGeomModel = std::make_shared<sphereBody>(mesh_,startPosition,radius,thrSurf);
        }
        else
        {
            word stlPath(timePath + "/stlFiles/"+bodyId+".stl");
            InfoH << iB_Info << "bodyGeom: " << bodyGeom
                << " not supported, using bodyGeom nonConvex" << endl;
            bodyGeom = "nonConvex";
            bodyGeomModel = std::make_shared<nonConvexBody>(mesh_,stlPath,thrSurf);
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
    }
}
//---------------------------------------------------------------------------//
