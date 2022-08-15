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
//#include "SVD.H"
#include "scalarMatrices.H"
#include "OFstream.H"
#include "defineExternVars.H"

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
    dictionary demDic = HFDIBDEMDict_.subDict("DEM");
    dictionary materialsDic = demDic.subDict("materials");
    List<word> materialsNames = materialsDic.toc();
    forAll(materialsNames, matI)
    {
        dictionary matIDic = materialsDic.subDict(materialsNames[matI]);
        materialInfos_.insert(
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

            matInterAdh_.insert(
                interKey,
                readScalar(interDicI.lookup("value"))
            );
        }
    }

    dictionary patchDic = demDic.subDict("collisionPatches");
    List<word> patchNames = patchDic.toc();
    forAll(patchNames, patchI)
    {
        word patchMaterial = patchDic.lookup(patchNames[patchI]);
        wallInfos_.insert(
            patchNames[patchI],
            materialInfos_[patchMaterial]
        );
    }

    if(demDic.found("cyclicPatches"))
    {
        cyclicPatches_ = wordList(demDic.lookup("cyclicPatches"));
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

    recordOutDir_ = mesh_.time().rootPath() + "/" + mesh_.time().globalCaseName() + "/bodiesInfo";
}
//---------------------------------------------------------------------------//
openHFDIBDEM::~openHFDIBDEM()
{
}
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
        InfoH.setOutput(basicOutput,iBoutput,DEMoutput,addModelOutput);
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
            autoPtr<geomModel> bodyGeomModel(addModels_[modelI].addBody(body, immersedBodies_));
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
                        body,
                        HFDIBDEMDict_,
                        transportProperties_,
                        addIBPos,
                        recomputeM0_,
                        bodyGeomModel.ptr(),
                        ibInterp_,
                        cellPoints_,
                        materialInfos_,
                        wallInfos_,
                        matInterAdh_
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
    mkDir(curOutDir);
    mkDir(curOutDir +"/stlFiles");

    wordList bodyNames;
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
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
    if (cyclicPatches_.size() > 0)
    {
        forAll (immersedBodies_,bodyId)
        {
            if (!immersedBodies_[bodyId].getGeomModel().isCluster())
            {
                vector newPos = vector::zero;
                bool inContact = detectCyclicContact(mesh_, cyclicPatches_, immersedBodies_[bodyId].getWallCntInfo(), newPos);
                if (inContact)
                {
                    verletList_.removeBodyFromVList(immersedBodies_[bodyId]);

                    scalar thrSurf(readScalar(HFDIBDEMDict_.lookup("surfaceThreshold")));
                    autoPtr<periodicBody> newClusterBody(new periodicBody(mesh_, thrSurf));
                    autoPtr<geomModel> iBcopy(immersedBodies_[bodyId].getGeomModel().getGeomModel());
                    vector transVec = newPos - iBcopy().getCoM();
                    iBcopy->bodyMovePoints(transVec);
                    newClusterBody->addBodyToCluster(immersedBodies_[bodyId].getGeomModelPtr());
                    newClusterBody->addBodyToCluster(iBcopy);
                    geomModel* clusterGeomModel = newClusterBody.ptr();
                    immersedBodies_[bodyId].getGeomModelPtr().set(clusterGeomModel);

                    verletList_.addBodyToVList(immersedBodies_[bodyId]);
                }
            }
            else
            {
                periodicBody& cBody = dynamic_cast<periodicBody&>(immersedBodies_[bodyId].getGeomModel());

                if(cBody.shouldBeUnclustered())
                {
                    verletList_.removeBodyFromVList(immersedBodies_[bodyId]);

                    immersedBodies_[bodyId].getGeomModelPtr().reset(cBody.getRemGeomModel().ptr());

                    verletList_.addBodyToVList(immersedBodies_[bodyId]);
                }
            }
        }
    }

    scalar deltaTime(mesh_.time().deltaT().value());
    scalar pos(0.0);
    scalar step(stepDEM_);

    while( pos < 1)
    {
        InfoH << DEM_Info << " Start DEM pos: " << pos
            << " DEM step: " << step << endl;

        verletList_.update(immersedBodies_);

        forAll (immersedBodies_,bodyId)
        {
            immersedBody& cIb(immersedBodies_[bodyId]);
            if (cIb.getIsActive())
            {
                // set F_ and T_ to zero.
                cIb.resetForces();
                // add fluid coupling force to F and T.
                cIb.updateFAndT(
                    cIb.getHistoryCouplingF(),
                    cIb.getHistoryCouplingT()
                );

                if(cIb.getbodyOperation() != 0)
                {
                    // detect wall contact
                    if(detectWallContact
                    (
                        mesh_,
                        cIb.getibContactClass()
                    ))
                    {
                        cIb.getibContactClass().setWallContact(true);
                        cIb.getibContactClass().inContactWithStatic(true);
                    }
                    // solve wall contact
                    if (cIb.checkWallContact())
                    {
                        solveWallContact
                        (
                            mesh_,
                            cIb.getWallCntInfo(),
                            deltaTime*step
                        );

                        cIb.updateFAndT
                        (
                            cIb.getWallCntInfo().getOutForce().F,
                            cIb.getWallCntInfo().getOutForce().T
                        );
                    }
                }
            }
        }

        // check only pairs whose bounding boxes are intersected for the contact
        for (auto it = verletList_.begin(); it != verletList_.end(); ++it)
        {
            const Tuple2<label, label>& cPair = it.key();

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

                ibContactClass& cClass(immersedBodies_[cInd].getibContactClass());
                ibContactClass& tClass(immersedBodies_[tInd].getibContactClass());

                if(detectPrtPrtContact(
                    mesh_,
                    cClass,
                    tClass
                ))
                {
                    prtContactInfo& prtcInfo(getPrtcInfo(
                        cPair)
                    );

                    if(solvePrtContact(
                        mesh_,
                        prtcInfo,
                        deltaTime*step
                    ))
                    {
                        immersedBodies_[cInd].updateFAndT
                        (
                            prtcInfo.getOutForce().first().F,
                            prtcInfo.getOutForce().first().T
                        );
                        immersedBodies_[tInd].updateFAndT
                        (
                            prtcInfo.getOutForce().second().F,
                            prtcInfo.getOutForce().second().T
                        );
                    }
                }
                else
                {
                    if(prtcInfoTable_.found(cPair))
                    {
                        prtcInfoTable_.erase(cPair);
                    }
                }
            }
        }

        forAll (immersedBodies_,ib)
        {
            immersedBodies_[ib].updateMovement(deltaTime*step);
            immersedBodies_[ib].moveImmersedBody(deltaTime*step);
        }

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
        prtcInfoTable_.insert(cPair, prtContactInfo(
            immersedBodies_[cPair.first()].getibContactClass(),
            immersedBodies_[cPair.first()].getContactVars(),
            immersedBodies_[cPair.second()].getibContactClass(),
            immersedBodies_[cPair.second()].getContactVars(),
            matInterAdh_
        ));
    }

    return prtcInfoTable_[cPair];
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
            autoPtr<geomModel> bodyGeomModel(addModels_[modelI].addBody(body, immersedBodies_));

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
                        body,
                        HFDIBDEMDict_,
                        transportProperties_,
                        addIBPos,
                        recomputeM0_,
                        bodyGeomModel.ptr(),
                        ibInterp_,
                        cellPoints_,
                        materialInfos_,
                        wallInfos_,
                        matInterAdh_
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

        autoPtr<geomModel> bodyGeomModel;
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
            bodyGeomModel.set(new convexBody(mesh_,stlPath,thrSurf));
        }
        else if(bodyGeom == "nonConvex")
        {
            word stlPath(timePath + "/stlFiles/"+bodyId+".stl");
            bodyGeomModel.set(new nonConvexBody(mesh_,stlPath,thrSurf));
        }
        else if(bodyGeom == "sphere")
        {
            vector startPosition = bodyDict.subDict("sphere").lookup("position");
            scalar radius = readScalar(bodyDict.subDict("sphere").lookup("radius"));

            bodyGeomModel.set(new sphereBody(mesh_,startPosition,radius,thrSurf));
        }
        else
        {
            word stlPath(timePath + "/stlFiles/"+bodyId+".stl");
            InfoH << iB_Info << "bodyGeom: " << bodyGeom
                << " not supported, using bodyGeom nonConvex" << endl;
            bodyGeom = "nonConvex";
            bodyGeomModel.set(new nonConvexBody(mesh_,stlPath,thrSurf));
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
                body,
                HFDIBDEMDict_,
                transportProperties_,
                addIBPos,
                recomputeM0_,
                bodyGeomModel.ptr(),
                ibInterp_,
                cellPoints_,
                materialInfos_,
                wallInfos_,
                matInterAdh_
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
