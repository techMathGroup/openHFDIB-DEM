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
    Martin Isoz (2019-*), Martin Šourek (2019-*)
\*---------------------------------------------------------------------------*/
#include "openHFDIBDEM.H"
#include "polyMesh.H"
#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"

#include "interpolationCellPoint.H"
#include "interpolationCell.H"
#include "SVD.H"
#include "scalarMatrices.H"
#include "OFstream.H"

#define ORDER 2

using namespace Foam;
using namespace contactModel;

//---------------------------------------------------------------------------//
openHFDIBDEM::openHFDIBDEM(const Foam::dynamicFvMesh& mesh )
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
minDEMloops_(readScalar(HFDIBDEMDict_.lookup("minDEMloops"))),
minDEMtimeStep_(readScalar(HFDIBDEMDict_.lookup("minDEMtimeStep"))),
recordSimulation_(readBool(HFDIBDEMDict_.lookup("recordSimulation")))
{
    scalar adhN(0.0);
    if (HFDIBDEMDict_.subDict("wallProps").found("adhN"))
    {
        adhN = readScalar(HFDIBDEMDict_.subDict("wallProps").lookup("adhN"));
    }
    scalar adhEqui(0.0);
    if (HFDIBDEMDict_.subDict("wallProps").found("adhEqui"))
    {
        adhEqui = readScalar(HFDIBDEMDict_.subDict("wallProps").lookup("adhEqui"));
    }

    wallInfo_.set(
        new wallInfo(
            readScalar(HFDIBDEMDict_.subDict("wallProps").lookup("kN")),
            readScalar(HFDIBDEMDict_.subDict("wallProps").lookup("kt")),
            readScalar(HFDIBDEMDict_.subDict("wallProps").lookup("gammaN")),
            readScalar(HFDIBDEMDict_.subDict("wallProps").lookup("gammat")),
            readScalar(HFDIBDEMDict_.subDict("wallProps").lookup("mu")),
            adhN,
            adhEqui
                  )
        );

    if (HFDIBDEMDict_.found("geometricD"))
    {
        geometricD_ = HFDIBDEMDict_.lookup("geometricD");
    }
    else
    {
        geometricD_ = mesh_.geometricD();
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
    // get data from HFDIBDEMDict
    HFDIBinterpDict_ = HFDIBDEMDict_.subDict("interpolationSchemes");
    preCalculateCellPoints();

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
        Info << "Creating immersed body based on: " << bodyName << endl;

        label maxAdditions(1000);
        label cAddition(0);

        while (addModels_[modelI].shouldAddBody(body) and cAddition < maxAdditions)
        {
            Info << "addModel invoked action, trying to add new body" << endl;
            autoPtr<geomModel> bodyGeomModel(addModels_[modelI].addBody(body));
            cAddition++;

            // initialize the immersed bodies
            if (addModels_[modelI].getBodyAdded())
            {
                label newIBSize(immersedBodies_.size()+1);
                label addIBPos(newIBSize - 1);
                immersedBodies_.setSize(newIBSize);

                Info << "Trying to set immersedBodies" << endl;
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
                        geometricD_,
                        bodyGeomModel.ptr(),
                        cellPoints_
                    )
                );
                immersedBodies_[addIBPos].createImmersedBody(body,refineF);
                immersedBodies_[addIBPos].computeBodyCharPars();
                if (immersedBodies_[addIBPos].getStartSynced())
                {
                    immersedBodies_[addIBPos].initSyncWithFlow(U);
                }
                Info << "Body based on: " << bodyName << " successfully added" << endl;
                cAddition = 0;
            }
            else
            {
                Info << "Body based on: " << bodyName << " should have been added but was not "
                     << "(probably overlap with an already existing body)"
                     << endl;
            }
        }
    }

    const dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/pow(average(mesh_.V()), 1.0/3.0)
    );
    surfNorm_ = -fvc::grad(body);
    surfNorm_ /= (mag(surfNorm_)+deltaN.value());

    // initialize list for neighbour list method based on number of IBs
    boundValueNeighbourList_.setSize(3);
    boundLabelNeighbourList_.setSize(3);
    contactInCoordNeighbourList_.setSize(3);
    numberOfDEMloops_.setSize(Pstream::nProcs());
    for (label i = 0; i < 3; i = i + 1)
    {
        boundValueNeighbourList_[i].setSize(2 * immersedBodies_.size());
        boundLabelNeighbourList_[i].setSize(2 * immersedBodies_.size());
    }

    forAll (immersedBodies_,bodyI)
    {
        // get references to bounding points and label list
        List<List<label>> IBboundList = immersedBodies_[bodyI].getBoundIndList();
        vector IBboundMinPoint = immersedBodies_[bodyI].getMinBoundPoint();
        vector IBboundMaxPoint = immersedBodies_[bodyI].getMaxBoundPoint();
        // prepere label for bounding point identification. This is only list where bodyID is increased by
        // minimal bounding point has negative value
        // maximal bounding point has positive value
        label bodyIdInc(bodyI+1);

        // iterate over dimensions and assined proper bounding position with its label
        for (label coord = 0; coord < 3; coord = coord + 1)
        {
            boundValueNeighbourList_[coord][IBboundList[coord][0]] = IBboundMinPoint[coord];
            boundLabelNeighbourList_[coord][IBboundList[coord][0]] = -1 * bodyIdInc;
            boundValueNeighbourList_[coord][IBboundList[coord][1]] = IBboundMaxPoint[coord];
            boundLabelNeighbourList_[coord][IBboundList[coord][1]] = bodyIdInc;
        }
    }
    
    sortBoundingListPrtContact();
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::preUpdateBodies
(
    volScalarField& body,
    volVectorField& f
)
{
    // clear particle-particle contact data
    ibContactList_.clear();
    prtContactIBList_.clear();

    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            // create body or compute body-fluid coupling and estimate
            // potential contacts with walls
            immersedBodies_[bodyId].inContactWithStatic(false);

            immersedBodies_[bodyId].preContactUpdateImmersedBody(body,f);

            if(immersedBodies_[bodyId].shouldDetectWallContact())
                detectWallContact(mesh_,immersedBodies_[bodyId].getContactInfo());

            immersedBodies_[bodyId].printStats();

            // check for body-wall contact
            if (immersedBodies_[bodyId].checkWallContact())
            {
                ibContactList_.append(immersedBodies_[bodyId].getBodyId());
            }
        }
    }

    // Update neighbour list and detect possible prt-prt contact
    updateNeighbourLists();
    detectPrtContact();
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
            immersedBodies_[bodyId].postContactUpdateImmersedBody(body,f);
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::moveBodies
(
    volScalarField& body,
    volScalarField& refineF
)
{
    refineF *= 0;
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            if (findIndex(ibContactList_, bodyId) == -1)
            {
                immersedBodies_[bodyId].moveImmersedBody();
            }
        }
    }
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            if (findIndex(ibContactList_, bodyId) == -1)
            {
                immersedBodies_[bodyId].postContactUpdateBodyField(body,refineF);
            }
        }
    }
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            if (findIndex(ibContactList_, bodyId) == -1)
            {
                immersedBodies_[bodyId].syncCreateImmersedBody(body,refineF);
                immersedBodies_[bodyId].checkIfInDomain(body);
                immersedBodies_[bodyId].updateOldMovementVars();
                immersedBodies_[bodyId].printStats();
            }
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
            Info << "-- body " << immersedBodies_[bodyId].getBodyId() << " Re-created" << endl;
        }
    }
    const dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/pow(average(mesh_.V()), 1.0/3.0)
    );
    surfNorm_ = -fvc::grad(body);
    surfNorm_ /= (mag(surfNorm_)+deltaN.value());
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::interpolateIB( volVectorField & V
                              ,volVectorField & Vs
                              ,volScalarField & body)
{

    // create interpolator
    autoPtr<interpolation<vector>> interpV =
                   interpolation<vector>::New(HFDIBinterpDict_, V);
    // reset imposed field
    Vs *= scalar(0);

    // loop over all the immersed bodies
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            // update imposed field according to body
            immersedBodies_[bodyId].updateVectorField(Vs, V.name(),body, surfNorm_);

            if(immersedBodies_[bodyId].getUseInterpolation())
            {
                // get interpolation info a request list
                const List<DynamicList<immersedBody::intVecRequest>>& intVecReqList = immersedBodies_[bodyId].getinterpolationVecReqs();;
                List<DynamicList<immersedBody::interpolationInfo>>& intInfoList  = immersedBodies_[bodyId].getInterpolationInfo();

                List<DynamicPointList> intPointsToSend;
                List<DynamicLabelList> intCellsToSend;
                intPointsToSend.setSize(Pstream::nProcs());
                intCellsToSend.setSize(Pstream::nProcs());
                // deal with all requests and ask processors for interpolation
                // pstreamBuffers works only with list so requests are rearenged
                for (label proci = 0; proci < Pstream::nProcs(); proci++)
                {
                    if (proci != Pstream::myProcNo())
                    {
                        for (label i = 0; i < intVecReqList[proci].size(); i++)
                        {
                            intPointsToSend[proci].append(intVecReqList[proci][i].intPoint_);
                            intCellsToSend[proci].append(intVecReqList[proci][i].intCell_);
                        }
                    }
                }

                PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
                PstreamBuffers pBufs2(Pstream::commsTypes::nonBlocking);

                for (label proci = 0; proci < Pstream::nProcs(); proci++)
                {
                    if (proci != Pstream::myProcNo())
                    {
                        UOPstream send(proci, pBufs);
                        send << intPointsToSend[proci];

                        UOPstream send2(proci, pBufs2);
                        send2 << intCellsToSend[proci];
                    }
                }

                pBufs.finishedSends();
                pBufs2.finishedSends();

                List<DynamicPointList> intPointsRcv;
                intPointsRcv.setSize(Pstream::nProcs());
                List<DynamicLabelList> intCellsRcv;
                intCellsRcv.setSize(Pstream::nProcs());

                List<DynamicVectorList> intVecToReturn;
                intVecToReturn.setSize(Pstream::nProcs());
                // receive interpolation requests from other processors
                for (label proci = 0; proci < Pstream::nProcs(); proci++)
                {
                    if (proci != Pstream::myProcNo())
                    {
                        UIPstream recv(proci, pBufs);
                        DynamicPointList recIntPoints (recv);
                        intPointsRcv[proci].append(recIntPoints);

                        UIPstream recv2(proci, pBufs2);
                        DynamicLabelList recIntCells (recv2);
                        intCellsRcv[proci].append(recIntCells);
                    }
                }
                // compute interpolation for other processors
                for (label otherProci = 0; otherProci < intPointsRcv.size(); otherProci++)
                {
                    for (label intReqI = 0; intReqI < intPointsRcv[otherProci].size(); intReqI++)
                    {
                        vector VP1 =  interpV->interpolate(intPointsRcv[otherProci][intReqI],
                                                            intCellsRcv[otherProci][intReqI]
                                                            );
                        intVecToReturn[otherProci].append(VP1);
                    }
                }

                pBufs.clear();
                // send computed data back
                for (label proci = 0; proci < Pstream::nProcs(); proci++)
                {
                    if (proci != Pstream::myProcNo())
                    {
                        UOPstream send(proci, pBufs);
                        send << intVecToReturn[proci];
                    }
                }

                pBufs.finishedSends();

                List<DynamicVectorList> intVecRcv;
                intVecRcv.setSize(Pstream::nProcs());
                // receive interpolation data from other processors
                for (label proci = 0; proci < Pstream::nProcs(); proci++)
                {
                    if (proci != Pstream::myProcNo())
                    {
                        UIPstream recv(proci, pBufs);
                        DynamicVectorList recIntVec (recv);
                        intVecRcv[proci].append(recIntVec);
                    }
                }
                // assign received data
                for (label otherProci = 0; otherProci < intVecRcv.size(); otherProci++)
                {
                    for (label intVecI = 0; intVecI < intVecRcv[otherProci].size(); intVecI++)
                    {
                        intInfoList[Pstream::myProcNo()][intVecReqList[otherProci][intVecI].requestLabel_].intVec_[intVecReqList[otherProci][intVecI].vecLabel_] = intVecRcv[otherProci][intVecI];
                    }
                }
                // compute interpolation for cells found on this processor
                forAll (intInfoList[Pstream::myProcNo()],infoI)
                {
                    switch(intInfoList[Pstream::myProcNo()][infoI].order_)
                    {
                        case 1:
                        {
                            if (intInfoList[Pstream::myProcNo()][infoI].procWithIntCells_[0] == Pstream::myProcNo())
                            {
                                vector VP1 =  interpV->interpolate(intInfoList[Pstream::myProcNo()][infoI].intPoints_[1], intInfoList[Pstream::myProcNo()][infoI].intCells_[0]
                                                            );
                                intInfoList[Pstream::myProcNo()][infoI].intVec_[0] = VP1;
                            }
                            break;
                        }
                        case 2:
                        {
                            if (intInfoList[Pstream::myProcNo()][infoI].procWithIntCells_[0] == Pstream::myProcNo())
                            {
                                vector VP1 =  interpV->interpolate(intInfoList[Pstream::myProcNo()][infoI].intPoints_[1], intInfoList[Pstream::myProcNo()][infoI].intCells_[0]
                                                            );
                                intInfoList[Pstream::myProcNo()][infoI].intVec_[0] = VP1;
                            }

                            if (intInfoList[Pstream::myProcNo()][infoI].procWithIntCells_[1] == Pstream::myProcNo())
                            {
                                vector VP2 =  interpV->interpolate(intInfoList[Pstream::myProcNo()][infoI].intPoints_[2], intInfoList[Pstream::myProcNo()][infoI].intCells_[1]
                                                            );
                                intInfoList[Pstream::myProcNo()][infoI].intVec_[1] = VP2;
                            }
                            break;
                        }
                    }
                }

                // loop over all interpolation info
                forAll (intInfoList[Pstream::myProcNo()],infoI)
                {
                    label cellI = intInfoList[Pstream::myProcNo()][infoI].surfCell_;
                    // based on order calculate Vs and assign
                    switch(intInfoList[Pstream::myProcNo()][infoI].order_)
                    {
                        case 0:
                        {
                            Vs[cellI] = body[cellI]*Vs[cellI]  + (1.0-body[cellI])*V[cellI];
                            break;
                        }
                        case 1:
                        {
                            vector VP1 = intInfoList[Pstream::myProcNo()][infoI].intVec_[0] - Vs[cellI];

                            // distance between interpolation points
                            scalar deltaR = mag(intInfoList[Pstream::myProcNo()][infoI].intPoints_[1] -intInfoList[Pstream::myProcNo()][infoI].intPoints_[0]);

                            // cell center to surface distance
                            scalar ds(0.0);
                            if (immersedBodies_[bodyId].getSDBasedLambda())
                            {
                                scalar minMaxBody(max(min(body[cellI],1.0-SMALL),SMALL));
                                ds = Foam::atanh(2.0*minMaxBody - 1.0)*Foam::pow(mesh_.V()[cellI],0.333);
                                ds/= immersedBodies_[bodyId].getIntSpan();
                            }
                            else
                            {
                                ds = deltaR*(0.5-body[cellI]) ;
                            }

                            vector linCoeff = VP1/(deltaR+SMALL);

                            Vs[cellI] = linCoeff*ds + Vs[cellI];
                            break;
                        }
                        case 2:
                        {
                            vector VP1 =  intInfoList[Pstream::myProcNo()][infoI].intVec_[0] - Vs[cellI];

                            vector VP2 =  intInfoList[Pstream::myProcNo()][infoI].intVec_[1] - Vs[cellI];


                            // distance between interpolation points
                            scalar deltaR1 = mag(intInfoList[Pstream::myProcNo()][infoI].intPoints_[2] -intInfoList[Pstream::myProcNo()][infoI].intPoints_[1]);
                            scalar deltaR2 = mag(intInfoList[Pstream::myProcNo()][infoI].intPoints_[1] -intInfoList[Pstream::myProcNo()][infoI].intPoints_[0]);

                            // cell center to surface distance
                            scalar ds(0.0);
                            if (immersedBodies_[bodyId].getSDBasedLambda())
                            {
                                scalar minMaxBody(max(min(body[cellI],1.0-SMALL),SMALL));
                                ds = Foam::atanh(2.0*minMaxBody - 1.0)*Foam::pow(mesh_.V()[cellI],0.333);
                                ds/= immersedBodies_[bodyId].getIntSpan();
                            }
                            else
                            {
                                if((0.5-body[cellI]) < 0)
                                {
                                    ds = -1*mag(mesh_.C()[cellI] - intInfoList[Pstream::myProcNo()][infoI].intPoints_[0]);
                                }
                                else
                                {
                                    ds = mag(mesh_.C()[cellI] - intInfoList[Pstream::myProcNo()][infoI].intPoints_[0]);
                                }
                            }

                            vector quadCoeff = (VP2 - VP1)*deltaR1 - VP1*deltaR2;
                            quadCoeff       /= (deltaR1*deltaR2*(deltaR1 + deltaR2)+SMALL);

                            vector linCoeff  = (VP1-VP2)*Foam::pow(deltaR1,2.0);
                            linCoeff        += 2.0*VP1*deltaR1*deltaR2;
                            linCoeff        += VP1*Foam::pow(deltaR2,2.0);
                            linCoeff        /= (deltaR1*deltaR2*(deltaR1 + deltaR2)+SMALL);

                            Vs[cellI] = quadCoeff*ds*ds + linCoeff*ds + Vs[cellI];
                            break;
                        }
                    }
                }
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
            IOobject outInfo
                (
                    path,
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                );
            IOdictionary outDict(outInfo);

            outDict.writeHeader(ofStream);
            immersedBodies_[bodyId].recordBodyInfo(outDict,curOutDir);
            outDict.writeData(ofStream);
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::correctContact(volScalarField& body,volScalarField& refineF)
{
    // first detect prt contact to be sure that all Contact during CFD step will be solved during DEM inner loops
    detectPrtContact();
    // get new contacts with wall and return their position
    getContactListAndReturnPositions(body);
    scalar deltaTime(mesh_.time().deltaT().value());
    scalar maxDemStep(1.0/minDEMloops_);
    scalar pos(0.0);
    scalar step(maxDemStep);
    scalar historyPos(0.0);
    bool doubleMinStep(false);
    bool doubleMinStep2(false);

    DynamicLabelList ibPrtContactList;

    forAll (prtContactIBList_,pair)
    {
        if (findIndex(ibPrtContactList, prtContactIBList_[pair].first()) == -1)
        {
            ibPrtContactList.append(immersedBodies_[prtContactIBList_[pair].first()].getBodyId());
            immersedBodies_[prtContactIBList_[pair].first()].returnPosition();
        }
        if (findIndex(ibPrtContactList, prtContactIBList_[pair].second()) == -1)
        {
            ibPrtContactList.append(immersedBodies_[prtContactIBList_[pair].second()].getBodyId());
            immersedBodies_[prtContactIBList_[pair].second()].returnPosition();
        }

        if(immersedBodies_[prtContactIBList_[pair].first()].getbodyOperation() == 0)
            immersedBodies_[prtContactIBList_[pair].second()].inContactWithStatic(true);

        if(immersedBodies_[prtContactIBList_[pair].second()].getbodyOperation() == 0)
            immersedBodies_[prtContactIBList_[pair].first()].inContactWithStatic(true);
    }

    forAll (ibPrtContactList,ib)
    {
        label ibIndex(-1);
        ibIndex = findIndex(ibContactList_, ibPrtContactList[ib]);
        if (ibIndex > -1)
        {
            DynamicLabelList helpList;

            for (label i = 0; i < ibIndex; i = i + 1)
            {
                helpList.append(ibContactList_[i]);
            }

            for (label i = ibIndex + 1; i < ibContactList_.size(); i = i + 1)
            {
                helpList.append(ibContactList_[i]);
            }
            ibContactList_ = helpList;
        }
    }

    // resolve prt-prt contacts
    while(prtContactIBList_.size() > 0)
    {
        pos = 0.0;
        step = maxDemStep;
        historyPos = 0.0;
        doubleMinStep = false;
        doubleMinStep2 = false;

        DynamicList<Tuple2<label, label>> pairsToResolve;
        // we have to resolve whole contact chain not only two prts in contact
        findIndexesOfPairWithSomeIb(prtContactIBList_,prtContactIBList_[0], pairsToResolve);

        DynamicLabelList ibToResolve;
        forAll (pairsToResolve,pair)
        {
            if (findIndex(ibToResolve, pairsToResolve[pair].first()) == -1)
            {
                ibToResolve.append(immersedBodies_[pairsToResolve[pair].first()].getBodyId());
            }
            if (findIndex(ibToResolve, pairsToResolve[pair].second()) == -1)
            {
                ibToResolve.append(immersedBodies_[pairsToResolve[pair].second()].getBodyId());
            }
        }

        // return IBs to the start position
        forAll (ibToResolve, ib)
        {
            immersedBodies_[ibToResolve[ib]].initializeVarHistory(true);
            immersedBodies_[ibToResolve[ib]].resetBody(body, false);
            immersedBodies_[ibToResolve[ib]].createImmersedBody(body, refineF);
        }

        // compute motion over time. particles position is reversed when bad motion occurs
        while(true)
        {
            Info << " Start DEM pos: " << pos << " DEM step: " << step << endl;
            forAll (ibToResolve, ib)
            {
                // set F_ and T_ to zero. Do not assign history values
                immersedBodies_[ibToResolve[ib]].initializeVarHistory(false);
                // add fluid coupling force to F and T. This is still same for whole DEM inner loop
                immersedBodies_[ibToResolve[ib]].updateFAndT(immersedBodies_[ibToResolve[ib]].getHistoryCouplingF(), immersedBodies_[ibToResolve[ib]].getHistoryCouplingT());
            }
            // find prt-prt contact info
            forAll (pairsToResolve, pair)
            {
                label cInd(pairsToResolve[pair].first());
                label tInd(pairsToResolve[pair].second());

                Tuple2<Tuple2<vector,vector>,Tuple2<vector,vector>> outVars;

                solvePrtContact(
                    mesh_,
                    immersedBodies_[cInd].getContactInfo(),
                    immersedBodies_[cInd].getContactVars(),
                    immersedBodies_[tInd].getContactInfo(),
                    immersedBodies_[tInd].getContactVars(),
                    geometricD_,
                    deltaTime*step,
                    outVars
                );

                immersedBodies_[cInd].updateFAndT(outVars.first().first(),outVars.first().second());
                immersedBodies_[tInd].updateFAndT(outVars.second().first(),outVars.second().second());
            }

            DynamicList<bool> contactOk;
            forAll (ibToResolve,ib)
            {
                // detect wall contact and solve it
                // update movement and move bodies
                if(immersedBodies_[ibToResolve[ib]].shouldDetectWallContact())
                    detectWallContact(mesh_,immersedBodies_[ibToResolve[ib]].getContactInfo());
                if (immersedBodies_[ibToResolve[ib]].checkWallContact())
                {
                    Info << "-- Body " << immersedBodies_[ibToResolve[ib]].getBodyId() << " is in contact with wall" << endl;
                    Tuple2<vector,vector> outVars;
                    solveWallContact(
                        mesh_,
                        wallInfo_(),
                        immersedBodies_[ibToResolve[ib]].getContactInfo(),
                        immersedBodies_[ibToResolve[ib]].getContactVars(),
                        geometricD_,
                        deltaTime*step,
                        outVars
                        );
                    immersedBodies_[ibToResolve[ib]].updateFAndT(outVars.first(),outVars.second());
                }

                contactOk.append(immersedBodies_[ibToResolve[ib]].checkContactMovement(deltaTime*step));
            }
            // if there is any bad motion reverse all currently moving particles
            bool contactIsOk(true);
            forAll (contactOk,isOK)
            {
                if (!contactOk[isOK])
                {
                    contactIsOk = false;
                    break;
                }
            }

            reduce(contactIsOk, orOp<bool>());

            // move or reverse particles
            forAll (ibToResolve,ib)
            {
                if (contactIsOk || step == minDEMtimeStep_)
                {
                    immersedBodies_[ibToResolve[ib]].assignFullHistory();
                    immersedBodies_[ibToResolve[ib]].updateMovement(deltaTime*step);
                    immersedBodies_[ibToResolve[ib]].resetBody(body);
                    immersedBodies_[ibToResolve[ib]].moveImmersedBody(deltaTime*step);
                    immersedBodies_[ibToResolve[ib]].createImmersedBody(body, refineF);
                }
                else
                {
                    immersedBodies_[ibToResolve[ib]].initializeVarHistory(true);
                    immersedBodies_[ibToResolve[ib]].returnPosition();//return to history position
                    immersedBodies_[ibToResolve[ib]].resetBody(body, false);
                    immersedBodies_[ibToResolve[ib]].createImmersedBody(body, refineF);
                }
            }
            // change dem step. If Step is already minimal wait to avoid infinite loop
            if (doubleMinStep)
            {
                if (doubleMinStep2)
                {
                    doubleMinStep2 = false;
                    doubleMinStep = false;
                }
                else
                    doubleMinStep2 = true;
                historyPos = pos;
                pos += step;
            }
            else if (contactIsOk)
            {
                historyPos = pos;
                pos += step;
                step *= 1.2;
                if (step > maxDemStep)
                    step = maxDemStep;
            }
            else if (!contactIsOk && step == minDEMtimeStep_)
            {
                historyPos = pos;
                pos += step;
            }
            else
            {
                pos = historyPos;
                step /= 2;
                if (step < minDEMtimeStep_)
                {
                    step = minDEMtimeStep_;
                    doubleMinStep = true;
                }
            }

            if (pos + step > 1)
                step = 1 - pos;
            // if moved to end time remove from list or terminate loop
            if (pos >= 1)
            {
                forAll (pairsToResolve, pair)
                {
                    label indexOfPair(-1);
                    indexOfPair = findIndOfPairInNeighbourList(prtContactIBList_, pairsToResolve[pair]);
                    // remove at index
                    DynamicList<Tuple2<label, label>> helpList;

                    for (label i = 0; i < indexOfPair; i = i + 1)
                    {
                        helpList.append(prtContactIBList_[i]);
                    }

                    for (label i = indexOfPair + 1; i < prtContactIBList_.size(); i = i + 1)
                    {
                        helpList.append(prtContactIBList_[i]);
                    }
                    prtContactIBList_ = helpList;
                }

                break;
            }
        }
    }
    // resolve wall contacts
    while(ibContactList_.size() > 0)
    {
        pos = 0.0;
        step = maxDemStep;
        historyPos = 0.0;
        doubleMinStep = false;
        doubleMinStep2 = false;
        immersedBodies_[ibContactList_[0]].initializeVarHistory(true);
        immersedBodies_[ibContactList_[0]].resetBody(body, false);
        immersedBodies_[ibContactList_[0]].createImmersedBody(body, refineF);
        while(true)
        {
            Info << " Start wall DEM pos: " << pos << " DEM step: " << step << endl;
            // set F_ and T_ to zero. Do not assign history values
            immersedBodies_[ibContactList_[0]].initializeVarHistory(false);
            // add fluid coupling force to F and T. This is still same for whole DEM inner loop
            immersedBodies_[ibContactList_[0]].updateFAndT(immersedBodies_[ibContactList_[0]].getHistoryCouplingF(), immersedBodies_[ibContactList_[0]].getHistoryCouplingT());

            // detect wall contact and solve it
            // update movement and move bodies
            if(immersedBodies_[ibContactList_[0]].shouldDetectWallContact())
                    detectWallContact(mesh_,immersedBodies_[ibContactList_[0]].getContactInfo());
            if (immersedBodies_[ibContactList_[0]].checkWallContact())
            {
                Info << "-- Body " << immersedBodies_[ibContactList_[0]].getBodyId() << " is in contact with wall" << endl;
                Tuple2<vector,vector> outVars;
                solveWallContact(
                    mesh_,
                    wallInfo_(),
                    immersedBodies_[ibContactList_[0]].getContactInfo(),
                    immersedBodies_[ibContactList_[0]].getContactVars(),
                    geometricD_,
                    deltaTime*step,
                    outVars
                    );
                immersedBodies_[ibContactList_[0]].updateFAndT(outVars.first(),outVars.second());
            }

            bool contactIsOk(immersedBodies_[ibContactList_[0]].checkContactMovement(deltaTime*step));

            reduce(contactIsOk, orOp<bool>());

            if (contactIsOk || step == minDEMtimeStep_)
            {
                immersedBodies_[ibContactList_[0]].assignFullHistory();
                immersedBodies_[ibContactList_[0]].updateMovement(deltaTime*step);
                immersedBodies_[ibContactList_[0]].resetBody(body);
                immersedBodies_[ibContactList_[0]].moveImmersedBody(deltaTime*step);
                immersedBodies_[ibContactList_[0]].createImmersedBody(body, refineF);
            }
            else
            {
                immersedBodies_[ibContactList_[0]].initializeVarHistory(true);
                immersedBodies_[ibContactList_[0]].returnPosition();
                immersedBodies_[ibContactList_[0]].resetBody(body, false);
                immersedBodies_[ibContactList_[0]].createImmersedBody(body, refineF);
            }

            if (doubleMinStep)
            {
                if (doubleMinStep2)
                {
                    doubleMinStep2 = false;
                    doubleMinStep = false;
                }
                else
                    doubleMinStep2 = true;
                historyPos = pos;
                pos += step;
            }
            else if (contactIsOk)
            {
                historyPos = pos;
                pos += step;
                step *= 1.2;
                if (step > maxDemStep)
                    step = maxDemStep;
            }
            else if (!contactIsOk && step == minDEMtimeStep_)
            {
                historyPos = pos;
                pos += step;
            }
            else
            {
                pos = historyPos;
                step /= 2;
                if (step < minDEMtimeStep_)
                {
                    step = minDEMtimeStep_;
                    doubleMinStep = true;
                }
            }
            if (pos + step > 1)
                step = 1 - pos;

            if (pos >= 1)
            {
                //remove at idex
                DynamicLabelList helpList;

                for (label i = 1; i < ibContactList_.size(); i = i + 1)
                {
                    helpList.append(ibContactList_[i]);
                }
                ibContactList_ = helpList;
                break;
            }
        }
    }

    const dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/pow(average(mesh_.V()), 1.0/3.0)
    );
    surfNorm_ = -fvc::grad(body);
    surfNorm_ /= (mag(surfNorm_)+deltaN.value());

    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].chceckBodyOp();
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::getContactListAndReturnPositions(volScalarField& body)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            if (findIndex(ibContactList_, bodyId) == -1)
            {
                // detect wall contact and assign the body to contact list and return its position
                if(immersedBodies_[bodyId].shouldDetectWallContact())
                    detectWallContact(mesh_,immersedBodies_[bodyId].getContactInfo());
                if (immersedBodies_[bodyId].checkWallContact())
                {
                    ibContactList_.append(immersedBodies_[bodyId].getBodyId());
                    immersedBodies_[bodyId].returnPosition();
                }
            }
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::sortBoundingListPrtContact()
{
    // use insertion sort to sort the neighbour list
    // besides the first sorting, this is O(N) efficient
    for (label coord = 0; coord < 3; coord = coord + 1)
    {
        label i(1);
        while(i < boundValueNeighbourList_[coord].size())
        {
            label j(i);
            while(j > 0 && boundValueNeighbourList_[coord][j-1] > boundValueNeighbourList_[coord][j])
            {
                // if we should swap two bounding points resolve this action
                swapInBoundingListPrtContact(coord, j);
                j = j -1;
            }
            i = i + 1;
        }
    }
    forAll (possibleContactNeighbourList_,value)
    {
        Info << "possible contact: " << possibleContactNeighbourList_[value].first() << " - " << possibleContactNeighbourList_[value].second() << endl;
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::swapInBoundingListPrtContact(label coord, label j)
{
    // swap actual values of bounding points
    scalar scratchScalar(boundValueNeighbourList_[coord][j-1]);
    boundValueNeighbourList_[coord][j-1] = boundValueNeighbourList_[coord][j];
    boundValueNeighbourList_[coord][j] = scratchScalar;

    // update identification list for both bodies
    if (boundLabelNeighbourList_[coord][j-1] > 0)
    {
        immersedBodies_[boundLabelNeighbourList_[coord][j-1] - 1].getBoundIndList()[coord][1] = j;
    }
    else
    {
        immersedBodies_[-1 * boundLabelNeighbourList_[coord][j-1] - 1].getBoundIndList()[coord][0] = j;
    }

    if (boundLabelNeighbourList_[coord][j] > 0)
    {
        immersedBodies_[boundLabelNeighbourList_[coord][j] - 1].getBoundIndList()[coord][1] = j-1;
    }
    else
    {
        immersedBodies_[-1 * boundLabelNeighbourList_[coord][j] - 1].getBoundIndList()[coord][0] = j-1;
    }

    if (boundLabelNeighbourList_[coord][j-1] * boundLabelNeighbourList_[coord][j] < 0)
    {
        // if maximal point had lower index there is new intersection
        if (boundLabelNeighbourList_[coord][j-1] > 0)
        {
            //Create a Tuple2 for bodies that are newly intersected always keep first the body with lower bodyID
            Tuple2<label, label> newPair(min(boundLabelNeighbourList_[coord][j-1] - 1,-1*boundLabelNeighbourList_[coord][j] - 1), max(boundLabelNeighbourList_[coord][j-1] - 1,-1*boundLabelNeighbourList_[coord][j] - 1));
            // if this intersection is not already in possible contact list for this coord assign it
            // Note: This should never happen but it prevent multiple assignment
            if (findIndOfPairInNeighbourList(contactInCoordNeighbourList_[coord], newPair) == -1)
            {
                contactInCoordNeighbourList_[coord].append(newPair);

                bool appendToPossibleContactList(true);
                // check if this pair is already assigned for the remaining dimensions
                // if so, add this pair to possible contact list because their bounding boxes are intersected
                for (label coordi = 0; coordi < 3; coordi = coordi + 1)
                {
                    if (coordi != coord)
                    {
                        if (findIndOfPairInNeighbourList(contactInCoordNeighbourList_[coordi],newPair) == -1)
                        {
                            appendToPossibleContactList = false;
                        }
                    }
                }

                if (appendToPossibleContactList)
                {
                    possibleContactNeighbourList_.append(newPair);
                }
            }
        }
        else
        {
            // there is no intersection for this coord any more for this situation so the pair should be removed from lists
            Tuple2<label, label> removePair(min(-1*boundLabelNeighbourList_[coord][j-1] - 1,boundLabelNeighbourList_[coord][j] - 1), max(-1*boundLabelNeighbourList_[coord][j-1] - 1,boundLabelNeighbourList_[coord][j] - 1));

            label indexOfpair(-1);
            // remove the pair from coord contact list
            indexOfpair = findIndOfPairInNeighbourList(contactInCoordNeighbourList_[coord], removePair);
            if (indexOfpair != -1)
            {
                DynamicList<Tuple2<label, label>> newContactInCoordNeighbourList(0);

                for (label i = 0; i < indexOfpair; i = i + 1)
                {
                    newContactInCoordNeighbourList.append(contactInCoordNeighbourList_[coord][i]);
                }

                for (label i = indexOfpair + 1; i < contactInCoordNeighbourList_[coord].size(); i = i + 1)
                {
                    newContactInCoordNeighbourList.append(contactInCoordNeighbourList_[coord][i]);
                }
                contactInCoordNeighbourList_[coord] = newContactInCoordNeighbourList;
            }

            // same situation for contact list. If there is not instersection in one dimension the
            // bounding boxes are not intersected so remove this pair from contact list
            indexOfpair = findIndOfPairInNeighbourList(possibleContactNeighbourList_, removePair);
            if (indexOfpair != -1)
            {
                DynamicList<Tuple2<label, label>> newPossibleContactNeighbourList(0);

                for (label i = 0; i < indexOfpair; i = i + 1)
                {
                    newPossibleContactNeighbourList.append(possibleContactNeighbourList_[i]);
                }

                for (label i = indexOfpair + 1; i < possibleContactNeighbourList_.size(); i = i + 1)
                {
                    newPossibleContactNeighbourList.append(possibleContactNeighbourList_[i]);
                }

                possibleContactNeighbourList_ = newPossibleContactNeighbourList;
            }
        }
    }

    // actually swap the labels
    label scratchLabel(boundLabelNeighbourList_[coord][j-1]);
    boundLabelNeighbourList_[coord][j-1] = boundLabelNeighbourList_[coord][j];
    boundLabelNeighbourList_[coord][j] = scratchLabel;
}
//---------------------------------------------------------------------------//
// find index of pair in given list. Return -1 if the pair is not in the list
label openHFDIBDEM::findIndOfPairInNeighbourList(DynamicList<Tuple2<label, label>>& listToSearch, Tuple2<label, label> pair)
{
    label returnValue(-1);

    for (label i = 0; i < listToSearch.size(); i = i +1)
    {
        if (listToSearch[i].first() == pair.first() && listToSearch[i].second() == pair.second())
        {
            returnValue = i;
            break;
        }
    }

    return returnValue;
}
//---------------------------------------------------------------------------//
// find index of pairs that contains at least one ib from given pair
void openHFDIBDEM::findIndexesOfPairWithSomeIb(DynamicList<Tuple2<label, label>>& listToSearch, Tuple2<label, label> pair, DynamicList<Tuple2<label, label>>& listToAppend)
{
    DynamicList<Tuple2<label, label>> helpList;

    for (label i = 0; i < listToSearch.size(); i = i +1)
    {
        if (listToSearch[i].first() == pair.first() || listToSearch[i].first() == pair.second() || listToSearch[i].second() == pair.first() || listToSearch[i].second() == pair.second())
        {
            if (findIndOfPairInNeighbourList(listToAppend,listToSearch[i]) == -1)
            {
                listToAppend.append(listToSearch[i]);
                helpList.append(listToSearch[i]);
            }
        }
    }

    for (label i = 0; i < helpList.size(); i++)
    {
        findIndexesOfPairWithSomeIb(listToSearch,helpList[i],listToAppend);
    }
}
//---------------------------------------------------------------------------//
// based on identification list of Bodies update the values of minimal and maximal bounding points
// and sort the neighbour list
void openHFDIBDEM::updateNeighbourLists()
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            List<List<label>> IBboundList = immersedBodies_[bodyId].getBoundIndList();
            vector IBboundMinPoint = immersedBodies_[bodyId].getMinBoundPoint();
            vector IBboundMaxPoint = immersedBodies_[bodyId].getMaxBoundPoint();

            for (label coord = 0; coord < 3; coord = coord + 1)
            {
                boundValueNeighbourList_[coord][IBboundList[coord][0]] = IBboundMinPoint[coord];
                boundValueNeighbourList_[coord][IBboundList[coord][1]] = IBboundMaxPoint[coord];
            }
        }
    }

    sortBoundingListPrtContact();
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

        label maxAdditions(1000);
        label cAddition(0);

        while (addModels_[modelI].shouldAddBody(body) and cAddition < maxAdditions)
        {
            Info << "addModel invoked action, trying to add new body" << endl;
            autoPtr<geomModel> bodyGeomModel(addModels_[modelI].addBody(body));

            cAddition++;

            if (addModels_[modelI].getBodyAdded())
            {
                Info << "STL file correctly generated, registering the new body" << endl;

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
                        geometricD_,
                        bodyGeomModel.ptr(),
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
                nBody.assignFullHistory();

                // update the contact stuff
                for (label i = 0; i < 3; i = i + 1)
                {
                    boundValueNeighbourList_[i].setSize(2 * newIBSize);
                    boundLabelNeighbourList_[i].setSize(2 * newIBSize);
                }

                List<List<label>> IBboundList = nBody.getBoundIndList();
                vector IBboundMinPoint = nBody.getMinBoundPoint();
                vector IBboundMaxPoint = nBody.getMaxBoundPoint();
                label bodyIdInc(addIBPos+1);

                for (label coord = 0; coord < 3; coord = coord + 1)
                {
                    boundValueNeighbourList_[coord][IBboundList[coord][0]] = IBboundMinPoint[coord];
                    boundLabelNeighbourList_[coord][IBboundList[coord][0]] = -1 * bodyIdInc;
                    boundValueNeighbourList_[coord][IBboundList[coord][1]] = IBboundMaxPoint[coord];
                    boundLabelNeighbourList_[coord][IBboundList[coord][1]] = bodyIdInc;
                }

                Info << "new body included into the simulation" << endl;
                cAddition = 0;
            }
            else
            {
                Info << "new body should have been added but was not "
                     << "(probably overlap with an existing body)"
                     << endl;
            }
        }
    }

    sortBoundingListPrtContact();
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::detectPrtContact()
{
    // check only pairs whose bounding boxes are intersected for the contact
    forAll (possibleContactNeighbourList_,possiblePair)
    {
        // check only if the pair is not alredy in contactList
        if (findIndOfPairInNeighbourList(prtContactIBList_, possibleContactNeighbourList_[possiblePair]) == -1)
        {
            label cInd(possibleContactNeighbourList_[possiblePair].first());
            contactInfo& cInfo(immersedBodies_[cInd].getContactInfo());

            label tInd(possibleContactNeighbourList_[possiblePair].second());
            contactInfo& tInfo(immersedBodies_[tInd].getContactInfo());

            if(detectPrtPrtContact(
                mesh_,
                cInfo,
                tInfo
            ))
            {
                prtContactIBList_.append(possibleContactNeighbourList_[possiblePair]);
            };
        }
    }
}
//---------------------------------------------------------------------------//
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
            Info << "Found bodyGeom for " << bodyName << ", the body is: " << bodyGeom << endl;
        }
        else
        {
            bodyGeom = "convex";
            Info << "Did not find bodyGeom for " << bodyName << ", using bodyGeom: " << bodyGeom << endl;
        }

        if(bodyGeom == "convex")
        {
            word stlPath(timePath + "/stlFiles/"+bodyId+".stl");
            bodyGeomModel.set(new convexBody(mesh_,stlPath,thrSurf,geometricD_));
        }
        else if(bodyGeom == "nonConvex")
        {
            word stlPath(timePath + "/stlFiles/"+bodyId+".stl");
            bodyGeomModel.set(new nonConvexBody(mesh_,stlPath,thrSurf,geometricD_));
        }
        else if(bodyGeom == "sphere")
        {
            vector startPosition = bodyDict.subDict("sphere").lookup("position");
            scalar radius = readScalar(bodyDict.subDict("sphere").lookup("radius"));

            bodyGeomModel.set(new sphereBody(mesh_,startPosition,radius,thrSurf,geometricD_));
        }
        else
        {
            word stlPath(timePath + "/stlFiles/"+bodyId+".stl");
            Info << "bodyGeom: " << bodyGeom << " not supported, using bodyGeom nonConvex" << endl;
            bodyGeom = "nonConvex";
            bodyGeomModel.set(new nonConvexBody(mesh_,stlPath,thrSurf,geometricD_));
        }

        label newIBSize(immersedBodies_.size()+1);
        label addIBPos(newIBSize - 1);
        immersedBodies_.setSize(newIBSize);

        Info << "Restarting body: " << bodyId << " as " << addIBPos << " bodyName: " << bodyName << endl;
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
                geometricD_,
                bodyGeomModel.ptr(),
                cellPoints_
            )
        );

        immersedBodies_[addIBPos].createImmersedBody(body,refineF);
        immersedBodies_[addIBPos].computeBodyCharPars();
        immersedBodies_[addIBPos].setRestartSim(Vel,omega,Axis,isStatic,timeStepsInContWStatic);
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::preCalculateCellPoints()
{
    cellPoints_.clear();
    cellPoints_.setSize(mesh_.nCells());
    const pointField& pp = mesh_.points();
    forAll(mesh_.C(), cellI)
    {
        labelList vertexLabels = mesh_.cellPoints()[cellI];
        cellPoints_[cellI] = filterField(pp,vertexLabels);
    }

    forAll (immersedBodies_,bodyId)
    {
        immersedBodies_[bodyId].getGeomModel().resetHashTable();
    }
}
//---------------------------------------------------------------------------//
