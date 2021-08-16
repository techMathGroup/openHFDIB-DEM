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
Description
    class for immersed bodies representation.
SourceFiles
    immersedBodies.C
Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Šourek (2019-*)
\*---------------------------------------------------------------------------*/
#include "immersedBody.H"
#include "polyMesh.H"
#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"

#include "interpolationCellPoint.H"
#include "interpolationCell.H"
#include "meshSearch.H"
#include "List.H"
#include "ListOps.H"

#include "OFstream.H"

#include "cellSet.H"
#include "addToRunTimeSelectionTable.H"
#include "PstreamReduceOps.H"

#include "fvcSmooth.H"
#include "fvMeshSubset.H"

#define ORDER 2

using namespace Foam;

//---------------------------------------------------------------------------//
immersedBody::immersedBody
(
    word bodyName,
    const Foam::dynamicFvMesh& mesh,
    dictionary& HFDIBDEMDict,
    dictionary& transportProperties,
    label bodyId,
    label recomputeM0,
    vector geometricD,
    geomModel* bodyGeomModel,
    List<pointField>& cellPoints
)
:
bodyName_(bodyName),
isActive_(true),
immersedDict_(HFDIBDEMDict.subDict(bodyName_)),
mesh_(mesh),
transportProperties_(transportProperties),
geomModel_(bodyGeomModel),
cellPoints_(cellPoints),
Axis_(vector::one),
AxisOld_(vector::one),
omega_(0.0),
omegaOld_(0.0),
Vel_(vector::zero),
VelOld_(vector::zero),
a_(vector::zero),
alpha_(vector::zero),
totalAngle_(vector::zero),
F_(vector::zero),
T_(vector::zero),
CoNum_(0.0),
rhoF_(transportProperties_.lookup("rho")),
bodyId_(bodyId),
updateTorque_(false),
bodyOperation_(0),
maxDistInDEMloop_(readScalar(HFDIBDEMDict.lookup("maxDistInDEMloop"))),
boundIndList_(3),
owner_(false),
historyCouplingF_(vector::zero),
historyCouplingT_(vector::zero),
historyAxis_(vector::one/mag(vector::one)),
historyOmega_(0.0),
historyVel_(vector::zero),
historya_(vector::zero),
historyAlpha_(vector::zero),
historyTotalAngle_(vector::zero),
octreeField_(mesh_.nCells(), 0),
cellToStartInCreateIB_(0),
totRotMatrix_(tensor::I),
sdBasedLambda_(false),
intSpan_(2.0),
charCellSize_(1e3),
refineBuffers_(0),
useInterpolation_(true),
recomputeM0_(recomputeM0),
geometricD_(geometricD),
created_(false),
timesToSetStatic_(-1)
{
    initializeIB();
}
//---------------------------------------------------------------------------//
immersedBody::~immersedBody()
{
}
//---------------------------------------------------------------------------//
void immersedBody::initializeIB()
{
    #include "initializeIB.H"

    Info << "Finished body initialization" << endl;
    Info << "New bodyID: " << bodyId_ << " name: " << bodyName_ << " rhoS: " << geomModel_->getRhoS() << " dC: " << getDC() << endl;
}
//---------------------------------------------------------------------------//
// Create immersed body info
void immersedBody::createImmersedBody(volScalarField& body, volScalarField& refineF, bool synchCreation)
{
    if(!created_ || bodyOperation_ != 0)
    {
        geomModel_->createImmersedBody(body,octreeField_,surfCells_,intCells_,cellPoints_);
    }

    if(synchCreation)
        syncCreateImmersedBody(body, refineF);
}
//---------------------------------------------------------------------------//
void immersedBody::syncCreateImmersedBody(volScalarField& body, volScalarField& refineF)
{
    if(!created_ || bodyOperation_ != 0)
    {
        owner_ = geomModel_->getOwner();
        Info << "body: " << bodyId_ << " owner: " << owner_ << endl;

        //refine body as stated in the dictionary
        if(useInterpolation_)
        {
            Info << "Computing interpolation points" << endl;
            calculateInterpolationPoints(body);
        }
        Info << "Computing geometrical properties" << endl;
        geomModel_->calculateGeometricalProperties(body,surfCells_,intCells_);
        Info << "-- body " << bodyId_ << " current center of mass position: " << getCoM() << endl;
        created_ = true;
    }

    DynamicLabelList zeroList(surfCells_[Pstream::myProcNo()].size(), 0);
    constructRefineField(body, refineF, surfCells_[Pstream::myProcNo()], zeroList);
    scalarList charCellSizeL(Pstream::nProcs(),1e4);
    forAll (surfCells_[Pstream::myProcNo()],sCellI)
    {
        label cellI = surfCells_[Pstream::myProcNo()][sCellI];
        charCellSizeL[Pstream::myProcNo()] = min(charCellSizeL[Pstream::myProcNo()],Foam::pow(mesh_.V()[cellI],0.3333));
    }
    forAll(charCellSizeL,indl)
    {
        if(charCellSizeL[indl] > 5e3)
            charCellSizeL[indl] = -1.0;
    }

    charCellSize_ = gMax(charCellSizeL);
    Info << "Body characteristic cell size: " << charCellSize_ << endl;
}
//---------------------------------------------------------------------------//
void immersedBody::constructRefineField
(
    volScalarField& body,
    volScalarField& refineF,
    DynamicLabelList cellsToIterate,
    DynamicLabelList startLevel
)
{
    if(refineBuffers_ == 0)
        return;

    DynamicLabelList cellsToIterateC;
    DynamicLabelList cellsToIterateF;

    List<DynamicLabelList> cellsToSendToProcs;
    cellsToSendToProcs.setSize(Pstream::nProcs());
    List<DynamicLabelList> cellsToSendToProcsLevel;
    cellsToSendToProcsLevel.setSize(Pstream::nProcs());

    for(label i = 0; i < refineBuffers_; i++)
    {
        forAll(cellsToIterate, cellI)
        {
            if(startLevel[cellI] == i)
                cellsToIterateC.append(cellsToIterate[cellI]);
        }

        forAll(cellsToIterateC, cellI)
        {
            labelList cellFaces(mesh_.cells()[cellsToIterateC[cellI]]);
            forAll(cellFaces, faceI)
            {
                if (mesh_.isInternalFace(cellFaces[faceI]))
                {
                    label nCell(mesh_.owner()[cellFaces[faceI]]);
                    if(nCell == cellsToIterateC[cellI])
                    {
                        nCell = mesh_.neighbour()[cellFaces[faceI]];
                    }

                    if(refineF[nCell] == 0)
                    {
                        if(i > 0)
                        {
                            if(body[nCell] < SMALL)
                            {
                                refineF[nCell] = 1;
                                cellsToIterateF.append(nCell);
                            }
                        }
                        else
                        {
                            refineF[nCell] = 1;
                            cellsToIterateF.append(nCell);
                        }
                    }
                }
                else
                {
                    label facePatchId(mesh_.boundaryMesh().whichPatch(cellFaces[faceI]));
                    const polyPatch& cPatch = mesh_.boundaryMesh()[facePatchId];
                    if (cPatch.type() == "processor")
                    {
                        const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(cPatch);
                        if (procPatch.myProcNo() == Pstream::myProcNo())
                        {
                            cellsToSendToProcs[procPatch.neighbProcNo()].append(cPatch.whichFace(cellFaces[faceI]));
                            cellsToSendToProcsLevel[procPatch.neighbProcNo()].append(i+1);
                        }
                        else
                        {
                            cellsToSendToProcs[procPatch.myProcNo()].append(cPatch.whichFace(cellFaces[faceI]));
                            cellsToSendToProcsLevel[procPatch.myProcNo()].append(i+1);
                        }
                    }
                }
            }
        }
        cellsToIterateC = cellsToIterateF;
        cellsToIterateF.clear();
    }

    List<DynamicLabelList> facesReceivedFromProcs;
    List<DynamicLabelList> cellsReceivedFromProcsLevel;
    
    // send points that are not on this proc to other proc
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UOPstream send(proci, pBufs);
            send << cellsToSendToProcs[proci];
        }
    }
    pBufs.finishedSends();
    // recieve points from other procs
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recv(proci, pBufs);
            DynamicLabelList recList (recv);
            facesReceivedFromProcs.append(recList);
        }
        else
        {
            DynamicLabelList recList;
            facesReceivedFromProcs.append(recList);
        }
    }
    pBufs.clear();

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {

            UOPstream send(proci, pBufs);
            send << cellsToSendToProcsLevel[proci];
        }
    }
    pBufs.finishedSends();

    // recieve points from other procs
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recv(proci, pBufs);
            DynamicLabelList recList (recv);
            cellsReceivedFromProcsLevel.append(recList);
        }
        else
        {
            DynamicLabelList recList;
            cellsReceivedFromProcsLevel.append(recList);
        }
    }
    pBufs.clear();


    DynamicLabelList newCellsToIterate;
    DynamicLabelList newCellsToIterateStartLevel;

    // check if some point from other proc is on this processor
    for (label otherProci = 0; otherProci < facesReceivedFromProcs.size(); otherProci++)
    {
        for (label faceI = 0; faceI < facesReceivedFromProcs[otherProci].size(); faceI++)
        {
            label cellProcI(0);
            forAll (mesh_.boundaryMesh(), patchi)
            {
                const polyPatch& cPatch = mesh_.boundaryMesh()[patchi];
                if (cPatch.type() == "processor")
                {
                    const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(cPatch);
                    if (procPatch.myProcNo() == Pstream::myProcNo() && procPatch.neighbProcNo() == otherProci)
                    {
                        cellProcI = mesh_.faceOwner()[cPatch.start() + facesReceivedFromProcs[otherProci][faceI]];
                        if(refineF[cellProcI] == 0)
                        {
                            newCellsToIterate.append(cellProcI);
                            newCellsToIterateStartLevel.append(cellsReceivedFromProcsLevel[otherProci][faceI]);
                        }
                        break;
                    }
                    else if (procPatch.myProcNo() == otherProci && procPatch.neighbProcNo() == Pstream::myProcNo())
                    {
                        cellProcI = mesh_.faceOwner()[cPatch.start() + facesReceivedFromProcs[otherProci][faceI]];
                        if(refineF[cellProcI] == 0)
                        {
                            newCellsToIterate.append(cellProcI);
                            newCellsToIterateStartLevel.append(cellsReceivedFromProcsLevel[otherProci][faceI]);
                        }
                        break;
                    }
                }
            }
        }
    }

    bool contBool(false);
    if(newCellsToIterate.size() > 0)
    {
        contBool = true;
    }

    reduce(contBool, orOp<bool>());

    if(contBool)
    {
        constructRefineField(body, refineF, newCellsToIterate, newCellsToIterateStartLevel);
    }
}
//---------------------------------------------------------------------------//
bool immersedBody::shouldDetectWallContact()
{
    if((getM0()-getM()) < 0)
    {
        contactInfo_->setWallContact(false);
        return false;
    }
    return true;
}
//---------------------------------------------------------------------------//
// update immersed body info (pre-contact, ends with contact detection)
void immersedBody::preContactUpdateImmersedBody
(
    volScalarField& body,
    volVectorField& f
)
{
    F_*=0.0;
    T_*=0.0;
    setWallContact(false);

    updateOldMovementVars();

    // assigned variables for potential contact correction
    historyAxis_ = Axis_;
    historyOmega_ = omega_;
    historyVel_ = Vel_;
    historya_ = a_;
    historyAlpha_ = alpha_;
    historyTotalAngle_ = totalAngle_;
    geomModel_->recordHistory();
}
//---------------------------------------------------------------------------//
// update immersed body info (post-contact, contact forces need to be included)
void immersedBody::postContactUpdateImmersedBody
(
    volScalarField& body,
    volVectorField& f
)
{
    // update Vel_, Axis_ and omega_
    updateCoupling(body,f);
    updateMovement(VelOld_,AxisOld_,omegaOld_);

    // update body courant number
    computeBodyCoNumber();

    Info << "-- body: " << bodyId_ << " current center of mass position: " << getCoM() << endl;
}
void immersedBody::postContactUpdateImmersedBody(scalar deltaT)
{
    // update Vel_, Axis_ and omega_
    updateMovement(deltaT);

    // update body courant number
    computeBodyCoNumber();

    Info << "-- body: " << bodyId_ << " current center of mass position: " << getCoM() << endl;
}
//---------------------------------------------------------------------------//
void immersedBody::updateCoupling
(
    volScalarField& body,
    volVectorField& f
)
{
  const uniformDimensionedVectorField& g =
        mesh_.lookupObject<uniformDimensionedVectorField>("g");

  vector FV(vector::zero);
  vector FG(getM()*(1.0-rhoF_.value()/geomModel_->getRhoS().value())*g.value());
  vector TA(vector::zero);

  // calcualate viscous force and torque
  forAll (surfCells_[Pstream::myProcNo()],sCellI)
  {
     label cellI = surfCells_[Pstream::myProcNo()][sCellI];

     FV -=  f[cellI]*mesh_.V()[cellI];
     TA +=  ((mesh_.C()[cellI] - getCoM())^f[cellI])*mesh_.V()[cellI];
  }
  forAll (intCells_[Pstream::myProcNo()],iCellI)
  {
     label cellI = intCells_[Pstream::myProcNo()][iCellI];

     FV -=  f[cellI]*mesh_.V()[cellI];//viscosity?
     TA +=  ((mesh_.C()[cellI] - getCoM())^f[cellI])*mesh_.V()[cellI];
  }

  reduce(FV, sumOp<vector>());
  reduce(TA, sumOp<vector>());

  FV *= rhoF_.value();
  TA *= rhoF_.value();

  historyCouplingF_ = FV+FG;
  historyCouplingT_ = TA;

  F_ += FV+FG;
  T_ -= TA;

  printForcesAndTorques();
}
//---------------------------------------------------------------------------//
// update movement variables of the body
void immersedBody::updateMovement()
{
    scalar deltaT = mesh_.time().deltaT().value();

    updateMovementComp(deltaT,Vel_,Axis_,omega_,1);
}
void immersedBody::updateMovement
(
    scalar deltaT
)
{
    updateMovementComp(deltaT,Vel_,Axis_,omega_,1);
}
void immersedBody::updateMovement
(
    vector Vel,
    vector Axis,
    scalar omega
)
{
    scalar deltaT = mesh_.time().deltaT().value();
    updateMovementComp(deltaT,Vel,Axis,omega,1);
}
void immersedBody::updateMovement
(
    vector Vel,
    vector Axis,
    scalar omega,
    scalar velRelaxFac
)
{
    scalar deltaT = mesh_.time().deltaT().value();
    updateMovementComp(deltaT,Vel,Axis,omega,velRelaxFac);
}
void immersedBody::updateMovement
(
    scalar deltaT,
    vector Vel,
    vector Axis,
    scalar omega
)
{
    updateMovementComp(deltaT,Vel,Axis,omega,1);
}
void immersedBody::updateMovement
(
    scalar deltaT,
    vector Vel,
    vector Axis,
    scalar omega,
    scalar velRelaxFac
)
{
    updateMovementComp(deltaT,Vel,Axis,omega,velRelaxFac);
}
void immersedBody::updateMovementComp
(
    scalar deltaT,
    vector Vel,
    vector Axis,
    scalar omega,
    scalar velRelaxFac
)
{
    auto updateTranslation = [&]()
    {
        // compute current acceleration (assume constant over timeStep)
        a_  = F_/(getM0()+SMALL);
        Info << "--// Body trans update F: " << F_ << endl;

        // update body linear velocity
        Vel_ = (Vel + deltaT*a_) * velRelaxFac;
    };

    auto updateRotation = [&]()
    {
        // update body angular acceleration
        alpha_ = inv(getI()) & T_;

        // update body angular velocity
        vector Omega((Axis*omega + deltaT*alpha_) * velRelaxFac);

        // split Omega into Axis_ and omega_
        omega_ = mag(Omega);

        if (omega_ < SMALL)
        {
            Axis_ = vector::one;
            if (mesh_.nGeometricD() < 3)
            {
                const vector validDirs = (mesh_.geometricD() + Vector<label>::one)/2;
                Axis_ -= validDirs;
            }
        }
        else
        {
            Axis_ =  Omega/(omega_+SMALL);

            if (mesh_.nGeometricD() < 3)
            {// in 2D, I need to keep only the correct part of the rotation axis
                const vector validDirs = (mesh_.geometricD() + Vector<label>::one)/2;
                Axis_ = cmptMultiply(vector::one-validDirs,Axis_);
            }
        }
        Axis_ /= mag(Axis_);
    };

    auto updateRotationFixedAxis = [&]()
    {
        // update body angular velocity
        vector Omega((Axis*omega + deltaT * ( inv(getI()) & T_ )) * velRelaxFac);

        // split Omega into Axis_ and omega_
        omega_ = mag(Omega);

        vector newAxis = Omega/(omega_+SMALL);
        if ((newAxis & Axis_) < 0) Axis_ *= (-1.0);;
    };

    if (bodyOperation_ == 0 or bodyOperation_ == 3)
    {
        F_ *= 0.0;
        T_ *= 0.0;
        return;
    }
    else if (bodyOperation_ == 1)
    {
        updateRotation();

        F_ *= 0.0;
        T_ *= 0.0;
        return;
    }
    else if (bodyOperation_ == 2)
    {
        updateTranslation();

        F_ *= 0.0;
        T_ *= 0.0;
        return;
    }
    else if (bodyOperation_ == 4)
    {
        updateRotationFixedAxis();

        F_ *= 0.0;
        T_ *= 0.0;
        return;
    }

    updateTranslation();
    if (updateTorque_) updateRotation();

    F_ *= 0.0;
    T_ *= 0.0;

    return;
}

//---------------------------------------------------------------------------//
void immersedBody::calculateInterpolationPoints
(
    volScalarField& body
)
{
    meshSearch search(mesh_);

    // stabilisation for normalisation of the interface normal
    const dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/pow(average(mesh_.V()), 1.0/3.0)
    );
    
    // create temporary unit surface normals
    vectorField surfNorm(-fvc::grad(body));
    surfNorm /= (mag(surfNorm)+deltaN.value());

    // clear the old interpolation data
    interpolationInfo_[Pstream::myProcNo()].clear();
    interpolationVecReqs_.clear();
    interpolationVecReqs_.setSize(Pstream::nProcs());

    // variables to find in other processors
    List<DynamicLabelList> surfCellRef;
    surfCellRef.setSize(Pstream::nProcs());
    List<DynamicLabelList> orderRef;
    orderRef.setSize(Pstream::nProcs());
    List<DynamicLabelList> surfCellPointRef;
    surfCellPointRef.setSize(Pstream::nProcs());
    List<DynamicLabelList> surfCellLabelRef;
    surfCellLabelRef.setSize(Pstream::nProcs());
    List<DynamicPointList> surfCellLabelRefP;
    surfCellLabelRefP.setSize(Pstream::nProcs());

    forAll (surfCells_[Pstream::myProcNo()],cell)
    {
        // get surface cell label
        label scell = surfCells_[Pstream::myProcNo()][cell];

        // estimate intDist
        scalar intDist(0);
        label cellI(scell);
        point surfPoint(mesh_.C()[scell]);

        if (sdBasedLambda_)
        {
            scalar minMaxBody(max(min(body[scell],1.0-SMALL),SMALL));
            intDist = Foam::atanh(2.0*minMaxBody - 1.0)*Foam::pow(mesh_.V()[scell],0.333)/intSpan_;
            surfPoint += surfNorm[scell]*intDist;
        }
        else
        {
            scalar dotProd(-GREAT);
            labelList cellNb(mesh_.cellCells()[scell]);//list of neighbours
            forAll (cellNb,nbCellI)
            {
                vector rVec(mesh_.C()[cellNb[nbCellI]] - mesh_.C()[scell]);
                intDist += mag(rVec);
                scalar auxDotProd(rVec & surfNorm[scell]);
                if (auxDotProd > dotProd)
                {
                    dotProd = auxDotProd;
                    cellI   = cellNb[nbCellI];
                }
            }
            intDist /= (scalar(cellNb.size())+SMALL);
            surfPoint -= intDist*surfNorm[scell]*(0.5-body[scell]);
        }

        // create vector for points and cells and add to main vectors
        DynamicPointList intPoints;
        DynamicLabelList intCells;
        DynamicLabelList procOfCells;

        intDist = Foam::pow(mesh_.V()[scell],0.333);
        intDist*= 0.5;

        // add to list
        intPoints.append(surfPoint);
        vector surfNormToSend(vector::zero);
        if (mag(surfNorm[scell]) > SMALL)
        {
            surfNormToSend = surfNorm[scell]/mag(surfNorm[scell]);
            Tuple2<label,label> helpTup(scell,-1);
            Tuple2<vector,Tuple2<label,label>> startCell(vector::zero,helpTup);
            // add other interpolation points
            for (int order=0;order<ORDER;order++)
            {
                startCell = findCellCustom(startCell.first(),startCell.second().first(),startCell.second().second(),surfNormToSend, intDist);
                surfPoint = startCell.first();
                cellI = startCell.second().first();

                if (startCell.second().second() == -1)
                {
                    if (startCell.second().first() != -1)
                    {
                        if (body[cellI] > SMALL)
                        {
                            order--;
                            continue;
                        }
                    }
                    intPoints.append(surfPoint);
                    intCells.append(cellI);
                    procOfCells.append(Pstream::myProcNo());
                }
                else
                {
                    intPoints.append(surfPoint);
                    intCells.append(-1);
                    procOfCells.append(startCell.second().second());
                    surfCellRef[startCell.second().second()].append(cell);
                    orderRef[startCell.second().second()].append(order);
                    surfCellLabelRef[startCell.second().second()].append(cellI);
                    surfCellPointRef[startCell.second().second()].append(order);
                    surfCellLabelRefP[startCell.second().second()].append(surfPoint);
                }
            }
        }
        else
        {
            for (int order=0;order<ORDER;order++)
            {
                intPoints.append(surfPoint);
                intCells.append(-1);
                procOfCells.append(Pstream::myProcNo());
            }
        }
        
        // assign to global variables
        immersedBody::interpolationInfo intInfo;
        intInfo.surfCell_ = scell;
        intInfo.intPoints_ = intPoints;
        intInfo.intCells_ = intCells;
        intInfo.procWithIntCells_ = procOfCells;
        intInfo.intVec_.setSize(2);
        interpolationInfo_[Pstream::myProcNo()].append(intInfo);
    }
    
    // send points that are not on this proc to other proc
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    PstreamBuffers pBufs2(Pstream::commsTypes::nonBlocking);

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UOPstream send(proci, pBufs);
            send << surfCellLabelRef[proci];

            UOPstream send2(proci, pBufs2);
            send2 << surfCellPointRef[proci];
        }
    }
    pBufs.finishedSends();
    pBufs2.finishedSends();

    List<DynamicLabelList> surfPointRefFromProcs;
    List<DynamicLabelList> surfPointLabelRefFromProcs;
    List<DynamicLabelList> surfCellRefToReturn;
    // recieve points from other procs
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recv(proci, pBufs);
            DynamicLabelList recSurfPointLabelRef (recv);
            surfPointLabelRefFromProcs.append(recSurfPointLabelRef);
            UIPstream recv2(proci, pBufs2);
            DynamicLabelList recSurfPointRef (recv2);
            surfPointRefFromProcs.append(recSurfPointRef);
        }
        else
        {
            DynamicLabelList recSurfPointLabelRef;
            surfPointLabelRefFromProcs.append(recSurfPointLabelRef);
            DynamicLabelList recSurfPointRef;
            surfPointRefFromProcs.append(recSurfPointRef);
        }
    }

    // check if some point from other proc is on this processor
    for (label otherProci = 0; otherProci < surfPointRefFromProcs.size(); otherProci++)
    {
        DynamicLabelList surfCellRefFromProc;
        for (label surfPointI = 0; surfPointI < surfPointRefFromProcs[otherProci].size(); surfPointI++)
        {
            label cellProcI(0);
            forAll (mesh_.boundaryMesh(), patchi)
            {
                const polyPatch& cPatch = mesh_.boundaryMesh()[patchi];
                if (cPatch.type() == "processor")
                {
                    const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(cPatch);
                    if (procPatch.myProcNo() == Pstream::myProcNo() && procPatch.neighbProcNo() == otherProci)
                    {
                        cellProcI = mesh_.faceOwner()[cPatch.start() + surfPointLabelRefFromProcs[otherProci][surfPointI]];
                        if (surfPointRefFromProcs[otherProci][surfPointI] == 1)
                        {
                            vector sfUnit(mesh_.Sf()[cPatch.start() + surfPointLabelRefFromProcs[otherProci][surfPointI]]/mag(mesh_.Sf()[cPatch.start() + surfPointLabelRefFromProcs[otherProci][surfPointI]]));

                            label bestCell(0);

                            scalar dotProd(-GREAT);
                            labelList cellNb(mesh_.cellCells()[cellProcI]);//list of neighbours
                            forAll (cellNb,nbCellI)
                            {
                                vector vecI(mesh_.C()[cellNb[nbCellI]] - mesh_.C()[cellProcI]);
                                vecI /= mag(vecI);
                                scalar auxDotProd(vecI & sfUnit);
                                if (auxDotProd > dotProd)
                                {
                                    dotProd = auxDotProd;
                                    bestCell = cellNb[nbCellI];
                                }
                            }
                            cellProcI = bestCell;
                        }
                        break;
                    }
                    else if (procPatch.myProcNo() == otherProci && procPatch.neighbProcNo() == Pstream::myProcNo())
                    {
                        cellProcI = mesh_.faceOwner()[cPatch.start() + surfPointLabelRefFromProcs[otherProci][surfPointI]];
                        if (surfPointRefFromProcs[otherProci][surfPointI] == 1)
                        {
                            vector sfUnit(mesh_.Sf()[cPatch.start() + surfPointLabelRefFromProcs[otherProci][surfPointI]]/mag(mesh_.Sf()[cPatch.start() + surfPointLabelRefFromProcs[otherProci][surfPointI]]));

                            label bestCell(0);

                            scalar dotProd(-GREAT);
                            labelList cellNb(mesh_.cellCells()[cellProcI]);//list of neighbours
                            forAll (cellNb,nbCellI)
                            {
                                vector vecI(mesh_.C()[cellNb[nbCellI]] - mesh_.C()[cellProcI]);
                                vecI /= mag(vecI);
                                scalar auxDotProd(vecI & sfUnit);
                                if (auxDotProd > dotProd)
                                {
                                    dotProd = auxDotProd;
                                    bestCell = cellNb[nbCellI];
                                }
                            }
                            cellProcI = bestCell;
                        }
                        break;
                    }
                }
            }

            if (cellProcI >= 0)
            {
                if (body[cellProcI] > SMALL)
                    cellProcI = -1;
            }
            surfCellRefFromProc.append(cellProcI);
        }
        // assing found points to send them back
        surfCellRefToReturn.append(surfCellRefFromProc);
    }

    pBufs.clear();
    // send points back
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UOPstream send(proci, pBufs);
            send << surfCellRefToReturn[proci];
        }
    }

    pBufs.finishedSends();

    List<DynamicLabelList> surfCellRefRcv;
    // receive points from other processors
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recv(proci, pBufs);
            DynamicLabelList recsurfCellRef (recv);
            surfCellRefRcv.append(recsurfCellRef);
        }
        else
        {
            DynamicLabelList recsurfCellRef;
            surfCellRefRcv.append(recsurfCellRef);
        }
    }

    // check if some points were found and assign them to info
    for (label otherProci = 0; otherProci < surfCellRefRcv.size(); otherProci++)
    {
        for (label intCellI = 0; intCellI < surfCellRefRcv[otherProci].size(); intCellI++)
        {
            label cellProcI(surfCellRefRcv[otherProci][intCellI]);
            if (cellProcI > -1)
            {
                interpolationInfo_[Pstream::myProcNo()][surfCellRef[otherProci][intCellI]].intCells_[orderRef[otherProci][intCellI]] = cellProcI;
            }
        }
    }
    // decide which order should be used
    for (label infoI = 0; infoI < interpolationInfo_[Pstream::myProcNo()].size(); infoI++)
    {
        List<bool> allowedOrder;
        allowedOrder.setSize(ORDER);

        for (int intPoint=0;intPoint<ORDER;intPoint++)
        {
            if (interpolationInfo_[Pstream::myProcNo()][infoI].intCells_[intPoint] == -1)
            {
                allowedOrder[intPoint] = false;
            }
            else
            {
                allowedOrder[intPoint] = true;
            }
        }

        interpolationInfo_[Pstream::myProcNo()][infoI].order_ = 2;
        if ( allowedOrder[1] == false)
        {
            interpolationInfo_[Pstream::myProcNo()][infoI].order_ = 1;
        }

        // check if first order is possible
        if ( allowedOrder[0] == false)
        {
            interpolationInfo_[Pstream::myProcNo()][infoI].order_ = 0;
        }

        if (interpolationInfo_[Pstream::myProcNo()][infoI].order_ == 2)
        {
            if (interpolationInfo_[Pstream::myProcNo()][infoI].intCells_[0] == interpolationInfo_[Pstream::myProcNo()][infoI].intCells_[1])
                interpolationInfo_[Pstream::myProcNo()][infoI].order_ = 1;
        }
        // if some cell on other proc prepare requst for interpolation
        switch(interpolationInfo_[Pstream::myProcNo()][infoI].order_)
        {
            case 1:
                if (interpolationInfo_[Pstream::myProcNo()][infoI].procWithIntCells_[0] != Pstream::myProcNo())
                {
                    immersedBody::intVecRequest vecReq;
                    vecReq.requestLabel_ = infoI;
                    vecReq.vecLabel_ = 0;
                    vecReq.intPoint_ = interpolationInfo_[Pstream::myProcNo()][infoI].intPoints_[1];
                    vecReq.intCell_ = interpolationInfo_[Pstream::myProcNo()][infoI].intCells_[0];
                    interpolationVecReqs_[interpolationInfo_[Pstream::myProcNo()][infoI].procWithIntCells_[0]].append(vecReq);
                }
                break;
            case 2:
                if (interpolationInfo_[Pstream::myProcNo()][infoI].procWithIntCells_[0] != Pstream::myProcNo())
                {
                    immersedBody::intVecRequest vecReq;
                    vecReq.requestLabel_ = infoI;
                    vecReq.vecLabel_ = 0;
                    vecReq.intPoint_ = interpolationInfo_[Pstream::myProcNo()][infoI].intPoints_[1];
                    vecReq.intCell_ = interpolationInfo_[Pstream::myProcNo()][infoI].intCells_[0];
                    interpolationVecReqs_[interpolationInfo_[Pstream::myProcNo()][infoI].procWithIntCells_[0]].append(vecReq);
                }

                if (interpolationInfo_[Pstream::myProcNo()][infoI].procWithIntCells_[1] != Pstream::myProcNo())
                {
                    immersedBody::intVecRequest vecReq;
                    vecReq.requestLabel_ = infoI;
                    vecReq.vecLabel_ = 1;
                    vecReq.intPoint_ = interpolationInfo_[Pstream::myProcNo()][infoI].intPoints_[2];
                    vecReq.intCell_ = interpolationInfo_[Pstream::myProcNo()][infoI].intCells_[1];
                    interpolationVecReqs_[interpolationInfo_[Pstream::myProcNo()][infoI].procWithIntCells_[1]].append(vecReq);
                }
                break;
        }
    }
}
//---------------------------------------------------------------------------//
// Custom function to find cell containing point
// Note: this is much cheaper than standard OF functions
Tuple2<vector,Tuple2<label,label>> immersedBody::findCellCustom
(
    vector& prevPoint,
    label& startCell,
    label& startProc,
    vector& gradToBody,
    scalar& intDist
)
{
    if (startProc != -1)
    {

        Tuple2<label,label> helpTup(startCell,startProc);
        Tuple2<vector,Tuple2<label,label>> tupleToReturn(prevPoint + intDist*gradToBody,helpTup);
        return tupleToReturn;
    }
    else if (startCell == -1)
    {
        Tuple2<label,label> helpTup(-1,-1);
        Tuple2<vector,Tuple2<label,label>> tupleToReturn(prevPoint + intDist*gradToBody,helpTup);
        return tupleToReturn;
    }

    labelList cellFaces(mesh_.cells()[startCell]);
    label bestFace(0);

    scalar dotProd(-GREAT);
    forAll (cellFaces,faceI)
    {
        vector vecI(mesh_.Cf()[cellFaces[faceI]] - mesh_.C()[startCell]);
        vecI /= mag(vecI);
        scalar auxDotProd(vecI & gradToBody);
        if (auxDotProd > dotProd)
        {
            dotProd = auxDotProd;
            bestFace = cellFaces[faceI];
        }
    }
    labelList cellPoints(mesh_.faces()[bestFace]);
    DynamicList<Tuple2<vector,Tuple2<label,label>>> pointsToCheck;
    forAll (cellPoints, pointI)
    {
        labelList pointFaces(mesh_.pointFaces()[cellPoints[pointI]]);
        forAll (pointFaces, faceI)
        {
            if (mesh_.isInternalFace(pointFaces[faceI]))
            {
                label owner(mesh_.owner()[pointFaces[faceI]]);
                label neighbour(mesh_.neighbour()[pointFaces[faceI]]);
                if (owner != startCell)
                {
                    bool add(true);
                    forAll (pointsToCheck, i)
                    {
                        if (pointsToCheck[i].second().first() == owner)
                        {
                            add = false;
                            break;
                        }
                    }
                    if (add)
                    {
                        Tuple2<label,label> helpTup(owner,-1);
                        Tuple2<vector,Tuple2<label,label>> tupleToAdd(mesh_.C()[owner],helpTup);
                        pointsToCheck.append(tupleToAdd);
                    }
                }
                if (neighbour != startCell)
                {
                    bool add(true);
                    forAll (pointsToCheck, i)
                    {
                        if (pointsToCheck[i].second().first() == neighbour)
                        {
                            add = false;
                            break;
                        }
                    }
                    if (add)
                    {
                        Tuple2<label,label> helpTup(neighbour,-1);
                        Tuple2<vector,Tuple2<label,label>> tupleToAdd(mesh_.C()[neighbour],helpTup);
                        pointsToCheck.append(tupleToAdd);
                    }
                }
            }
            else
            {
                label owner(mesh_.faceOwner()[pointFaces[faceI]]);
                vector distToFace(mesh_.Cf()[pointFaces[faceI]] - mesh_.C()[owner]);
                vector sfUnit(mesh_.Sf()[pointFaces[faceI]]/mag(mesh_.Sf()[pointFaces[faceI]]));
                vector pointToAppend(mesh_.C()[owner] + mag(distToFace)*sfUnit);

                label facePatchId(mesh_.boundaryMesh().whichPatch(pointFaces[faceI]));
                const polyPatch& cPatch = mesh_.boundaryMesh()[facePatchId];
                if (cPatch.type() == "processor")
                {
                    const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(cPatch);
                    if (procPatch.myProcNo() == Pstream::myProcNo())
                    {
                        Tuple2<label,label> helpTup(cPatch.whichFace(pointFaces[faceI]),procPatch.neighbProcNo());
                        Tuple2<vector,Tuple2<label,label>> tupleToAdd(pointToAppend,helpTup);
                        pointsToCheck.append(tupleToAdd);
                    }
                    else
                    {
                        Tuple2<label,label> helpTup(cPatch.whichFace(pointFaces[faceI]),procPatch.myProcNo());
                        Tuple2<vector,Tuple2<label,label>> tupleToAdd(pointToAppend,helpTup);
                        pointsToCheck.append(tupleToAdd);
                    }
                }
                else
                {
                    bool add(true);
                    forAll (pointsToCheck, i)
                    {
                        if (pointsToCheck[i].first() == pointToAppend)
                        {
                            add = false;
                            break;
                        }
                    }
                    if (add)
                    {
                        Tuple2<label,label> helpTup(-1,-1);
                        Tuple2<vector,Tuple2<label,label>> tupleToAdd(pointToAppend,helpTup);
                        pointsToCheck.append(tupleToAdd);
                    }
                }
            }
        }
    }

    dotProd = -GREAT;
    Tuple2<label,label> helpTup(-1,-1);
    Tuple2<vector,Tuple2<label,label>> tupleToReturn(vector::zero,helpTup);
    forAll (pointsToCheck,pointI)
    {
        vector vecI(pointsToCheck[pointI].first() - mesh_.C()[startCell]);
        vecI /= (mag(vecI)+SMALL);
        scalar auxDotProd(vecI & gradToBody);
        if (auxDotProd > dotProd)
        {
            dotProd = auxDotProd;
            tupleToReturn   = pointsToCheck[pointI];
        }
    }

    return tupleToReturn;
}
//---------------------------------------------------------------------------//
// move immersed body according to body operation
void immersedBody::moveImmersedBody
(
    scalar deltaT
)
{
    if (bodyOperation_ == 0) return;

    if (owner_ == Pstream::myProcNo())
    {
        if (mag(deltaT + 1.0) < SMALL) deltaT = mesh_.time().deltaT().value();

        // incremental rotation angle
        scalar angle     = omega_*deltaT - 0.5*mag(alpha_)*deltaT*deltaT;

        // translation increment
        vector transIncr = Vel_*deltaT - 0.5*a_*deltaT*deltaT;
        
        // rotation matrix
        tensor rotMatrix(Foam::cos(angle)*tensor::I);
        rotMatrix += Foam::sin(angle)*tensor(
            0.0,      -Axis_.z(),  Axis_.y(),
            Axis_.z(), 0.0,       -Axis_.x(),
            -Axis_.y(), Axis_.x(),  0.0
        );
        rotMatrix += (1.0-Foam::cos(angle))*(Axis_ * Axis_);

        // update total rotation matrix
        totRotMatrix_ = rotMatrix & totRotMatrix_;
        vector eulerAngles;
        scalar sy = Foam::sqrt(totRotMatrix_.xx()*totRotMatrix_.xx() + totRotMatrix_.yy()*totRotMatrix_.yy());
        if (sy > SMALL)
        {
            eulerAngles.x() = Foam::atan2(totRotMatrix_.zy(),totRotMatrix_.zz());
            eulerAngles.y() = Foam::atan2(-totRotMatrix_.zx(),sy);
            eulerAngles.z() = Foam::atan2(totRotMatrix_.yx(),totRotMatrix_.xx());
        }
        else
        {
            eulerAngles.x() = Foam::atan2(-totRotMatrix_.yz(),totRotMatrix_.yy());
            eulerAngles.y() = Foam::atan2(-totRotMatrix_.zx(),sy);
            eulerAngles.z() = 0.0;
        }

        Info << "-- body " << bodyId_ << " linear velocity      :  " << Vel_ << endl;
        Info << "-- body " << bodyId_ << " angluar velocity     : " << omega_ << endl;
        Info << "-- body " << bodyId_ << " axis of rotation     : " << Axis_ << endl;
        Info << "-- body " << bodyId_ << " total rotation matrix: " << totRotMatrix_ << endl;
        Info << "-- body " << bodyId_ << " total euler angles   : " << eulerAngles << endl;


        geomModel_->bodyRotatePoints(angle,Axis_);
        geomModel_->bodyMovePoints(transIncr);
    }

    geomModel_->synchronPos();

    // update bounds of the body
    boundBox bound(geomModel_->getBounds());
    minBoundPoint_ = bound.min();
    maxBoundPoint_ = bound.max();
}
//---------------------------------------------------------------------------//
void immersedBody::updateVectorField(volVectorField& VS, word VName,volScalarField& body, vectorField surfNorm)
{
    // check dictionary for parameters (only noSlip allowed)
    word BC = immersedDict_.subDict(VName).lookup("BC");

    if (BC=="noSlip")
    {
        // if STATICBODY set to zero
        if ( bodyOperation_==0)
        {
            forAll (surfCells_[Pstream::myProcNo()],cell)
            {
                label cellI = surfCells_[Pstream::myProcNo()][cell];
                VS[cellI]   = Vel_;
            }
            forAll (intCells_[Pstream::myProcNo()],cell)
            {
                label cellI = intCells_[Pstream::myProcNo()][cell];
                VS[cellI] = Vel_;
            }
        }
        else
        {
            label cellI;
            forAll (surfCells_[Pstream::myProcNo()],cell)
            {
                cellI=surfCells_[Pstream::myProcNo()][cell];

                // correction for estimated interface position
                scalar intDist(0);
                point surfPoint(mesh_.C()[cellI]);

                if (sdBasedLambda_)
                {
                    scalar minMaxBody(max(min(body[cellI],1.0-SMALL),SMALL));
                    intDist = Foam::atanh(2.0*minMaxBody - 1.0)*Foam::pow(mesh_.V()[cellI],0.333)/intSpan_;
                    surfPoint += surfNorm[cellI]*intDist;
                }
                else
                {
                    labelList cellNb(mesh_.cellCells()[cellI]);//list of neighbours
                    forAll (cellNb,nbCellI)
                    {
                        vector rVec(mesh_.C()[cellNb[nbCellI]] - mesh_.C()[cellI]);
                        intDist += mag(rVec);
                    }
                    intDist /= (scalar(cellNb.size())+SMALL);
                    surfPoint -= intDist*surfNorm[cellI]*(0.5-body[cellI]);
                }

                vector planarVec       =  surfPoint - getCoM()
                                     - Axis_*(
                                          (surfPoint-getCoM())&Axis_
                                         );

                vector VSvalue = -(planarVec^Axis_)*omega_ + Vel_;
                VS[cellI] = VSvalue;
            }
            forAll (intCells_[Pstream::myProcNo()],cell)
            {
                cellI=intCells_[Pstream::myProcNo()][cell];

                vector planarVec       =  mesh_.C()[cellI] - getCoM()
                                     - Axis_*(
                                          (mesh_.C()[cellI]-getCoM())&Axis_
                                         );

                vector VSvalue = -(planarVec^Axis_)*omega_ + Vel_;
                VS[cellI] = VSvalue;
            }
        }

    }

}
//---------------------------------------------------------------------------//
// reset body field for this immersed object
void immersedBody::resetBody(volScalarField& body, bool resethistoryFt)
{
    if(!created_ || bodyOperation_ != 0)
    {
        forAll (intCells_[Pstream::myProcNo()],cellI)
        {
            body[intCells_[Pstream::myProcNo()][cellI]] = 0;
        }
        forAll (surfCells_[Pstream::myProcNo()],cellI)
        {
            body[surfCells_[Pstream::myProcNo()][cellI]] = 0;
        }
    }

    setWallContact(false);

    if(!created_ || bodyOperation_ != 0)
    {
        interpolationInfo_[Pstream::myProcNo()].clear();
        interpolationVecReqs_[Pstream::myProcNo()].clear();

        surfCells_[Pstream::myProcNo()].clear();
        intCells_[Pstream::myProcNo()].clear();
    }

    //keep last Ft force only if it was updated in this time step
    if (resethistoryFt)
    {
        DynamicList<Tuple2<label,Tuple2<label,vector>>> historyFtNew;

        forAll (historyFt_,Fti)
        {
            if (historyFt_[Fti].second().first())
            {
                Tuple2<label,vector> help(0,historyFt_[Fti].second().second());
                Tuple2<label,Tuple2<label,vector>> help2(historyFt_[Fti].first(), help);

                historyFtNew.append(help2);
            }
        }

        historyFt_ = historyFtNew;
    }
}
//---------------------------------------------------------------------------//
// function to compute local particle radii w.r.t. CoM_
DynamicList<scalar> immersedBody::getLocPartRad(DynamicLabelList& cellsOfInt)
{
    DynamicScalarList aLst;

    forAll (cellsOfInt,cellI)
    {
        label cCell(cellsOfInt[cellI]);
        aLst.append(mag(mesh_.C()[cCell] - getCoM()));
    }
    return aLst;
}
//---------------------------------------------------------------------------//
// function to move the body after the contact
void immersedBody::postContactUpdateBodyField(volScalarField& body, volScalarField& refineF)
{
    // move the body due to the contact
    resetBody(body);
    createImmersedBody(body,refineF,false);
}
//---------------------------------------------------------------------------//
void immersedBody::recreateBodyField(volScalarField& body, volScalarField& refineF)
{
    created_ = false;

    octreeField_ = Field<label>(mesh_.nCells(), 0);
    interpolationInfo_[Pstream::myProcNo()].clear();
    interpolationVecReqs_[Pstream::myProcNo()].clear();
    surfCells_[Pstream::myProcNo()].clear();
    intCells_[Pstream::myProcNo()].clear();
    createImmersedBody(body,refineF,false);
}
//---------------------------------------------------------------------------//
// function to compute maximal and mean courant number of the body
void immersedBody::computeBodyCoNumber()
{
    label auxCntr(0);
    scalar VelMag(mag(Vel_));

    meanCoNum_ = 0.0;

    // rotation body courant number
    scalar rotCoNumB(omega_*getDC()*0.5*mesh_.time().deltaT().value());


    forAll (surfCells_[Pstream::myProcNo()],sCellI)
    {
        label cellI(surfCells_[Pstream::myProcNo()][sCellI]);

        scalar dCell(Foam::pow(mesh_.V()[cellI],0.3333));
        scalar CoNumCell(VelMag*mesh_.time().deltaT().value()/dCell);

        CoNumCell+=rotCoNumB/dCell;

        CoNum_      = max(CoNum_,CoNumCell);
        meanCoNum_ += CoNumCell;
        auxCntr    += 1;
    }
    forAll (intCells_[Pstream::myProcNo()],iCellI)
    {
        label cellI(intCells_[Pstream::myProcNo()][iCellI]);

        scalar dCell(Foam::pow(mesh_.V()[cellI],0.3333));
        scalar CoNumCell(VelMag*mesh_.time().deltaT().value()/dCell);

        CoNum_      = max(CoNum_,CoNumCell);
        meanCoNum_ += CoNumCell;
        auxCntr    += 1;
    }

    reduce(meanCoNum_, sumOp<scalar>());
    reduce(auxCntr, sumOp<scalar>());
    reduce(CoNum_, maxOp<scalar>());

    if(auxCntr > 0)
    {
        meanCoNum_ /= auxCntr;
    }

    Info << "-- body " << bodyId_ << " Courant Number mean: " << meanCoNum_
         << " max: " << CoNum_ << endl;

}

//---------------------------------------------------------------------------//
// print out body linear and angular momentum
void immersedBody::printMomentum()
{
    vector L(getI()&(Axis_*omega_));
    vector p(getM()*Vel_);

    Info << "-- body " << bodyId_ << "  linear momentum:" << p
         << " magnitude: " << mag(p) <<endl;
    Info << "-- body " << bodyId_ << " angular momentum:" << L
         << " magnitude: " << mag(L) <<endl;
}
//---------------------------------------------------------------------------//
// print out body statistics
void immersedBody::printStats()
{
    vector L(getI()&(Axis_*omega_));
    vector p(getM()*Vel_);

    Info << "-- body " << bodyId_ << "  linear momentum:" << p
         << " magnitude: " << mag(p) <<endl;
    Info << "-- body " << bodyId_ << " angular momentum:" << L
         << " magnitude: " << mag(L) <<endl;
    Info << "-- body " << bodyId_ << "  linear velocity:" << Vel_
         << " magnitude: " << mag(Vel_) <<endl;
    Info << "-- body " << bodyId_ << " angular velocity:" << omega_
         << " magnitude: " << mag(omega_) <<endl;
    Info << "-- body " << bodyId_ << "    rotation axis:" << Axis_
         << " magnitude: " << mag(Axis_) <<endl;
}
//---------------------------------------------------------------------------//
// print out eulerian forces acting on the body
void immersedBody::printForcesAndTorques()
{
    const uniformDimensionedVectorField& g =
        mesh_.lookupObject<uniformDimensionedVectorField>("g");
    vector FG(getM()*(1.0-rhoF_.value()/geomModel_->getRhoS().value())*g.value());

    Info << "-- body " << bodyId_ << "     linear force:" << F_
         << " magnitude: " << mag(F_) <<endl;
    Info << "-- body " << bodyId_ << "   grav/buy force:" << FG
         << " magnitude: " << mag(FG) <<endl;
    Info << "-- body " << bodyId_ << "    viscous force:" << F_ - FG
         << " magnitude: " << mag(F_ - FG) <<endl;
    Info << "-- body " << bodyId_ << "           torque:" << T_
         << " magnitude: " << mag(T_) <<endl;
}
//---------------------------------------------------------------------------//
// return to history position when the particle gets in contact during the time step
void immersedBody::returnPosition()
{
    geomModel_->returnHistory();
    boundBox bound(geomModel_->getBounds());
    minBoundPoint_ = bound.min();
    maxBoundPoint_ = bound.max();
}
//---------------------------------------------------------------------------//
// initialize variables base to history. Set history variables only when needed
void immersedBody::initializeVarHistory(bool setHistory)
{
    if (setHistory && bodyOperation_ != 0)
    {
        Axis_ = historyAxis_;
        omega_ = historyOmega_;
        Vel_ = historyVel_;
        a_ = historya_;
        alpha_ = historyAlpha_;
        totalAngle_ = historyTotalAngle_;
    }
    F_*=0.0;
    T_*=0.0;
}
//---------------------------------------------------------------------------//
// switch the particle off (remove it from the simulation)
void immersedBody::switchActiveOff
(
    volScalarField& body
)
{
    // turn of the particle
    isActive_ = false;

    // rewrite the body field
    resetBody(body);
}
//---------------------------------------------------------------------------//
//update movement variables of the body
bool immersedBody::checkContactMovement
(
    scalar deltaT
)
{
    if (mag(deltaT + 1.0) < SMALL) deltaT = mesh_.time().deltaT().value();

    vector ai(vector::zero);
    vector Veli(Vel_);
    vector Axisi = Axis_;
    scalar omegai = omega_;
    vector alphai = alpha_;

    // auxiliary (nested) functions
    auto updateTranslation = [&]()
    {
        // compute current acceleration (assume constant over timeStep)
        ai  = F_/(getM0()+SMALL);

        // update body linear velocity
        Veli += deltaT*ai;
    };

    auto updateRotation = [&]()
    {
        // update body angular acceleration
        alphai = inv(getI()) & T_;

        // update body angular velocity
        vector Omega(Axis_*omegai + deltaT*alphai);

        // split Omega into Axis_ and omega_
        omegai = mag(Omega);

        if (omegai < SMALL)
        {
            Axisi = vector::one;
        }
        else
        {
            Axisi =  Omega/(omegai+SMALL);
            forAll (Axisi,axElI)
            {
                if (mag(Axisi[axElI]) < 1.0e-08) Axisi[axElI] = 0.0;
            }
        }
        Axisi /= mag(Axisi);
    };

    auto updateRotationFixedAxis = [&]()
    {
        // update body angular velocity
        vector Omega(Axisi*omegai + deltaT * ( inv(getI()) & T_ ));

        // split Omega into Axis_ and omega_
        omegai = mag(Omega);

        vector newAxis = Omega/(omegai+SMALL);
        if ((newAxis & Axisi) < 0) Axisi *= (-1.0);;
    };

    if (bodyOperation_ == 0 or bodyOperation_ == 3)
    {
        if (mag(Veli)/(mag(Vel_)+SMALL) > maxDistInDEMloop_ || (omegai)/(omega_+SMALL) > maxDistInDEMloop_) return false;
        return true;
    }
    else if (bodyOperation_ == 1)
    {
        updateRotation();

        if (mag(Veli)/(mag(Vel_)+SMALL) > maxDistInDEMloop_ || (omegai)/(omega_+SMALL) > maxDistInDEMloop_) return false;
        return true;
    }
    else if (bodyOperation_ == 2)
    {
        updateTranslation();

        if (mag(Veli)/(mag(Vel_)+SMALL) > maxDistInDEMloop_ || (omegai)/(omega_+SMALL) > maxDistInDEMloop_) return false;
        return true;
    }
    else if (bodyOperation_ == 4)
    {
        updateRotationFixedAxis();

        if (mag(Veli)/(mag(Vel_)+SMALL) > maxDistInDEMloop_ || (omegai)/(omega_+SMALL) > maxDistInDEMloop_) return false;
        return true;
    }

    updateTranslation();
    if (updateTorque_) updateRotation();

    if (mag(Veli)/(mag(Vel_)+SMALL) > maxDistInDEMloop_ || (omegai)/(omega_+SMALL) > maxDistInDEMloop_) return false;
    return true;
}
//---------------------------------------------------------------------------//
void immersedBody::assignFullHistory()
{
    // assigned variables for potential contact correction
    historyAxis_ = Axis_;
    historyOmega_ = omega_;
    historyVel_ = Vel_;
    historya_ = a_;
    historyAlpha_ = alpha_;
    historyTotalAngle_ = totalAngle_;
    geomModel_->recordHistory();
}
//---------------------------------------------------------------------------//
void immersedBody::initSyncWithFlow(const volVectorField& U)
{
    // auxiliary computation (unnecessarily expensive)
    volVectorField curlU(fvc::curl(U));
    // Note (MI): if this initialization proves OK, than this needs to
    //            be computed only ONCE for all the bodies and re-used

    // computation itself
    vector meanV(vector::zero);
    scalar totVol(0);
    vector meanC(vector::zero);
    label  cellI;
    forAll (intCells_[Pstream::myProcNo()],iCellI)
    {
        cellI   = intCells_[Pstream::myProcNo()][iCellI];
        meanV  += U[cellI]*mesh_.V()[cellI];
        meanC  += curlU[cellI]*mesh_.V()[cellI];
        totVol += mesh_.V()[cellI];
    }
    reduce(meanV, sumOp<vector>());
    reduce(meanC, sumOp<vector>());
    reduce(totVol, sumOp<scalar>());
    Vel_ = meanV/(totVol+SMALL);
    meanC/=(totVol+SMALL);
    vector Omega(0.5*meanC);
    if(updateTorque_)
    {
        omega_ = mag(Omega);
        if (omega_ < SMALL)
        {
            Axis_ = vector::one;
            if (mesh_.nGeometricD() < 3)
            {
                const vector validDirs = (mesh_.geometricD() + Vector<label>::one)/2;
                Axis_ -= validDirs;
            }
        }
        else
        {
            Axis_ =  Omega/(omega_+SMALL);
            if (mesh_.nGeometricD() < 3)
            {
                const vector validDirs = (mesh_.geometricD() + Vector<label>::one)/2;
                Axis_ = cmptMultiply(vector::one-validDirs,Axis_);
            }
        }
        Axis_ /= mag(Axis_);
    }
    // update old storage
    VelOld_     = Vel_;
    omegaOld_   = omega_;
    AxisOld_    = Axis_;
    // print data:
    Info << "-- body " << bodyId_ << "initial movement variables:" << endl;
    printStats();
}
//---------------------------------------------------------------------------//
void immersedBody::pimpleUpdate
(
    volScalarField& body,
    volVectorField& f
)
{
    updateCoupling(body,f);
    updateMovement(VelOld_,AxisOld_,omegaOld_,velRelaxFac_);
}
//---------------------------------------------------------------------------//
void immersedBody::checkIfInDomain(volScalarField& body)
{
    if(getM0() < SMALL)
    {
        switchActiveOff(body);
    }

    Info << "M0: " << getM0() << endl;
    Info << "-- body " << bodyId_ << " current M/M0: " << getM()/getM0() << endl;
    // if only 1% of the initial particle mass remains in the domain, switch it off
    if (getM()/(getM0()+SMALL) < 1e-2)
    {
        switchActiveOff(body);
    }
}
//---------------------------------------------------------------------------//
void immersedBody::setRestartSim(vector vel, scalar angVel, vector axisRot, bool setStatic, label timesInContact)
{
    Vel_ = vel;
    omega_ = angVel;
    Axis_ = axisRot;
    contactInfo_->setTimeStepsInContWStatic(timesInContact);
    Info << "-- body " << bodyId_ << " timeStepsInContWStatic_: " << getTimeStepsInContWStatic() << endl;
    if(setStatic)
    {
        bodyOperation_ = 0;
        omega_ = 0;
        Vel_ *= 0;
        Info << "-- body " << bodyId_ << " set as Static" << endl;
    }
}
//---------------------------------------------------------------------------//
void immersedBody::chceckBodyOp()
{

    if(bodyOperation_ != 5 || timesToSetStatic_ == -1)
        return;

    if(!contactInfo_->checkInContactWithStatic() && getTimeStepsInContWStatic() > 0)
    {
    contactInfo_->setTimeStepsInContWStatic(0);
        return;
    }

    if(contactInfo_->checkInContactWithStatic())
    {
        contactInfo_->setTimeStepsInContWStatic(getTimeStepsInContWStatic() + 1);
        Info << "-- body " << bodyId_ << " timeStepsInContWStatic_: " << getTimeStepsInContWStatic() << endl;

        if(getTimeStepsInContWStatic() >= timesToSetStatic_)
        {
            bodyOperation_ = 0;
            omega_ = 0;
            Vel_ *= 0;
            Info << "-- body " << bodyId_ << " set as Static" << endl;
        }
    }

    contactInfo_->inContactWithStatic(false);
}
