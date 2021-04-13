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

    openHFDIB is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

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

//~ #include "addModel.H"
//~ #include "addModelOnce.H"

#define ORDER 2

using namespace Foam;

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
stlNames_(HFDIBDEMDict_.lookup("stlNames")),
kWN_(readScalar(HFDIBDEMDict_.subDict("wallProps").lookup("kN"))),
gammaWN_(readScalar(HFDIBDEMDict_.subDict("wallProps").lookup("gammaN"))),
kWt_(readScalar(HFDIBDEMDict_.subDict("wallProps").lookup("kt"))),
gammaWt_(readScalar(HFDIBDEMDict_.subDict("wallProps").lookup("gammat"))),
muW_(readScalar(HFDIBDEMDict_.subDict("wallProps").lookup("mu"))),
adhWN_(readScalar(HFDIBDEMDict_.subDict("wallProps").lookup("adhN"))),
minDEMloops_(readScalar(HFDIBDEMDict_.lookup("minDEMloops"))),
minDEMtimeStep_(readScalar(HFDIBDEMDict_.lookup("minDEMtimeStep")))
{
}
openHFDIBDEM::openHFDIBDEM(const Foam::dynamicFvMesh& mesh, const Foam::wordList stlNames )
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
stlNames_(stlNames),
kWN_(readScalar(HFDIBDEMDict_.subDict("wallProps").lookup("kN"))),
gammaWN_(readScalar(HFDIBDEMDict_.subDict("wallProps").lookup("gammaN"))),
kWt_(readScalar(HFDIBDEMDict_.subDict("wallProps").lookup("kt"))),
gammaWt_(readScalar(HFDIBDEMDict_.subDict("wallProps").lookup("gammat"))),
muW_(readScalar(HFDIBDEMDict_.subDict("wallProps").lookup("mu"))),
adhWN_(readScalar(HFDIBDEMDict_.subDict("wallProps").lookup("adhN"))),
minDEMloops_(readScalar(HFDIBDEMDict_.lookup("minDEMloops"))),
minDEMtimeStep_(readScalar(HFDIBDEMDict_.lookup("minDEMtimeStep")))
{
}
//---------------------------------------------------------------------------//
openHFDIBDEM::~openHFDIBDEM()
{
    //~ Info << "Calling openHFDIBDEM destructor" << endl;
    //~ forAll (immersedBodies_,bodyI)
    //~ {
        //~ Info << "Calling destructor for immersedBody " << bodyI << endl;
        //~ immersedBodies_[bodyI].~immersedBody();
        //~ Info << "Destroyed " << bodyI << endl;
    //~ }
    //~ Info << "Clearing out now empty pointer list" << endl;
    //~ immersedBodies_.~PtrList();
    //~ Info << "Cleared out" << endl;
    //~ forAll (addModels_,addModelI)
    //~ {
        //~ Info << "Calling destructor for addModel " << addModelI << endl;
        //~ addModels_[addModelI].~addModel();
    //~ }
    //~ addModels_.~PtrList();
    //~ Info << "Clearing out immersed bodies" << endl;
    //~ immersedBodies_.clear();
    //~ Info << "Clearing out auxiliary variables" << endl;
    //~ ibContactList_.clear();
    //~ Info << "About to exit openHFDIBDEM destructor" << endl;
}
//---------------------------------------------------------------------------//
//~ void openHFDIBDEM::initialize(volScalarField& body)
void openHFDIBDEM::initialize
(
    volScalarField& body,
    volVectorField& U
)
{
    // get data from HFDIBDEMDict
    wordList stlNames( HFDIBDEMDict_.lookup("stlNames") );
    HFDIBinterpDict_ = HFDIBDEMDict_.subDict("interpolationSchemes");
    
    // initialize addModels
    addModels_.setSize(stlNames_.size());
    #include "initializeAddModels.H"
    
    // register IBs as needed
    immersedBodies_.setSize(0);                                         //on the fly creation
    forAll (addModels_,modelI)
    {
        word stlName(stlNames_[modelI]);
        Info << "Creating immersed body based on: " << stlName << endl;
        
        label maxAdditions(100);
        label cAddition(0);
        // Note (MI): I hardcoded maxAdditions as our code is not
        //            efficient for high number of bodies
                
        while (addModels_[modelI].shouldAddBody(body) and cAddition < maxAdditions)
        {
            Info << "addModel invoked action, trying to add new body" << endl;
            triSurface bodySTL = addModels_[modelI].addBody(body);
            
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
                        bodySTL,
                        stlName,
                        mesh_,
                        HFDIBDEMDict_,
                        transportProperties_,
                        addIBPos
                    )
                );
                immersedBodies_[addIBPos].createImmersedBody(body);
                immersedBodies_[addIBPos].computeBodyCharPars();
                if (immersedBodies_[addIBPos].getStartSynced())
                {
                    immersedBodies_[addIBPos].initSyncWithFlow(U);
                }
                Info << "Body based on: " << stlName << " successfully added" << endl;
                cAddition = 0;
            }
            else
            {
                Info << "Body based on: " << stlName << " should have been added but was not "
                     << "(probably overlap with an already existing body)"
                     << endl;
            }
        }
    }
    
    // Initialize list for neighbour list method based on number of IBs
    boundValueNeighbourList_.setSize(3);
    boundLabelNeighbourList_.setSize(3);
    contactInCoordNeighbourList_.setSize(3);
    numberOfDEMloops_.setSize(Pstream::nProcs());
    for (label i = 0; i < 3; i = i + 1)
    {
        boundValueNeighbourList_[i].setSize(2 * immersedBodies_.size());
        boundLabelNeighbourList_[i].setSize(2 * immersedBodies_.size());
    }
    
    // Note: this needs to be done for all the actually registered bodies
    forAll (immersedBodies_,bodyI)
    {
        //Get references to bounding points and label list
        List<List<label>> IBboundList = immersedBodies_[bodyI].getBoundIndList();
        vector IBboundMinPoint = immersedBodies_[bodyI].getMinBoundPoint();
        vector IBboundMaxPoint = immersedBodies_[bodyI].getMaxBoundPoint();
        //Prepere label for bounding point identification. This is only list where bodyID is increased by
        //Minimal bounding point has negative value
        //Maximal bounding point has positive value
        label bodyIdInc(bodyI+1);

        //iterate over dimensions and assined proper bounding position with its label
        for (label coord = 0; coord < 3; coord = coord + 1)
        {
            boundValueNeighbourList_[coord][IBboundList[coord][0]] = IBboundMinPoint[coord];
            boundLabelNeighbourList_[coord][IBboundList[coord][0]] = -1 * bodyIdInc;
            boundValueNeighbourList_[coord][IBboundList[coord][1]] = IBboundMaxPoint[coord];
            boundLabelNeighbourList_[coord][IBboundList[coord][1]] = bodyIdInc;
        }
    }
    
    // After filling of the lists sort all bounding points.
    // This is only place where the sorting is not approximately O(N) efficient
    sortBoundingListPrtContact();
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::preUpdateBodies
(
    volScalarField& body,
    volVectorField& f
)
{    
    // - particles undergoing contact
    ibContactList_.clear();
    prtContactIBList_.clear();
    prtContactCenter_.clear();
    prtContactVolume_.clear();
    prtContactNormal_.clear();
    prtContactArea_.clear();
    
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            // print out linear and angular momentum (initial)
            //~ immersedBodies_[bodyId].printMomentum();
            //~ immersedBodies_[bodyId].printStats();
            
            // create body or compute body-fluid coupling and estimate
            // potential contacts with walls
            immersedBodies_[bodyId].preContactUpdateBodyField(body,f);
            
            immersedBodies_[bodyId].printStats();
            
            //~ Info << "-- body " << bodyId << "  linear velocity:  " << immersedBodies_[bodyId].getVel() << endl;
            //~ Info << "-- body " << bodyId << " angluar velocity: " << immersedBodies_[bodyId].getOmega() << endl;
            //~ Info << "-- body " << bodyId << " axis of rotation: " << immersedBodies_[bodyId].getAxis() << endl;
            
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
scalar openHFDIBDEM::postUpdateBodies
(
)
{    
    // compute max CoNum over the bodies
    scalar CoNum(0.0);
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].postContactUpdateImmersedBody();
            CoNum = max(CoNum,immersedBodies_[bodyId].getCoNum());
        }
    }
    return CoNum;
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::moveBodies
(
    volScalarField& body
)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            if (findIndex(ibContactList_, bodyId) == -1)
            {
                immersedBodies_[bodyId].postContactUpdateBodyField(body);
                // print out linear and angular momentum (final)
                //~ immersedBodies_[bodyId].printMomentum();
                immersedBodies_[bodyId].printStats();
            }
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::interpolateIB( volVectorField & V
                              ,volVectorField & Vs
                              ,volScalarField & body)
{

    //Create interpolator
    autoPtr<interpolation<vector>> interpV =
                   interpolation<vector>::New(HFDIBinterpDict_, V);
    //Reset imposed field
    Vs *= scalar(0);
    
    //Loop over all the immersed bodies
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            //Update imposed field according to body
            immersedBodies_[bodyId].updateVectorField(Vs, V.name(),body);
            
            //Get interpolation info a request list
            const List<DynamicList<immersedBody::intVecRequest>>& intVecReqList = immersedBodies_[bodyId].getinterpolationVecReqs();;
            List<DynamicList<immersedBody::interpolationInfo>>& intInfoList  = immersedBodies_[bodyId].getInterpolationInfo();
            
            List<DynamicPointList> intPointsToSend;
            List<DynamicLabelList> intCellsToSend;
            intPointsToSend.setSize(Pstream::nProcs());
            intCellsToSend.setSize(Pstream::nProcs());
            //Deal with all requests and ask processors for interpolation
            //pstreamBuffers works only with list so requests are rearenged
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
            //Receive interpolation requests from other processors
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
            //Compute interpolation for other processors
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
            //Send computed data back 
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
            //Receive interpolation data from other processors
            for (label proci = 0; proci < Pstream::nProcs(); proci++)
            {
                if (proci != Pstream::myProcNo())
                {
                    UIPstream recv(proci, pBufs);
                    DynamicVectorList recIntVec (recv);
                    intVecRcv[proci].append(recIntVec);
                }
            }
            //Assign received data
            for (label otherProci = 0; otherProci < intVecRcv.size(); otherProci++)
            {
                for (label intVecI = 0; intVecI < intVecRcv[otherProci].size(); intVecI++)
                {
                    intInfoList[Pstream::myProcNo()][intVecReqList[otherProci][intVecI].requestLabel_].intVec_[intVecReqList[otherProci][intVecI].vecLabel_] = intVecRcv[otherProci][intVecI];
                }
            }
            //Compute interpolation for cells found on this processor
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

            //loop over all interpolation info
            forAll (intInfoList[Pstream::myProcNo()],infoI)
            {
                label cellI = intInfoList[Pstream::myProcNo()][infoI].surfCell_;
                //Based on order calculat Vs and assign
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
                        
                        //distance between interpolation points
                        scalar deltaR = mag(intInfoList[Pstream::myProcNo()][infoI].intPoints_[1] -intInfoList[Pstream::myProcNo()][infoI].intPoints_[0]);
                        
                        //cell center to surface distance
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
                        // Note (MI): added +SMALL for robustness, effect on
                        //            solution not tested
                        
                        Vs[cellI] = linCoeff*ds + Vs[cellI];
                        break;
                    }
                    case 2:
                    {
                        vector VP1 =  intInfoList[Pstream::myProcNo()][infoI].intVec_[0] - Vs[cellI];
                        
                        vector VP2 =  intInfoList[Pstream::myProcNo()][infoI].intVec_[1] - Vs[cellI];
                        
                        
                        //distance between interpolation points
                        scalar deltaR1 = mag(intInfoList[Pstream::myProcNo()][infoI].intPoints_[2] -intInfoList[Pstream::myProcNo()][infoI].intPoints_[1]);
                        scalar deltaR2 = mag(intInfoList[Pstream::myProcNo()][infoI].intPoints_[1] -intInfoList[Pstream::myProcNo()][infoI].intPoints_[0]);
                        
                        //cell center to surface distance
                        scalar ds(0.0);
                        if (immersedBodies_[bodyId].getSDBasedLambda())
                        {
                            scalar minMaxBody(max(min(body[cellI],1.0-SMALL),SMALL));
                            ds = Foam::atanh(2.0*minMaxBody - 1.0)*Foam::pow(mesh_.V()[cellI],0.333);
                            ds/= immersedBodies_[bodyId].getIntSpan();
                        }
                        else
                        {
                            ds = deltaR1*(0.5-body[cellI]) ;
                        }
                        
                        vector quadCoeff = (VP2 - VP1)*deltaR1 - VP1*deltaR2;
                        quadCoeff       /= (deltaR1*deltaR2*(deltaR1 + deltaR2)+SMALL);
                        
                        vector linCoeff  = (VP1-VP2)*Foam::pow(deltaR1,2.0);
                        linCoeff        += 2.0*VP1*deltaR1*deltaR2;
                        linCoeff        += VP1*Foam::pow(deltaR2,2.0);
                        linCoeff        /= (deltaR1*deltaR2*(deltaR1 + deltaR2)+SMALL); 
        
                        //~ vector quadCoeff = 1.0/(deltaR*deltaR+SMALL) * ( VP2/2.0 - VP1 );
                        //~ vector linCoeff  = 1.0/(2.0*deltaR+SMALL) * ( 4.0*VP1 - VP2 );
                        // Note (MI): added +SMALL for the robustness, effects
                        //            on solution not tested
                        
                        Vs[cellI] = quadCoeff*ds*ds + linCoeff*ds + Vs[cellI];
                        break;
                    }
                }
            }
        }
    }
}
//---------------------------------------------------------------------------//
openHFDIBDEM::prtPrtContactInfo openHFDIBDEM::getPrtContactInfo(Tuple2<label, label> pairToTest)
{
    openHFDIBDEM::prtPrtContactInfo returnContactInfoValue;
    
    // reference to current body
    immersedBody& cBody(immersedBodies_[pairToTest.first()]);
    immersedBody& tBody(immersedBodies_[pairToTest.second()]);
    
    List<DynamicLabelList>    commonCells;
    commonCells.setSize(Pstream::nProcs());

    List<DynamicLabelList> cSurfCells(cBody.getSurfaceCellList());
    List<DynamicLabelList> cIntCells(cBody.getInternalCellList());
    List<DynamicLabelList> tSurfCells(tBody.getSurfaceCellList());
    
    //Bounding box of body t is used only for optimization of for loop
    boundBox tBox(tBody.getMinBoundPoint(),tBody.getMaxBoundPoint());
    // Inflate tBox to include cell center of surfCells
    tBox.inflate(2*sqrt(mesh_.magSf()[0])/tBox.mag());
    
    //Iterate over surfCells to find commen cells 
    forAll (cSurfCells[Pstream::myProcNo()],cSCellI)
    {
        if (tBox.contains(mesh_.C()[cSurfCells[Pstream::myProcNo()][cSCellI]]))
        {
            forAll (tSurfCells[Pstream::myProcNo()],tSCellI)
            {
                if (mag(cSurfCells[Pstream::myProcNo()][cSCellI]-tSurfCells[Pstream::myProcNo()][tSCellI]) < SMALL)
                {
                    commonCells[Pstream::myProcNo()].append(cSurfCells[Pstream::myProcNo()][cSCellI]);                    
                }
            }
        }
    }
    
    label numOfComCells(0);
    vector contactCenter(vector::zero);
    //If there are any common cells check if the surfaces are intersected
    if (commonCells[Pstream::myProcNo()].size() > SMALL)
    {
        //Evaluate center of the contact area
        forAll (commonCells[Pstream::myProcNo()], cCell)
        {
            contactCenter += mesh_.C()[commonCells[Pstream::myProcNo()][cCell]];
        }
        numOfComCells = commonCells[Pstream::myProcNo()].size();
    }
    
    sumReduce(contactCenter, numOfComCells);
    if (numOfComCells <= 0) 
    {
        returnContactInfoValue.prtsInContact_ = pairToTest;
        returnContactInfoValue.contactCenter_ = vector::zero;
        returnContactInfoValue.contactVolume_ = 0;
        returnContactInfoValue.contactNormal_ = vector::zero;
        returnContactInfoValue.contactArea_ = 0;
        returnContactInfoValue.inContact_ = false;
        return returnContactInfoValue;
    }
    
    contactCenter /= numOfComCells;
    
    vector cCoM(cBody.getCoM());        //center of mass
    vector tCoM(tBody.getCoM());        //center of mass
    scalar tDC(tBody.getDC());          //characteristic diameter
    
    //Create triSurfaceSearch
    const triSurface& tibTemp( tBody.getTriSurfMesh());
    triSurfaceSearch tibTriSurfSearch( tibTemp );
    
    //Get nearest point on surface from contact center
    pointIndexHit tibPointIndexHit = tibTriSurfSearch.nearest(contactCenter, vector::one * tDC);
    List<pointIndexHit> tibPointIndexHitList(1,tibPointIndexHit);
    vectorField normalVectorField;
    
    //Get contact normal direction
    const triSurfaceMesh& tibTempMesh( tBody.getTriSurfMesh());
    tibTempMesh.getNormal(tibPointIndexHitList,normalVectorField);
    // MS note: This normal direction is not correct for sharp edges! 
    // therefore it is used only if the normal direction cannot be set
    // based on surfcels (there is only one such cell
    
    vector normalVector(vector::zero);
    scalar contactArea(0);
    bool case3D(true);
    //Check if the case is 3D
    vector geomDirections(mesh_.geometricD());
    forAll (geomDirections, direction)
    {
        if (geomDirections[direction] == -1)
        {
            case3D = false;
            break;
        }
    }
    //Use edge cells to find contact area (better precision then surfcells)
    DynamicLabelList edgeCells;
        
    scalar intersectedVolume(0);
        
    if (commonCells[Pstream::myProcNo()].size() > SMALL)
    {

        //Create triSurfaceSearch
        const triSurface& cibTemp( cBody.getTriSurfMesh());
        triSurfaceSearch cibTriSurfSearch( cibTemp );
        const pointField& pp = mesh_.points();
        
        boolList tcenterInsideList = tibTriSurfSearch.calcInside( mesh_.C());
        boolList ccenterInsideList = cibTriSurfSearch.calcInside( mesh_.C());

        forAll (cSurfCells[Pstream::myProcNo()],cSCellI)
        {
            bool partiallyInT(false);
            
            if (tBox.contains(mesh_.C()[cSurfCells[Pstream::myProcNo()][cSCellI]]))
            {
                const labelList& vertexLabels = mesh_.cellPoints()[cSurfCells[Pstream::myProcNo()][cSCellI]];
                const pointField vertexPoints(pp,vertexLabels);
                boolList tvertexesInside = tibTriSurfSearch.calcInside( vertexPoints );
                boolList cvertexesInside = cibTriSurfSearch.calcInside( vertexPoints );
                bool tcenterInside(tcenterInsideList[cSurfCells[Pstream::myProcNo()][cSCellI]] );
                bool ccenterInside(ccenterInsideList[cSurfCells[Pstream::myProcNo()][cSCellI]] );
                //~ scalar rVInSize(1.0/(tvertexesInside.size()+1));// MI: WHIS IS THIS HERE?
                scalar rVInSize(0.5/tvertexesInside.size());
                // Note: weiggetContactNormalht of a single vertex in the cell
                
                scalar partialVolume(0);
                forAll (tvertexesInside, verIn)
                {
                    if (tvertexesInside[verIn]==true && cvertexesInside[verIn]==true)
                    {
                        partialVolume += rVInSize; //fraction of cell covered
                        partiallyInT = true;
                    }
                }

                if (tcenterInside==true && ccenterInside==true)
                {
                    partialVolume += 0.5; //fraction of cell covered
                    partiallyInT = true;
                }
                //Cells is edge cell when the cell is surfcell in both proccessors
                if (partialVolume + SMALL < 1 && partiallyInT)
                {
                    edgeCells.append(cSurfCells[Pstream::myProcNo()][cSCellI]);
                }
                
                
                intersectedVolume += mesh_.V()[cSurfCells[Pstream::myProcNo()][cSCellI]] * partialVolume;
            }
        }
        //Calculate remaining intersected volume
        forAll (cIntCells[Pstream::myProcNo()],cSCellI)
        {
            if (tBox.contains(mesh_.C()[cIntCells[Pstream::myProcNo()][cSCellI]]))
            {
                const labelList& vertexLabels = mesh_.cellPoints()[cIntCells[Pstream::myProcNo()][cSCellI]];
                const pointField vertexPoints(pp,vertexLabels);
                boolList tvertexesInside = tibTriSurfSearch.calcInside( vertexPoints );
                boolList cvertexesInside = cibTriSurfSearch.calcInside( vertexPoints );
                bool tcenterInside(tcenterInsideList[cIntCells[Pstream::myProcNo()][cSCellI]] );
                bool ccenterInside(ccenterInsideList[cIntCells[Pstream::myProcNo()][cSCellI]] );
                //~ scalar rVInSize(1.0/(tvertexesInside.size()+1));
                scalar rVInSize(0.5/tvertexesInside.size());
                // Note: weight of a single vertex in the cell
                
                scalar partialVolume(0);
                forAll (tvertexesInside, verIn)
                {
                    if (tvertexesInside[verIn]==true && cvertexesInside[verIn]==true)
                    {
                        partialVolume += rVInSize; //fraction of cell covered
                    }
                }
                
                if (tcenterInside==true && ccenterInside==true)
                {
                    partialVolume += 0.5; //fraction of cell covered
                }
                
                intersectedVolume += mesh_.V()[cIntCells[Pstream::myProcNo()][cSCellI]] * partialVolume;
            }
        }
    }
    
    reduce(intersectedVolume, sumOp<scalar>());
    
    if (intersectedVolume > 0)
    {
        if (case3D)
        {
            Tuple2<scalar,vector> returnTuple = get3DcontactInfo(edgeCells, normalVectorField[0], contactCenter, cBody.getOwner());
            contactArea = returnTuple.first();
            normalVector = returnTuple.second();
        }
        else
        {
            normalVector = getContactNormal(commonCells[Pstream::myProcNo()], normalVectorField[0], contactCenter);
                
            // Evaluate contact area
            contactArea = getContactArea(commonCells[Pstream::myProcNo()]);
        }
    }

    if (intersectedVolume > SMALL && contactArea > SMALL)
    {
        Info << "-- Particle-particle contact normal " << normalVector << endl;
        Info << "-- Particle-particle contact volume " << intersectedVolume << endl;
        Info << "-- Particle-particle contact area " << contactArea << endl;
        
        returnContactInfoValue.prtsInContact_ = pairToTest;
        returnContactInfoValue.contactCenter_ = contactCenter;
        returnContactInfoValue.contactVolume_ = intersectedVolume;
        returnContactInfoValue.contactNormal_ = normalVector;
        returnContactInfoValue.contactArea_ = contactArea;
        returnContactInfoValue.inContact_ = true;
        return returnContactInfoValue;
    }
    else
    {
        returnContactInfoValue.prtsInContact_ = pairToTest;
        returnContactInfoValue.contactCenter_ = vector::zero;
        returnContactInfoValue.contactVolume_ = 0;
        returnContactInfoValue.contactNormal_ = vector::zero;
        returnContactInfoValue.contactArea_ = 0;
        returnContactInfoValue.inContact_ = false;
        return returnContactInfoValue;
    }
}
//---------------------------------------------------------------------------//
Tuple2<scalar,vector> openHFDIBDEM::get3DcontactInfo(DynamicLabelList commonCells, vector normalVector, vector contactCenter, label owner)
{      
    //Collect edge positions from all processors
    List<DynamicPointList> commCellsPositionsProc;
    commCellsPositionsProc.setSize(Pstream::nProcs());
    DynamicPointList commCellsPositions;
    forAll (commonCells, cCell)
    {
        commCellsPositionsProc[Pstream::myProcNo()].append(mesh_.C()[commonCells[cCell]]);
    }
    
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);    
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UOPstream send(proci, pBufs);
            send << commCellsPositionsProc[Pstream::myProcNo()];
        }
    }
    
    pBufs.finishedSends();    
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recv(proci, pBufs);
            DynamicPointList commCellsPositionsi (recv);
            commCellsPositionsProc[proci].append(commCellsPositionsi);
        }
    }
    //First fill edge positions from owner (better stability for svd)
    for (label i = 0; i < commCellsPositionsProc[owner].size(); i++)
    {
        commCellsPositions.append(commCellsPositionsProc[owner][i]);
    }
    
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != owner)
        {
            for (label i = 0; i < commCellsPositionsProc[proci].size(); i++)
            {
                commCellsPositions.append(commCellsPositionsProc[proci][i]);
            }
        }
    }
    
    scalar area(0.0);
    vector normalVec(vector::zero);
    
    if (commCellsPositions.size() > 3)
    {    
        //Find best fitting plane for given positions
        RectangularMatrix<scalar> matrix(3,commCellsPositions.size());
        for (label i = 0; i < commCellsPositions.size(); i++)
        {
            vector subPoint(commCellsPositions[i] - contactCenter);
            
            matrix[0][i] = subPoint[0];
            matrix[1][i] = subPoint[1];
            matrix[2][i] = subPoint[2];
        }
        
        bool svdOk(true);
        vector normalVec(vector::zero);
        try{
            SVD svd(matrix);
            label vecInd(svd.U().n() - svd.nZeros() - 1);
            normalVec[0] = svd.U()[0][vecInd];
            normalVec[1] = svd.U()[1][vecInd];
            normalVec[2] = svd.U()[2][vecInd];
            if (mag(normalVec) < SMALL)
                normalVec = normalVector;
            
            if ((normalVec & normalVector) < 0)
                normalVec *= -1;
        }
        catch(...)
        {
            svdOk = false;
        }
        reduce(svdOk, andOp<bool>());        
        if (!svdOk)
            normalVec = normalVector;

        //Create best fitting plane
        plane bestFitPlane(contactCenter, normalVec);
        normalVec = bestFitPlane.normal();
        DynamicPointList commCellsPosInPlane;    
        forAll (commCellsPositions,cell)
            commCellsPosInPlane.append(bestFitPlane.nearestPoint(commCellsPositions[cell]));
        
        vector q1(1.0, 0.0, 0.0);
        vector q2(1.0, 0.0, 0.0);    
        if (abs(q1 & bestFitPlane.normal()) > abs(q2 & bestFitPlane.normal()))
            q1 = q2;
        
        vector u(bestFitPlane.normal() ^ q1);
        vector v(bestFitPlane.normal() ^ u); 
        
        DynamicList<plane> clockwisePlanes;
        List<scalar> helpList(6);
        helpList[0] = 0.0;
        helpList[1] = 1.0;
        helpList[2] = 0.0;
        helpList[3] = -1.0;
        helpList[4] = 0.0;
        helpList[5] = 1.0;
        
        //Loop over parts of plane to find and sort points
        DynamicVectorList commCellsInSections; 
        for (label i = 0; i < 4; i++)
        {
            scalar uStep(helpList[i + 1] -  helpList[i]);
            scalar vStep(helpList[i + 2] -  helpList[i + 1]);
            DynamicPointList pointsInSection;
            for (scalar j = 0.0; j < 3.0; j += 1.0)
            {
                plane uPlane(contactCenter, u*(helpList[i] + uStep*j/4.0) + v*(helpList[i + 1] + vStep*j/4.0));
                plane vPlane(contactCenter, u*(helpList[i] + uStep*(j+1)/4.0) + v*(helpList[i + 1] + vStep*(j+1)/4.0));
                
                forAll (commCellsPosInPlane, celli)
                {
                    if (uPlane.sideOfPlane(commCellsPosInPlane[celli]) == 0 && vPlane.sideOfPlane(commCellsPosInPlane[celli]) == 1)
                        pointsInSection.append(commCellsPosInPlane[celli]);
                }
                
                if (pointsInSection.size() > SMALL)
                {
                    vector average(vector::zero);
                    forAll (pointsInSection, pointI)
                        average += pointsInSection[pointI];                
                    average /= pointsInSection.size();
                    
                    commCellsInSections.append(average);
                }
            }
        }
        //Calculate contact area
        for (label i = 0; i + 1 < commCellsInSections.size(); i++)
        {
            vector AC(commCellsInSections[i] - contactCenter);
            vector BC(commCellsInSections[i + 1] - contactCenter);
            vector crossPr(AC ^ BC);
            area += mag(crossPr)/2;
        }
        
        if (commCellsInSections.size() > 2)
        {
            vector AC(commCellsInSections[commCellsInSections.size() - 1] - contactCenter);
            vector BC(commCellsInSections[0] - contactCenter);
            vector crossPr(AC ^ BC);
            area += mag(crossPr)/2;
        }
    }
    
    if (normalVec == vector::zero)
        normalVec = normalVector;
    
    reduce(area, sumOp<scalar>());
    area /= Pstream::nProcs();
    reduce(normalVec, sumOp<vector>());
    normalVec /= Pstream::nProcs();
    
    Tuple2<scalar,vector> returnValue(area,normalVec);

    return returnValue;
}
//---------------------------------------------------------------------------//
scalar openHFDIBDEM::getContactArea(DynamicLabelList commonCells)
{        
    DynamicPointList commCellsPositions;
    
    forAll (commonCells, cCell)
    {
        commCellsPositions.append(mesh_.C()[commonCells[cCell]]);
    }
    
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UOPstream send(proci, pBufs);
            send << commCellsPositions;
        }
    }
    
    pBufs.finishedSends();
    
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recv(proci, pBufs);
            DynamicPointList commCellsPositionsi (recv);
            commCellsPositions.append(commCellsPositionsi);
        }
    }
    
    // Evaluate contact area
    scalar contactArea(0);
        
    // MS note: this is dummy way how to calculate contactArea in 2D -> discuss!!
    scalar highestDistance(0);
    
    for (label i = 1; i < commCellsPositions.size(); i++)
    {
        if (mag(commCellsPositions[i-1] - commCellsPositions[i]) > highestDistance)
        {
            highestDistance = mag(commCellsPositions[i-1] - commCellsPositions[i]);
        }
    }
    
    reduce(highestDistance, maxOp<scalar>());
    if (commonCells.size() > 0)
        contactArea = sqrt(mag(mesh_.Sf()[commonCells[0]])) * highestDistance;
    reduce(contactArea, maxOp<scalar>());
    return contactArea;
}
//---------------------------------------------------------------------------//
vector openHFDIBDEM::getContactNormal(DynamicLabelList commonCells, vector normalVector, vector contactCenter)
{    
    // Ms note: mesh_.geometricD works only with empty type of boundary!
    vector geomDirections(mesh_.geometricD());
    forAll (geomDirections, direction)
    {
        if (geomDirections[direction] == 1)
            geomDirections[direction] = 0;
        if (geomDirections[direction] == -1)
            geomDirections[direction] = 1;
    }
    
    vector normalVector2(vector::zero);
    
    forAll (commonCells, cCell)
    {
        vector tempor((mesh_.C()[commonCells[cCell]] - contactCenter) ^ geomDirections);
        if ((tempor & normalVector) < 0)
            tempor *= -1;
        normalVector2 += tempor;
    }
    
    label numOfComCells(commonCells.size());
    sumReduce(normalVector2, numOfComCells);

    normalVector2 /= numOfComCells;
    
    if (mag(normalVector2) == 0) return normalVector;
    
    normalVector2 = normalVector2/mag(normalVector2);
        
    // Assure that the normal vector points out of the body
    if ((normalVector2 & normalVector) < 0)
        normalVector2 *= -1;
    
    return normalVector2;
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::solvePrtContact(openHFDIBDEM::prtPrtContactInfo contactInfo, scalar deltaT)
{
    // reference to current body
    immersedBody& cBody(immersedBodies_[contactInfo.prtsInContact_.first()]);
    
    // get data from the current body necessary for computation of
    // contact forces
    vector cCoM(cBody.getCoM());        //center of mass
    vector cVel(cBody.getVel());        //linear velocity
    scalar cKN(cBody.getKN());          //elastic normal stiffness
    scalar cGammaN(cBody.getGammaN());  //normal viscosity
    scalar cKt(cBody.getKt());          //elastic tangential stiffness
    scalar cGammat(cBody.getGammat());  //tangential viscosity
    scalar cmu(cBody.getmu());          //firction coef
    scalar cadhN(cBody.getadhN());      //adhesive force
    scalar cM(cBody.getM());            //mass

    immersedBody& tBody(immersedBodies_[contactInfo.prtsInContact_.second()]);
    
    vector tCoM(tBody.getCoM());        //center of mass
    vector tVel(tBody.getVel());        //linear velocity
    scalar tKN(tBody.getKN());          //elastic normal stiffness
    scalar tGammaN(tBody.getGammaN());  //normal viscosity
    scalar tKt(tBody.getKt());          //elastic tangential stiffness
    scalar tGammat(tBody.getGammat());  //tangential viscosity
    scalar tmu(tBody.getmu());          //firction coef
    scalar tadhN(tBody.getadhN());      //adhesive force
    scalar tM(tBody.getM());            //mass
    
    
    vector  FN(vector::zero);//placeholder for normal force
    vector  Ft(vector::zero);//placeholder for tangential force
    Info << "-- Detected particle-particle contact " << contactInfo.prtsInContact_.first() << " : " << contactInfo.prtsInContact_.second() << endl;
    
    // compute mean model parameters
    scalar aKN(0.5*(cKN+tKN));
    scalar aGammaN(0.5*(cGammaN+tGammaN));
    scalar aKt(0.5*(cKt+tKt));
    scalar aGammat(0.5*(cGammat+tGammat));
    scalar amu(0.5*(cmu+tmu));
    scalar aadhN(0.5*(cadhN+tadhN));
    
    vector cLVec(contactInfo.contactCenter_-cCoM);
    vector tLVec(contactInfo.contactCenter_-tCoM);
    
    // compute normal to movement and relative velocity
    vector nVec(contactInfo.contactNormal_);
    scalar Vn(-(cVel - tVel) & nVec);
    
    scalar Lc(4*mag(cLVec)*mag(tLVec)/(mag(cLVec)+mag(tLVec)));
    
    scalar reduceM(cM*tM/(cM+tM));

    // compute the normal force
    FN = (aKN*contactInfo.contactVolume_/(Lc+SMALL) + aGammaN*sqrt(aKN*reduceM/pow(Lc+SMALL,3))*(contactInfo.contactArea_ * Vn))*nVec;
    
    // compute adhesive force
    vector FA(aadhN*contactInfo.contactArea_*nVec);
    
    FN -= FA;
    
    vector cFtLast(vector::zero);
    vector tFtLast(vector::zero);
    
    // Find history of tangential force between these two particles
    DynamicList<Tuple2<label,Tuple2<label,vector>>> chistoryFt(cBody.getHistoryhistoryFt());
    DynamicList<Tuple2<label,Tuple2<label,vector>>> thistoryFt(tBody.getHistoryhistoryFt());
    
    bool cFtLastFinded(false);
    forAll (chistoryFt,cFti)
    {
        if (chistoryFt[cFti].first() == tBody.getBodyId())
        {
            cFtLastFinded = true;
            cFtLast = chistoryFt[cFti].second().second();
            break;
        }
    }
    
    bool tFtLastFinded(false);
    forAll (thistoryFt,tFti)
    {
        if (thistoryFt[tFti].first() == cBody.getBodyId())
        {
            tFtLastFinded = true;
            tFtLast = thistoryFt[tFti].second().second();
            break;
        }
    }
    
    // Note the magnitude of last Ft should be same but only oposite
    vector FtLast((cFtLast - tFtLast)/2);
    //Project last Ft into a new direction
    vector FtLastP(FtLast - (FtLast & nVec) * nVec);
    //Scale projected Ft to have same magnitude as FtLast
    vector FtLastr(mag(FtLast) * (FtLastP/(mag(FtLastP)+SMALL)));
    // Compute relative tangential velocity
    vector Vt((cVel - tVel) - ((cVel - tVel) & nVec) * nVec);
    //Compute tangential force
    Ft = (FtLastr - aKt*Vt*deltaT - aGammat*Vt);
    
    if (mag(Ft) > amu * mag(FN))
    {
        Ft *= amu * mag(FN) / mag(Ft);
    }
    
    // add the computed force to the affected bodies
    vector cTN(cLVec ^  FN);
    vector tTN(tLVec ^ -FN);
    cBody.updateFAndT( FN+Ft,cTN);
    tBody.updateFAndT(-FN-Ft,tTN);
    
    //Update history of tangential force
    if (cFtLastFinded)
    {
        forAll (chistoryFt,cFti)
        {
            if (chistoryFt[cFti].first() == tBody.getBodyId())
            {
                Tuple2<label,vector> help(1,Ft);
                chistoryFt[cFti].second() = help;
                break;
            }
        }
    }
    else
    {
        Tuple2<label,vector> help(1,Ft);
        Tuple2<label,Tuple2<label,vector>> help2(tBody.getBodyId(), help);
                
        chistoryFt.append(help2);
    }
    
    if (tFtLastFinded)
    {
        forAll (thistoryFt,tFti)
        {
            if (thistoryFt[tFti].first() == cBody.getBodyId())
            {
                Tuple2<label,vector> help(1,-Ft);
                thistoryFt[tFti].second() = help;
                break;
            }
        }
    }
    else
    {
        Tuple2<label,vector> help(1,-Ft);
        Tuple2<label,Tuple2<label,vector>> help2(cBody.getBodyId(), help);
                
        thistoryFt.append(help2);
    }
}
//---------------------------------------------------------------------------//
//~ void openHFDIBDEM::solveWallContact(const DynamicList<label>& wallContactLst)
//~ {
//~ }
// Note (MI): moved to immersedBody class by MS (as it should be there)

//---------------------------------------------------------------------------//
void openHFDIBDEM::writeBodySurfMeshes()
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].writeBodySurfMesh();
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::correctContact(volScalarField& body)
{
    //First detect prt contact to be sure that all Contact during CFD step will be solved during DEM inner loops
    detectPrtContact();
    // Get new contacts with wall and return their position
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
    }
    
    forAll (ibPrtContactList,ib)
    {
        label ibIndex(-1);
        ibIndex = findIndex(ibContactList_, ibPrtContactList[ib]);
        if (ibIndex > -1)
        {
            //remove at idex
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
    //Iterate over IBs not in time
    //First resolve prt-prt contacts
    while(prtContactIBList_.size() > 0)
    {
        pos = 0.0;
        step = maxDemStep;
        historyPos = 0.0;
        doubleMinStep = false;
        doubleMinStep2 = false;
        
        DynamicList<Tuple2<label, label>> pairsToResolve;
        //We have to resolve whole contact chain not only two prts in contact
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
        
        //Return IBs to begging position
        forAll (ibToResolve, ib)
        {
            immersedBodies_[ibToResolve[ib]].initializeVarHistory(true);
            immersedBodies_[ibToResolve[ib]].resetBody(body, false);
            immersedBodies_[ibToResolve[ib]].createImmersedBody(body, true);
        }
        
        //Compute motion over time. Particles are returned when bad motion occurs
        while(true)
        {
            Info << " Start DEM pos: " << pos << " DEM step: " << step << endl;
            forAll (ibToResolve, ib)
            {                
                // Set F_ and T_ to zero. Do not assign history values
                immersedBodies_[ibToResolve[ib]].initializeVarHistory(false);
                // Add fluid coupling force to F and T. This is still same for whole DEM inner loop
                immersedBodies_[ibToResolve[ib]].updateFAndT(immersedBodies_[ibToResolve[ib]].getHistoryCouplingF(), immersedBodies_[ibToResolve[ib]].getHistoryCouplingT());
            }
            //Find prt-prt contact info
            forAll (pairsToResolve, pair)
            {
                openHFDIBDEM::prtPrtContactInfo pairContactInfoValue;
                pairContactInfoValue = getPrtContactInfo(pairsToResolve[pair]);
                if (pairContactInfoValue.inContact_)
                {
                    solvePrtContact(pairContactInfoValue, deltaTime*step);
                }
            }
            
            DynamicList<bool> contactOk;            
            forAll (ibToResolve,ib)
            {
                // Detect wall contact and solve it
                // Update movement and move bodies
                immersedBodies_[ibToResolve[ib]].detectWallContact(body);
                if (immersedBodies_[ibToResolve[ib]].checkWallContact()) 
                {
                    immersedBodies_[ibToResolve[ib]].solveWallContact(kWN_, gammaWN_, kWt_, gammaWt_, muW_, adhWN_, deltaTime*step);
                }

                contactOk.append(immersedBodies_[ibToResolve[ib]].checkContactMovement(deltaTime*step));
            }
            //If there is any bad motion return all particles currently moving
            bool contactIsOk(true);
            forAll (contactOk,isOK)
            {
                if (!contactOk[isOK])
                {
                    contactIsOk = false;
                    break;
                }
            }
            
            //Move or return particles
            forAll (ibToResolve,ib)
            {
                if (contactIsOk || step == minDEMtimeStep_)
                {
                    immersedBodies_[ibToResolve[ib]].assignFullHistory();
                    immersedBodies_[ibToResolve[ib]].updateMovement(deltaTime*step);
                    immersedBodies_[ibToResolve[ib]].resetBody(body);
                    immersedBodies_[ibToResolve[ib]].moveImmersedBody(deltaTime*step);
                    immersedBodies_[ibToResolve[ib]].createImmersedBody(body, false, deltaTime*step);
                }
                else
                {
                    immersedBodies_[ibToResolve[ib]].initializeVarHistory(true);
                    immersedBodies_[ibToResolve[ib]].returnPosition();
                    immersedBodies_[ibToResolve[ib]].resetBody(body, false);
                    immersedBodies_[ibToResolve[ib]].createImmersedBody(body, true);
                }                
            }
            //Change dem step. If Step is already minimal wait to avoid circulation
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
            //If moved to end time remove from list or terminate loop
            if (pos >= 1)
            {
                forAll (pairsToResolve, pair)
                {
                    label indexOfPair(-1);
                    indexOfPair = findIndexOfPairInNeighbourList(prtContactIBList_, pairsToResolve[pair]);
                    //Remove at index
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
    //Resolve wall contacts 
    while(ibContactList_.size() > 0)
    {
        pos = 0.0;
        step = maxDemStep;
        historyPos = 0.0;
        doubleMinStep = false;
        doubleMinStep2 = false;
        
        immersedBodies_[ibContactList_[0]].initializeVarHistory(true);
        immersedBodies_[ibContactList_[0]].resetBody(body, false);
        immersedBodies_[ibContactList_[0]].createImmersedBody(body, true);
        
        while(true)
        {             
            Info << " Start DEM pos: " << pos << " DEM step: " << step << endl;
            // Set F_ and T_ to zero. Do not assign history values
            immersedBodies_[ibContactList_[0]].initializeVarHistory(false);
            // Add fluid coupling force to F and T. This is still same for whole DEM inner loop
            immersedBodies_[ibContactList_[0]].updateFAndT(immersedBodies_[ibContactList_[0]].getHistoryCouplingF(), immersedBodies_[ibContactList_[0]].getHistoryCouplingT());
            
            // Detect wall contact and solve it
            // Update movement and move bodies
            immersedBodies_[ibContactList_[0]].detectWallContact(body);
            if (immersedBodies_[ibContactList_[0]].checkWallContact()) 
            {
                immersedBodies_[ibContactList_[0]].solveWallContact(kWN_, gammaWN_, kWt_, gammaWt_, muW_, adhWN_, deltaTime*step);
            }

            
            bool contactIsOk(immersedBodies_[ibContactList_[0]].checkContactMovement(deltaTime*step));

            if (contactIsOk || step == minDEMtimeStep_)
            {
                immersedBodies_[ibContactList_[0]].assignFullHistory();
                immersedBodies_[ibContactList_[0]].updateMovement(deltaTime*step);
                immersedBodies_[ibContactList_[0]].resetBody(body);
                immersedBodies_[ibContactList_[0]].moveImmersedBody(deltaTime*step);
                immersedBodies_[ibContactList_[0]].createImmersedBody(body, false, deltaTime*step);
            }
            else
            {
                immersedBodies_[ibContactList_[0]].initializeVarHistory(true);
                immersedBodies_[ibContactList_[0]].returnPosition();
                immersedBodies_[ibContactList_[0]].resetBody(body, false);
                immersedBodies_[ibContactList_[0]].createImmersedBody(body, true);
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
}
//---------------------------------------------------------------------------//
// void openHFDIBDEM::correctContact(volScalarField& body)
// {
//     //First detect prt contact to be sure that all Contact during CFD step will be solved during DEM inner loops
//     detectPrtContact();
//     // Get new contacts with wall and return their position
//     getContactListAndReturnPositions(body);
//     scalar deltaTime(mesh_.time().deltaT().value());
//     scalar time(0.0);
//     //Create bool if there is any prt-prt contact
//     bool prtContactInTimeStep(prtContactIBList_.size() > SMALL);
//     // For each prt-prt contact pair assing IB to contact list and
//     // Return its position
//     forAll (prtContactIBList_,pair)
//     {
//         if (findIndex(ibContactList_, prtContactIBList_[pair].first()) == -1)
//         {
//             ibContactList_.append(immersedBodies_[prtContactIBList_[pair].first()].getBodyId());
//             immersedBodies_[prtContactIBList_[pair].first()].returnPosition();
//         }
//         if (findIndex(ibContactList_, prtContactIBList_[pair].second()) == -1)
//         {
//             ibContactList_.append(immersedBodies_[prtContactIBList_[pair].second()].getBodyId());
//             immersedBodies_[prtContactIBList_[pair].second()].returnPosition();
//         }
//     }
//     
//     scalar finesTimeRes(1);
//     // Assign history values createImmersed body on the history position
//     // Find out smallest resolution time overall to satisfy number of DEM loops
//     // Note (MS): All contacts are resolved in same number of DEM loops
//     //              If there is one contact which require a lot of DEM loops
//     //              But a lot of slow contact this is inefficient
//     //              Maybe we should think about adaptive time step in a way
//     //              that each contact will have its own number of DEM loops
//     forAll (ibContactList_,bodyId)
//     {
//         immersedBodies_[ibContactList_[bodyId]].initializeVarHistory(true);
//         immersedBodies_[ibContactList_[bodyId]].resetBody(body, false);
//         immersedBodies_[ibContactList_[bodyId]].createImmersedBody(body, true);
//         if (immersedBodies_[ibContactList_[bodyId]].getcontactTimeRes() < finesTimeRes)
//         {
//             finesTimeRes =immersedBodies_[ibContactList_[bodyId]].getcontactTimeRes();
//         }
//     }
//     // Estimate number of DEM loops
//     // If all contacts are slow we should use the minimal number of DEM loops
//     // Note (MS): Contacts may be resolved badly when there is only one DEM inner loop    
//     scalar numberOfDEMloopsi(ceil(deltaTime/finesTimeRes));
//     numberOfDEMloops_[Pstream::myProcNo()] = max(minDEMloops_,numberOfDEMloopsi);
//     
//     Pstream::gatherList(numberOfDEMloops_, 0);
//     Pstream::scatter(numberOfDEMloops_, 0);
//     
//     deltaTime /= max(numberOfDEMloops_);
//     
//     Info << "-- numberOfInsideLoops " << max(numberOfDEMloops_) << endl;
//     
//     for (int i = 0; i <= max(numberOfDEMloops_); i++)
//         //(time < mesh_.time().deltaT().value() + SMALL)
//     {
//         Info << "-- numberOfInsideLoops in loop: " << floor(time/deltaTime) << " : " << floor(mesh_.time().deltaT().value()/deltaTime) << endl;
//         forAll (ibContactList_,bodyId)
//         {
//             // Set F_ and T_ to zero. Do not assign history values
//             immersedBodies_[ibContactList_[bodyId]].initializeVarHistory(false);
//             // Add fluid coupling force to F and T. This is still same for whole DEM inner loop
//             immersedBodies_[ibContactList_[bodyId]].updateFAndT(immersedBodies_[bodyId].getHistoryCouplingF(), immersedBodies_[ibContactList_[bodyId]].getHistoryCouplingT());
//         }
//         
//         if (prtContactInTimeStep)
//         {
//             // Update neighbout lists for current positions; detect prt-prt contact and solve the contacts
//             updateNeighbourLists();
//             prtContactIBList_.clear();
//             prtContactCenter_.clear();
//             prtContactVolume_.clear();
//             prtContactNormal_.clear();
//             prtContactArea_.clear();
//             detectPrtContact();
//             solvePrtContact(deltaTime);
//         }
// 
//         forAll (ibContactList_,bodyId)
//         {
//             // Detect wall contact and solve it
//             // Update movement and move bodies
//             immersedBodies_[ibContactList_[bodyId]].detectWallContact(body);
//             if (immersedBodies_[ibContactList_[bodyId]].checkWallContact()) 
//             {
//                 immersedBodies_[ibContactList_[bodyId]].solveWallContact(kWN_, gammaWN_, kWt_, gammaWt_, muW_, adhWN_, deltaTime);
//             }
//             immersedBodies_[ibContactList_[bodyId]].updateMovement(deltaTime);
//             immersedBodies_[ibContactList_[bodyId]].resetBody(body);
//             immersedBodies_[ibContactList_[bodyId]].moveImmersedBody(deltaTime);
//             immersedBodies_[ibContactList_[bodyId]].createImmersedBody(body, false, deltaTime);
//             
//         }
//         time += deltaTime;
//     }
// }
//---------------------------------------------------------------------------//
void openHFDIBDEM::getContactListAndReturnPositions(volScalarField& body)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            if (findIndex(ibContactList_, bodyId) == -1)
            {   
                //Detect wall contact and assign the body to contact list and return its position
                immersedBodies_[bodyId].detectWallContact(body);
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
    //Use insertion sort to sort the neighbour list
    //Besides the first sorting, this is O(N) efficient
    for (label coord = 0; coord < 3; coord = coord + 1)
    {
        label i(1);
        while(i < boundValueNeighbourList_[coord].size())
        {
            label j(i);
            while(j > 0 && boundValueNeighbourList_[coord][j-1] > boundValueNeighbourList_[coord][j])
            {
                //If we should swap two bounding points resolve this action
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
    // Swap actual values of bounding points
    scalar scratchScalar(boundValueNeighbourList_[coord][j-1]);
    boundValueNeighbourList_[coord][j-1] = boundValueNeighbourList_[coord][j];
    boundValueNeighbourList_[coord][j] = scratchScalar;
    
    // Update identification list for both bodies
    // Note (MS): We have to resolve minimal and maximal point in different manner
    //              Therefore, if the label is negative (minimal point) we have to excess in list by its -value
    //              keep in mind that boundLabelNeighbourList_ does not keep labels by bodyID but bodyID+1
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
    
    //If we are actualy swapping minimal and maximal point there is change in intersection for these bodies
    //This means that multiplication of labels is negative (minimal is negative label and maximal is pozitive label)
    if (boundLabelNeighbourList_[coord][j-1] * boundLabelNeighbourList_[coord][j] < 0)
    {
        // If maximal point had lower index there is new intersection
        if (boundLabelNeighbourList_[coord][j-1] > 0)
        {
            //Create a Tuple2 for bodies that are newly intersected always keep first the body with lower bodyID
            Tuple2<label, label> newPair(min(boundLabelNeighbourList_[coord][j-1] - 1,-1*boundLabelNeighbourList_[coord][j] - 1), max(boundLabelNeighbourList_[coord][j-1] - 1,-1*boundLabelNeighbourList_[coord][j] - 1));
            // if this intersection is not already in possible contact list for this coord assign it
            // Note: This should never happen but it prevent multiple assignment
            if (findIndexOfPairInNeighbourList(contactInCoordNeighbourList_[coord], newPair) == -1)
            {
                contactInCoordNeighbourList_[coord].append(newPair);
                
                bool appendToPossibleContactList(true);
                //Check if this pair is already assigned for the remaining dimensions
                //If so, add this pair to possible contact list because their bounding boxes are intersected
                for (label coordi = 0; coordi < 3; coordi = coordi + 1)
                {
                    if (coordi != coord)
                    {
                        if (findIndexOfPairInNeighbourList(contactInCoordNeighbourList_[coordi],newPair) == -1)
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
            // There is no intersection for this coord any more for this situation so the pair should be removed from lists
            Tuple2<label, label> removePair(min(-1*boundLabelNeighbourList_[coord][j-1] - 1,boundLabelNeighbourList_[coord][j] - 1), max(-1*boundLabelNeighbourList_[coord][j-1] - 1,boundLabelNeighbourList_[coord][j] - 1));
            
            label indexOfpair(-1);
            // Remove the pair from coord contact list
            // Note: It is not so straightforward. Create help list and assign to it values before and after the removed pair
            indexOfpair = findIndexOfPairInNeighbourList(contactInCoordNeighbourList_[coord], removePair);
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
            
            // Same situation for contact list. If there is not instersection in one dimension the 
            // bounding boxes are not intersected so we should remove this pair from contact list
            indexOfpair = findIndexOfPairInNeighbourList(possibleContactNeighbourList_, removePair);
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
    
    // Actually swap the labels 
    label scratchLabel(boundLabelNeighbourList_[coord][j-1]);
    boundLabelNeighbourList_[coord][j-1] = boundLabelNeighbourList_[coord][j];
    boundLabelNeighbourList_[coord][j] = scratchLabel;
}
//---------------------------------------------------------------------------//
// Find index of pair in given list. Return -1 if the pair is not in the list
label openHFDIBDEM::findIndexOfPairInNeighbourList(DynamicList<Tuple2<label, label>>& listToSearch, Tuple2<label, label> pair)
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
// Find index of pairs that contains at least one ib from given pair
void openHFDIBDEM::findIndexesOfPairWithSomeIb(DynamicList<Tuple2<label, label>>& listToSearch, Tuple2<label, label> pair, DynamicList<Tuple2<label, label>>& listToAppend)
{
    DynamicList<Tuple2<label, label>> helpList;
    
    for (label i = 0; i < listToSearch.size(); i = i +1)        
    {
        if (listToSearch[i].first() == pair.first() || listToSearch[i].first() == pair.second() || listToSearch[i].second() == pair.first() || listToSearch[i].second() == pair.second())
        {
            if (findIndexOfPairInNeighbourList(listToAppend,listToSearch[i]) == -1)
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
//Based on identification list of Bodies update the values of minimal and maximal bounding points
//and sort the neighbour list
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
//~ void openHFDIBDEM::addRemoveBodies(volScalarField& body)
void openHFDIBDEM::addRemoveBodies
(
    volScalarField& body,
    volVectorField& U
)
{   
    forAll (addModels_,modelI)
    {
        word stlName(stlNames_[modelI]);
        
        label maxAdditions(100);
        label cAddition(0);
        // Note (MI): I hardcoded maxAdditions as our code is not
        //            efficient for high number of bodies
        
        while (addModels_[modelI].shouldAddBody(body) and cAddition < maxAdditions)
        {
            Info << "addModel invoked action, trying to add new body" << endl;
            triSurface bodySTL = addModels_[modelI].addBody(body);
            
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
                        bodySTL,
                        stlName,
                        mesh_,
                        HFDIBDEMDict_,
                        transportProperties_,
                        addIBPos
                    )
                );
                
                // get reference for further processing
                immersedBody& nBody(immersedBodies_[addIBPos]);
                nBody.createImmersedBody(body);
                nBody.computeBodyCharPars();
                if (nBody.getStartSynced())
                {
                    nBody.initSyncWithFlow(U);
                }
                
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
    //Check only pairs whose bounding boxes are intersected for the contact
    forAll (possibleContactNeighbourList_,possiblePair)
    {
        //Check only if the pair is not alredy in contactList
        if (findIndexOfPairInNeighbourList(prtContactIBList_,possibleContactNeighbourList_[possiblePair]) == -1)
        {
            // reference to current body
            immersedBody& cBody(immersedBodies_[possibleContactNeighbourList_[possiblePair].first()]);
            immersedBody& tBody(immersedBodies_[possibleContactNeighbourList_[possiblePair].second()]);
            
            List<DynamicLabelList>    commonCells;
            commonCells.setSize(Pstream::nProcs());

            List<DynamicLabelList> cSurfCells(cBody.getSurfaceCellList());
            List<DynamicLabelList> cIntCells(cBody.getInternalCellList());
            List<DynamicLabelList> tSurfCells(tBody.getSurfaceCellList());
            
            //Bounding box of body t is used only for optimization of for loop
            boundBox tBox(tBody.getMinBoundPoint(),tBody.getMaxBoundPoint());
            // Inflate tBox to include cell center of surfCells
            tBox.inflate(2*sqrt(mesh_.magSf()[0])/tBox.mag());
            
            //Iterate over surfCells to find commen cells 
            forAll (cSurfCells[Pstream::myProcNo()],cSCellI)
            {
                if (tBox.contains(mesh_.C()[cSurfCells[Pstream::myProcNo()][cSCellI]]))
                {
                    forAll (tSurfCells[Pstream::myProcNo()],tSCellI)
                    {
                        if (mag(cSurfCells[Pstream::myProcNo()][cSCellI]-tSurfCells[Pstream::myProcNo()][tSCellI]) < SMALL)
                        {
                            commonCells[Pstream::myProcNo()].append(cSurfCells[Pstream::myProcNo()][cSCellI]);                    
                        }
                    }
                }
            }
            
            scalar intersectedVolume(0);
                
            if (commonCells[Pstream::myProcNo()].size() > SMALL)
            {
                cBody.setRecomputeProjection(true);
                tBody.setRecomputeProjection(true);

                //Create triSurfaceSearch
                const triSurface& tibTemp( tBody.getTriSurfMesh());
                triSurfaceSearch tibTriSurfSearch( tibTemp );
                const triSurface& cibTemp( cBody.getTriSurfMesh());
                triSurfaceSearch cibTriSurfSearch( cibTemp );
                const pointField& pp = mesh_.points();
                
                boolList tcenterInsideList = tibTriSurfSearch.calcInside( mesh_.C());
                boolList ccenterInsideList = cibTriSurfSearch.calcInside( mesh_.C());

                // iterate over all surfCells and intCells to evaluate intersected volume
                forAll (cSurfCells[Pstream::myProcNo()],cSCellI)
                {
                    if (tBox.contains(mesh_.C()[cSurfCells[Pstream::myProcNo()][cSCellI]]))
                    {
                        const labelList& vertexLabels = mesh_.cellPoints()[cSurfCells[Pstream::myProcNo()][cSCellI]];
                        const pointField vertexPoints(pp,vertexLabels);
                        boolList tvertexesInside = tibTriSurfSearch.calcInside( vertexPoints );
                        boolList cvertexesInside = cibTriSurfSearch.calcInside( vertexPoints );
                        bool tcenterInside(tcenterInsideList[cSurfCells[Pstream::myProcNo()][cSCellI]] );
                        bool ccenterInside(ccenterInsideList[cSurfCells[Pstream::myProcNo()][cSCellI]] );
                        scalar rVInSize(1.0/(tvertexesInside.size()+1));
                        // Note: weiggetContactNormalht of a single vertex in the cell
                        
                        scalar partialVolume(0);
                        forAll (tvertexesInside, verIn)
                        {
                            if (tvertexesInside[verIn]==true && cvertexesInside[verIn]==true)
                            {
                                partialVolume += rVInSize; //fraction of cell covered
                            }
                        }

                        if (tcenterInside==true && ccenterInside==true)
                        {
                            partialVolume += rVInSize; //fraction of cell covered
                        }
                        
                        intersectedVolume += mesh_.V()[cSurfCells[Pstream::myProcNo()][cSCellI]] * partialVolume;
                        
                        if (intersectedVolume > SMALL) break;
                    }
                }

                forAll (cIntCells[Pstream::myProcNo()],cSCellI)
                {
                    if (tBox.contains(mesh_.C()[cIntCells[Pstream::myProcNo()][cSCellI]]))
                    {
                        const labelList& vertexLabels = mesh_.cellPoints()[cIntCells[Pstream::myProcNo()][cSCellI]];
                        const pointField vertexPoints(pp,vertexLabels);
                        boolList tvertexesInside = tibTriSurfSearch.calcInside( vertexPoints );
                        boolList cvertexesInside = cibTriSurfSearch.calcInside( vertexPoints );
                        bool tcenterInside(tcenterInsideList[cIntCells[Pstream::myProcNo()][cSCellI]] );
                        bool ccenterInside(ccenterInsideList[cIntCells[Pstream::myProcNo()][cSCellI]] );
                        scalar rVInSize(1.0/(tvertexesInside.size()+1));
                        // Note: weight of a single vertex in the cell
                        
                        scalar partialVolume(0);
                        forAll (tvertexesInside, verIn)
                        {
                            if (tvertexesInside[verIn]==true && cvertexesInside[verIn]==true)
                            {
                                partialVolume += rVInSize; //fraction of cell covered
                            }
                        }
                        
                        if (tcenterInside==true && ccenterInside==true)
                        {
                            partialVolume += rVInSize; //fraction of cell covered
                        }
                        
                        intersectedVolume += mesh_.V()[cIntCells[Pstream::myProcNo()][cSCellI]] * partialVolume;
                        
                        if (intersectedVolume > SMALL) break;
                    }
                }
            }
            
            reduce(intersectedVolume, sumOp<scalar>());

            if (intersectedVolume > SMALL)
            {                
                prtContactIBList_.append(possibleContactNeighbourList_[possiblePair]);
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
