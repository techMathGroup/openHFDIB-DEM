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
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*)
\*---------------------------------------------------------------------------*/
#include "immersedBody.H"
#include "polyMesh.H"
#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"


#include "processorPolyPatch.H"                                         //OF.com: required



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
#include "solverInfo.H" 

#include "interpolationTable.H"

#define ORDER 2

using namespace Foam;

//---------------------------------------------------------------------------//
immersedBody::immersedBody
(
    word bodyName,
    const Foam::fvMesh& mesh,
    dictionary& HFDIBDEMDict,
    dictionary& transportProperties,
    label bodyId,
    label recomputeM0,
    std::shared_ptr<geomModel> bodyGeomModel,
    autoPtr<ibInterpolation>& ibIntp,
    List<labelList>& cellPoints
)
:
bodyName_(bodyName),
isActive_(true),
immersedDict_(HFDIBDEMDict.subDict(bodyName_)),
mesh_(mesh),
transportProperties_(transportProperties),
geomModel_(std::move(bodyGeomModel)),
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
couplingHistCoef_(1.0),
CoNum_(0.0),
rhoF_(1.0),
bodyId_(bodyId),
bodyIdStr_(Foam::name(bodyId_)),
updateTorque_(false),
bodyOperation_(0),
octreeField_(mesh_.nCells(), 0),
cellToStartInCreateIB_(0),
startSynced_(false),
totRotMatrix_(tensor::I),
sdBasedLambda_(false),
overlayLambda_("add"),
intSpan_(2.0),
charCellSize_(1e3),
refineBuffers_(recomputeM0 > 0 ? recomputeM0 + 1 : -1),
recomputeM0_(recomputeM0),
timesToSetStatic_(-1),
staticContactPost_(vector::zero)
{
    #include "initializeIB.H"

    InfoH << iB_Info << "Finished body initialization" << endl;
    InfoH << basic_Info << "New bodyID: " << bodyIdStr_ << " name: "
        << bodyName_ << " rhoS: " << geomModel_->getRhoS()
        << " dC: " << getDC() << endl;
}
//---------------------------------------------------------------------------//
immersedBody::~immersedBody()
{
}
//---------------------------------------------------------------------------//
// Create immersed body info
void immersedBody::createImmersedBody
(
    volScalarField& body,
    volScalarField& refineF,
    bool synchCreation
)
{
    geomModel_->createImmersedBody(
        body,
        octreeField_,
        cellPoints_
    );
    
    if(synchCreation)
    {
        syncCreateImmersedBody(body, refineF);
    }

    computeCharCellSize();                                              //used in intpInfo_->setIntpInfo()
    // Note (MI): in theory, it should be enough to compute body
    //            characteristic cell size only once - after the first
    //            creation on a sufficiently refined mesh
    // => we should look into this in future
    // Note (MI): computeCharCellSize() has gMax in it - is it efficient?
    intpInfo_->setCharCellSize(charCellSize_);                          //set characteristic cell size to find interpolation points

    Info << "!! -- attempting to set interpolation info -- !!" << endl;
    intpInfo_->setIntpInfo();
    Info << "!! -- finished setting interpolation info -- !!" << endl;
}
//---------------------------------------------------------------------------//
void immersedBody::syncCreateImmersedBody                               //Note (MI): the name does not reflect the content
(
    volScalarField& body,
    volScalarField& refineF
)
{
    syncImmersedBodyGeometry(body, refineF);
    // Note (MI): during simplification, calculateGeometricalProperties
    //            was replaced by calculateGeometricalPropertiesParallel
    //            and this replacement WAS NOT tested (a problem might
    //            appear for clusterBodies)
    syncImmersedBodyRefinement(body, refineF);
}
//---------------------------------------------------------------------------//
void immersedBody::syncImmersedBodyGeometry
(
    volScalarField& body,
    volScalarField& refineF
)
{
    geomModel_->setOwner();
    InfoH << iB_Info << "body " << bodyIdStr_
        << " owner: " << geomModel_->getOwner() << endl;

    InfoH << iB_Info << "Computing geometrical properties" << endl;
    geomModel_->calculateGeometricalPropertiesParallel(body);
}
//---------------------------------------------------------------------------//
void immersedBody::syncImmersedBodyRefinement
(
    volScalarField& body,
    volScalarField& refineF
)
{
    // update body courant number
    // computeBodyCoNumber();

    InfoH << iB_Info << "-- body " << bodyIdStr_
        << " current center of mass position: " << geomModel_->getCoM() << endl;

    const List<DynamicLabelList>& surfCells = geomModel_->getSurfaceCellList();
    DynamicLabelList zeroList(surfCells[Pstream::myProcNo()].size(), 0);

    constructRefineField
    (
        body,
        refineF,
        surfCells[Pstream::myProcNo()],
        zeroList
    );

}
//---------------------------------------------------------------------------//
void immersedBody::computeCharCellSize()
{
    const List<DynamicLabelList>& surfCells = geomModel_->getSurfaceCellList();

    scalarList charCellSizeL(Pstream::nProcs(),1e4);
    forAll (surfCells[Pstream::myProcNo()],sCellI)
    {
        label cellI = surfCells[Pstream::myProcNo()][sCellI];
        
        scalar cellMeasure = mesh_.V()[cellI];
        label nGeometricD = mesh_.nGeometricD();

        if (!case3D)
        {
            scalar emptyThickness = 1.0;
            forAll(emptyDir,dirI)                                       //this is based on settings from HFDIBDEMDict, not actual mesh
            {                                                           //plus: adaptively refined meshes can be treated as 2D
                if (emptyDir[dirI])                                     //minus: forces user to check both mesh and HFDIBDEMDict
                {
                    emptyThickness *= mesh_.bounds().span()[dirI];
                    if (refineBuffers_ > 1)
                    {
                        emptyThickness /= (refineBuffers_ - 1);
                    }
                    nGeometricD--;
                }
            }     
            cellMeasure /= emptyThickness;
        }

        charCellSizeL[Pstream::myProcNo()] =
            min
            (
                charCellSizeL[Pstream::myProcNo()],
                Foam::pow(cellMeasure,1.0/nGeometricD)
            );
    }
    forAll(charCellSizeL,indl)
    {
        if(charCellSizeL[indl] > 5e3)
        {
            charCellSizeL[indl] = -1.0;
        }
    }

    charCellSize_ = gMax(charCellSizeL);                                //I want to be sure to always go to another cells
    InfoH << iB_Info << "-- body " << bodyIdStr_
        << " characteristic cell size: " << charCellSize_ << endl;
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

    List<DynamicLabelList> facesToSendToProcs;
    facesToSendToProcs.setSize(Pstream::nProcs());

    for(label i = 0; i < refineBuffers_; i++)
    {
        // filter only cells with current level
        forAll(cellsToIterate, cellI)
        {
            if(startLevel[cellI] == i)
            {
                cellsToIterateC.append(cellsToIterate[cellI]);
                refineF[cellsToIterate[cellI]] = 1;
            }
        }

        // iterate over cells with current level and search for neighbors
        forAll(cellsToIterateC, cellI)
        {
            // get cell faces
            labelList cellFaces(mesh_.cells()[cellsToIterateC[cellI]]);
            forAll(cellFaces, faceI)
            {
                // for internal faces assign next level to neighbors and save
                if (mesh_.isInternalFace(cellFaces[faceI]))
                {
                    // get cell label
                    label nCell(mesh_.owner()[cellFaces[faceI]]);

                    // switch from owner to neighbor if needed
                    if(nCell == cellsToIterateC[cellI])
                    {
                        nCell = mesh_.neighbour()[cellFaces[faceI]];
                    }

                    // if not included already
                    if(refineF[nCell] == 0)
                    {
                        // in futher iterations add only cells outside of body
                        if(i > 0)
                        {
                            if(body[nCell] < SMALL)
                            {
                                refineF[nCell] = 1;
                                cellsToIterateF.append(nCell);
                            }
                        }

                        // in first iteration add any neighbor
                        else
                        {
                            refineF[nCell] = 1;
                            cellsToIterateF.append(nCell);
                        }
                    }
                }

                // check for processor neighbors
                else
                {
                    label facePatchId(mesh_.boundaryMesh().whichPatch(
                        cellFaces[faceI]
                    ));

                    const polyPatch& cPatch = mesh_.boundaryMesh()[facePatchId];
                    if (cPatch.type() == "processor")
                    {
                        const processorPolyPatch& procPatch =
                            refCast<const processorPolyPatch>(cPatch);
                        if (procPatch.myProcNo() == Pstream::myProcNo())
                        {
                            // send patch local face id
                            facesToSendToProcs[procPatch.neighbProcNo()].append(
                                cPatch.whichFace(cellFaces[faceI])
                            );
                        }
                        else
                        {
                            facesToSendToProcs[procPatch.myProcNo()].append(
                                cPatch.whichFace(cellFaces[faceI])
                            );
                        }
                    }
                }
            }
        }

        // send faces to other procs 
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            if (proci != Pstream::myProcNo())
            {
                UOPstream send(proci, pBufs);
                send << facesToSendToProcs[proci];
                facesToSendToProcs[proci].clear();
            }
        }
        pBufs.finishedSends();

        // recieve faces from other procs
        List<DynamicLabelList> facesReceivedFromProcs(Pstream::nProcs());
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            if (proci != Pstream::myProcNo())
            {
                UIPstream recv(proci, pBufs);
                DynamicLabelList recList (recv);
                facesReceivedFromProcs[proci] = recList;
            }
        }
        pBufs.clear();

        // convert to cell labels
        List<DynamicLabelList> cellLabelRecv(Pstream::nProcs());
        forAll (mesh_.boundaryMesh(), patchi)
        {
            const polyPatch& cPatch = mesh_.boundaryMesh()[patchi];
            if (cPatch.type() == "processor")
            {
                const processorPolyPatch& procPatch
                    = refCast<const processorPolyPatch>(cPatch);

                label sProc = (Pstream::myProcNo() == procPatch.myProcNo())
                    ? procPatch.neighbProcNo() : procPatch.myProcNo();

                cellLabelRecv[sProc].setSize(facesReceivedFromProcs[sProc].size());
                forAll(facesReceivedFromProcs[sProc], faceI)
                {
                    cellLabelRecv[sProc][faceI]
                        = mesh_.faceOwner()[cPatch.start()
                        + facesReceivedFromProcs[sProc][faceI]];
                }
            }
        }

        forAll(cellLabelRecv, proci)
        {
            if (proci == Pstream::myProcNo())
            {
                continue;
            }

            // loop over received cells
            forAll(cellLabelRecv[proci], cellI)
            {
                // get cell label
                label nCell(cellLabelRecv[proci][cellI]);

                // if not included already
                if(refineF[nCell] == 0)
                {
                    // in futher iterations add only cells outside of body
                    if(i > 0)
                    {
                        if(body[nCell] < SMALL)
                        {
                            refineF[nCell] = 1;
                            cellsToIterateF.append(nCell);
                        }
                    }

                    // in first iteration add any neighbor
                    else
                    {
                        refineF[nCell] = 1;
                        cellsToIterateF.append(nCell);
                    }
                }
            }
        }

        // prepare for next iteration
        cellsToIterateC = cellsToIterateF;
        cellsToIterateF.clear();
    }
}
//---------------------------------------------------------------------------//
void immersedBody::postPimpleUpdateImmersedBody
(
    const volScalarField& body,
    const volVectorField& fPress,
    const volVectorField& fVisc,
    const bool applyAddedMass
)
{
    if(!solverInfo::getOnlyDEM())
    {
        updateCoupling(body, fPress, fVisc, applyAddedMass);
    }
    resetPostPimpleState();
}
//---------------------------------------------------------------------------//
void immersedBody::postPimpleUpdateImmersedBody
(
    const volScalarField& body,
    const volVectorField& f,
    const bool kinematicForce,
    const bool applyAddedMass
)
{
    if(!solverInfo::getOnlyDEM())
    {
        updateCoupling(body, f, kinematicForce, applyAddedMass);
    }
    resetPostPimpleState();
}
//---------------------------------------------------------------------------//
void immersedBody::postPimpleUpdateImmersedBody
(
    const volScalarField& body,
    const volVectorField& f,
    const volScalarField& rho,
    const bool kinematicForce,
    const bool applyAddedMass
)
{
    if(!solverInfo::getOnlyDEM())
    {
        updateCoupling(body, f, rho, kinematicForce, applyAddedMass);
    }
    resetPostPimpleState();
}
//---------------------------------------------------------------------------//
void immersedBody::resetPostPimpleState()
{
    Vel_ = VelOld_;
    Axis_ = AxisOld_;
    omega_ = omegaOld_;
    FCouplingOld_ = FCoupling_;
}
//---------------------------------------------------------------------------//
void immersedBody::updateCoupling
(
    const volScalarField& body,
    const volVectorField& fPress,
    const volVectorField& fVisc,
    const bool applyAddedMass
)
{
    vector FV(vector::zero);
    vector TA(vector::zero);
    vector FAdded(vector::zero);
    
    // calcualate viscous force and torque
    List<DynamicLabelList> intLists;
    List<DynamicLabelList> surfLists;
    List<DynamicLabelList> haloLists;
    DynamicVectorList refCoMList;

    geomModel_->getReferencedLists(
        intLists,
        surfLists,
        haloLists,
        refCoMList
    );

    forAll (intLists, i)
    {
        DynamicLabelList& intListI = intLists[i];
        forAll (intListI, intCell)
        {
            label cellI = intListI[intCell];

            vector fCellPress = fPress[cellI];

            FV -= (fCellPress)*mesh_.V()[cellI];
            FAdded -= (fPress.prevIter()[cellI] - fPress[cellI])
                *mesh_.V()[cellI];
            TA -= ((mesh_.C()[cellI] - refCoMList[i])^fCellPress)
                *mesh_.V()[cellI];
        }
    }

    const List<point>& ibPoints = intpInfo_->getIbPoints();             //get surface points
    forAll (surfLists, i)
    {
        DynamicLabelList& surfListI = surfLists[i];
        forAll (surfListI, surfCell)
        {
            label cellI = surfListI[surfCell];

            vector fCellVisc = body[cellI]*fVisc[cellI];
            vector fCellPress = body[cellI]*fPress[cellI];

            FV -= (fCellVisc + fCellPress)*mesh_.V()[cellI];
            TA -= ((ibPoints[i] - refCoMList[i])^(fCellVisc + fCellPress))
                *mesh_.V()[cellI];
            FAdded -= body[cellI]
                *((fVisc.prevIter()[cellI] - fVisc[cellI])
                + (fPress.prevIter()[cellI] - fPress[cellI]))
                *mesh_.V()[cellI];
        }
    }
    
    reduce(FV, sumOp<vector>());
    reduce(TA, sumOp<vector>());
    reduce(FAdded, sumOp<vector>());

    FV *= rhoF_.value();
    TA *= rhoF_.value();
    FAdded *= rhoF_.value();

    FAdded = FCouplingOld_.F - FV;
    
    FCoupling_ = couplingHistCoef_*forces(FV, TA) + (1.0-couplingHistCoef_)*FCouplingOld_;

    applyAddedMassScaling(FV, FAdded, applyAddedMass);

    couplingHistCoef_ = max(couplingHistCoef_*0.95, 0.5);
    
    InfoH << iB_Info << "-- body " << bodyIdStr_ << " coupling coefficient: " << couplingHistCoef_ << endl;
}
//---------------------------------------------------------------------------//
void immersedBody::updateCoupling                                       //full interface
(
    const volScalarField& body,
    const volVectorField& f,
    const bool kinematicForce,
    const bool applyAddedMass
)
{
    vector FV(vector::zero);
    vector TA(vector::zero);
    vector FAdded(vector::zero);


    List<DynamicLabelList> intLists;
    List<DynamicLabelList> surfLists;
    List<DynamicLabelList> haloLists;
    DynamicVectorList refCoMList;

    geomModel_->getReferencedLists(
        intLists,
        surfLists,
        haloLists,
        refCoMList
    );

    // forAll (intLists, i)
    // {
    //     DynamicLabelList& intListI = intLists[i];
    //     forAll (intListI, intCell)
    //     {
    //         label cellI = intListI[intCell];

    //         FV -=  f[cellI]*mesh_.V()[cellI];
    //         TA -=  ((mesh_.C()[cellI] - refCoMList[i])^f[cellI])
    //             *mesh_.V()[cellI];
    //         FAdded -= (f.prevIter()[cellI] - f[cellI])*mesh_.V()[cellI];
    //     }
    // }

    const List<point>& ibPoints = intpInfo_->getIbPoints();             //get surface points
    forAll (surfLists, i)
    {
        DynamicLabelList& surfListI = surfLists[i];
        forAll (surfListI, surfCell)
        {
            label cellI = surfListI[surfCell];

            // scalar fScale = 1.0*body[cellI]+0.5;
            // vector fCell = (1.0 - body[cellI])*f[cellI];
            scalar fScale = 1.0;
            vector fCell =  f[cellI]*mesh_.V()[cellI];
            fCell *= fScale; 

            // FV -=  fScale*f[cellI]*mesh_.V()[cellI];
            // TA -=  ((ibPoints[surfCell] - refCoMList[i])^(fScale*f[cellI])
            //     *mesh_.V()[cellI]);
            // FAdded -= (f.prevIter()[cellI] - f[cellI])*mesh_.V()[cellI];//under construction

            FV -=  fCell;
            TA -=  (ibPoints[surfCell] - refCoMList[i])^fCell;
            FAdded -= (f.prevIter()[cellI] - f[cellI])*mesh_.V()[cellI];//under construction
        }
    }

    reduce(FV, sumOp<vector>());
    reduce(TA, sumOp<vector>());
    reduce(FAdded, sumOp<vector>());

    if (kinematicForce)
    {
        FV *= rhoF_.value();
        TA *= rhoF_.value();
        FAdded *= rhoF_.value();
    }

    scalar rhoS = geomModel_->getRhoS().value();
    // FV /= rhoS;
    // TA /= rhoS;
    // FAdded /= rhoS;
    // FV *= rhoF_.value();
    // TA *= rhoF_.value();
    // FAdded *= rhoF_.value();

    // FAdded = FCouplingOld_.F - FV;

    FAdded *= rhoF_.value()/rhoS;

    FCoupling_ = couplingHistCoef_*forces(FV, TA) + (1.0-couplingHistCoef_)*FCouplingOld_;

    applyAddedMassScaling(FV, FAdded, applyAddedMass);

    couplingHistCoef_ = max(couplingHistCoef_*0.95, 0.5);
    
    InfoH << iB_Info << "-- body " << bodyIdStr_ << " coupling coefficient: " << couplingHistCoef_ << endl;
}
//---------------------------------------------------------------------------//
void immersedBody::updateCoupling                                       //full interface
(
    const volScalarField& body,
    const volVectorField& f,
    const volScalarField& rho,
    const bool kinematicForce,
    const bool applyAddedMass
)
{
    vector FV(vector::zero);
    vector TA(vector::zero);
    vector FAdded(vector::zero);


    List<DynamicLabelList> intLists;
    List<DynamicLabelList> surfLists;
    List<DynamicLabelList> haloLists;
    DynamicVectorList refCoMList;

    geomModel_->getReferencedLists(
        intLists,
        surfLists,
        haloLists,
        refCoMList
    );

    List<List<DynamicLabelList>>& surfToHaloAddressing = geomModel_->getSurfToHaloLabels();

    forAll (intLists, i)
    {
        DynamicLabelList& intListI = intLists[i];
        forAll (intListI, intCell)
        {
            label cellI = intListI[intCell];

            scalar fScale = rhoF_.value()/rho[cellI];
            vector fCell =  f[cellI]*mesh_.V()[cellI];
            fCell *= fScale;

            FV -=  fCell;
            TA -=  ((mesh_.C()[cellI] - refCoMList[i])^fCell);
            FAdded -= (f.prevIter()[cellI] - f[cellI])*mesh_.V()[cellI];
        }
    }

    const List<point>& ibPoints = intpInfo_->getIbPoints();             //get surface points
    forAll (surfLists, i)
    {
        DynamicLabelList& surfListI = surfLists[i];
        DynamicLabelList& haloListI = haloLists[i];
        forAll (surfListI, surfCell)
        {
            label cellI = surfListI[surfCell];
            DynamicLabelList& surfToHalo = surfToHaloAddressing[Pstream::myProcNo()][surfCell];

            scalar fluidMass(0);
            scalar fluidVol(0);
            forAll (surfToHalo, haloCell)
            {
                label sToHI = surfToHalo[haloCell];
                label cellH = haloListI[sToHI];
                fluidMass += rho[cellH]*mesh_.V()[cellH]*(1.0 - body[cellH]);
                fluidVol += mesh_.V()[cellH]*(1.0 - body[cellH]);
            }
            scalar surfFluidRho = fluidMass/(fluidVol + SMALL);

            // scalar fScale = 1.0*body[cellI]+0.5;
            // vector fCell = (1.0 - body[cellI])*f[cellI];
            scalar fScale = surfFluidRho/rho[cellI];
            vector fCell =  f[cellI]*mesh_.V()[cellI];
            fCell *= fScale; 

            // FV -=  fScale*f[cellI]*mesh_.V()[cellI];
            // TA -=  ((ibPoints[surfCell] - refCoMList[i])^(fScale*f[cellI])
            //     *mesh_.V()[cellI]);
            // FAdded -= (f.prevIter()[cellI] - f[cellI])*mesh_.V()[cellI];//under construction

            FV -=  fCell;
            TA -=  (ibPoints[surfCell] - refCoMList[i])^fCell;
            FAdded -= (f.prevIter()[cellI] - f[cellI])*mesh_.V()[cellI];//under construction
        }
    }

    reduce(FV, sumOp<vector>());
    reduce(TA, sumOp<vector>());
    reduce(FAdded, sumOp<vector>());

    if (kinematicForce)
    {
        FV *= rhoF_.value();
        TA *= rhoF_.value();
        FAdded *= rhoF_.value();
    }

    scalar rhoS = geomModel_->getRhoS().value();
    // FV /= rhoS;
    // TA /= rhoS;
    // FAdded /= rhoS;
    // FV *= rhoF_.value();
    // TA *= rhoF_.value();
    // FAdded *= rhoF_.value();

    // FAdded = FCouplingOld_.F - FV;

    FAdded *= rhoF_.value()/rhoS;

    FCoupling_ = couplingHistCoef_*forces(FV, TA) + (1.0-couplingHistCoef_)*FCouplingOld_;

    applyAddedMassScaling(FV, FAdded, applyAddedMass);

    couplingHistCoef_ = max(couplingHistCoef_*0.95, 0.5);
    
    InfoH << iB_Info << "-- body " << bodyIdStr_ << " coupling coefficient: " << couplingHistCoef_ << endl;
}
//---------------------------------------------------------------------------//
void immersedBody::applyAddedMassScaling
(
    const vector& FV,
    const vector& FAdded,
    const bool applyAddedMass
)
{
    if (!applyAddedMass)
    {
        return;
    }

    const scalar m0 = geomModel_->getM0();
    const scalar massSign = ((FV & FAdded) < 0.0) ? -1.0 : 1.0;
    // scalar massAdded = min(1.0*m0, mag(FAdded)/(mag(a_) + SMALL));
    scalar massAdded = mag(FAdded)/(mag(a_) + SMALL);
    massAdded *= massSign;
    InfoH << iB_Info << "-- body " << bodyIdStr_ << " massAdded: " << massAdded
        << " m0: " << m0 << endl;
    InfoH << iB_Info << "-- body " << bodyIdStr_ << " orig coupling force: " << FCoupling_.F << " orig coupling torque: " << FCoupling_.T << endl;
    const scalar scale = (m0 + massAdded)/m0;
    FCoupling_.F *= scale;
    FCoupling_.T *= scale;
    InfoH << iB_Info << "-- body " << bodyIdStr_ << " scld coupling force: " << FCoupling_.F << " scld coupling torque: " << FCoupling_.T << endl;
}
//---------------------------------------------------------------------------//
void immersedBody::updateLocalFluidDensity
(
    const volScalarField& body,
    volScalarField& rho
)
{
    List<DynamicLabelList> intLists;
    List<DynamicLabelList> surfLists;
    List<DynamicLabelList> haloLists;
    DynamicVectorList refCoMList;

    geomModel_->getReferencedLists(
        intLists,
        surfLists,
        haloLists,
        refCoMList
    );

    scalar rhoS = geomModel_->getRhoS().value();

    forAll (intLists, i)
    {
        DynamicLabelList& intListI = intLists[i];
        forAll (intListI, intCell)
        {
            label cellI = intListI[intCell];
            rho[cellI]  = rhoS;
        }
    }

    forAll (surfLists, i)
    {
        DynamicLabelList& surfListI = surfLists[i];
        forAll (surfListI, surfCell)
        {
            label cellI = surfListI[surfCell];
            rho[cellI]  = body[cellI] * rhoS + (1.0 - body[cellI]) * rho[cellI];
        }
    }
}
void immersedBody::updateLocalFluidDensity
(
    const volScalarField& body,
    volScalarField& rho,
    volScalarField& rhoS
)
{
    List<DynamicLabelList> intLists;
    List<DynamicLabelList> surfLists;
    List<DynamicLabelList> haloLists;
    DynamicVectorList refCoMList;

    geomModel_->getReferencedLists(
        intLists,
        surfLists,
        haloLists,
        refCoMList
    );

    scalar rhoSVal = geomModel_->getRhoS().value();

    forAll (intLists, i)
    {
        DynamicLabelList& intListI = intLists[i];
        forAll (intListI, intCell)
        {
            label cellI = intListI[intCell];
            rho[cellI]  = rhoSVal;
            rhoS[cellI] = rhoSVal;
        }
    }

    forAll (surfLists, i)
    {
        DynamicLabelList& surfListI = surfLists[i];
        forAll (surfListI, surfCell)
        {
            label cellI = surfListI[surfCell];
            rho[cellI]  = body[cellI] * rhoSVal + (1.0 - body[cellI]) * rho[cellI];
            rhoS[cellI] = rhoSVal;
        }
    }
}
//---------------------------------------------------------------------------//
// update movement variables of the body
void immersedBody::updateMovement
(
    scalar deltaT
)
{
    updateMovementComp(deltaT,Vel_,Axis_,omega_);
}
void immersedBody::updateMovement
(
    vector Vel,
    vector Axis,
    scalar omega
)
{
    scalar deltaT = mesh_.time().deltaT().value();
    updateMovementComp(deltaT,Vel,Axis,omega);
}
void immersedBody::updateMovementComp
(
    scalar deltaT,
    vector Vel,
    vector Axis,
    scalar omega
)
{
    auto updateTranslation = [&]()
    {

        const uniformDimensionedVectorField& g =
            mesh_.lookupObject<uniformDimensionedVectorField>("g");
            
        vector FG(vector::zero);
        if(!solverInfo::getOnlyDEM())
            FG = geomModel_->getM0()*(1.0-rhoF_.value()
            /geomModel_->getRhoS().value())*g.value();
        else
            FG = geomModel_->getM0()*g.value();

        vector F(FCoupling_.F);
        F += FContact_.F;
        F += FG;

        if(!case3D)
        {
            F[emptyDim] *= 0;
            FG[emptyDim] *= 0;
        }
        if(geomModel_->getM0() > 0)
        {
            // compute current acceleration (assume constant over timeStep)

            InfoH << iB_Info <<"-- body "<< bodyIdStr_ <<" mass            : " << geomModel_->getM0() << endl;
            InfoH << iB_Info <<"-- body "<< bodyIdStr_ <<" acting force    : " << F << endl;
            InfoH << iB_Info <<"-- body "<< bodyIdStr_ <<" coupling force  : " << FCoupling_.F << endl;
            InfoH << iB_Info <<"-- body "<< bodyIdStr_ <<" grav/buyo force : " << FG << endl;
            
            a_  = F/(geomModel_->getM0());
            
            // update body linear velocity
            Vel_ = Vel + deltaT*a_;
            InfoH << iB_Info <<"-- body "<< bodyIdStr_ <<" accelaration    : " << a_ << endl;
        }
    };

    auto updateRotation = [&]()
    {
        if(mag(geomModel_->getI()) > 0)
        {
            vector T(FCoupling_.T);
            T += FContact_.T;

            // update body angular acceleration
            alpha_ = inv(geomModel_->getI()) & T;
            // update body angular velocity
            vector Omega(Axis*omega + deltaT*alpha_);
            // split Omega into Axis_ and omega_
            omega_ = mag(Omega);

            if (omega_ < SMALL)
            {
                Axis_ = vector::one;
                if (!case3D)
                {
                    const vector validDirs = (geometricD + vector::one)/2;
                    Axis_ -= validDirs;
                }
            }
            else
            {
                Axis_ =  Omega/(omega_+SMALL);
                if (!case3D)
                {// in 2D, I need to keep only the part of the rotation axis
                    const vector validDirs = (geometricD + vector::one)/2;
                    Axis_ = cmptMultiply(vector::one-validDirs,Axis_);
                }
            }
            Axis_ /= mag(Axis_);
        }
    };

    auto updateRotationFixedAxis = [&]()
    {
        vector T(FCoupling_.T);
        T += FContact_.T;

        // update body angular velocity
        vector Omega(Axis*omega + deltaT * (inv(geomModel_->getI()) & T));

        // split Omega into Axis_ and omega_
        omega_ = mag(Omega);

        vector newAxis = Omega/(omega_+SMALL);
        if ((newAxis & Axis_) < 0) Axis_ *= (-1.0);;
    };

    auto updatePositionByTable = [&]()
    {
        vector F(FCoupling_.F);
        F *= 0.0; // no force

        scalar time = mesh_.time().value();

        dictionary functionDict = immersedDict_.subDict("prescribedPosTableBody");
        dictionary posIntTableDict = functionDict.subDict("posIntTableDict");
        interpolationTable<vector> posIntTable = interpolationTable<vector>(posIntTableDict);

        if(geomModel_->getM0() > 0)
        {
            vector posNew = posIntTable(time);
            vector posOld = posIntTable(time-deltaT);

            // update velocity
            Vel_ = (posNew - posOld)/deltaT;
        }
    };

    if (bodyOperation_ == 0 or bodyOperation_ == 3)
    {
        return;
    }
    else if (bodyOperation_ == 1)
    {
        updateRotation();
        return;
    }
    else if (bodyOperation_ == 2)
    {
        updateTranslation();
        return;
    }
    else if (bodyOperation_ == 4)
    {
        updateRotationFixedAxis();
        return;
    }
    else if (bodyOperation_ == 6)
    {
        updatePositionByTable();
        return;
    }

    updateTranslation();
    if (updateTorque_)
    {
        updateRotation();
    }

    return;
}
//---------------------------------------------------------------------------//
// move immersed body according to body operation
void immersedBody::moveImmersedBody
(
    scalar deltaT
)
{
    if (bodyOperation_ == 0) return;

    if (mag(deltaT + 1.0) < SMALL) deltaT = mesh_.time().deltaT().value();

    // incremental rotation angle
    scalar angle     = omega_*deltaT;

    // translation increment
    vector transIncr = Vel_*deltaT;

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
    scalar sy = Foam::sqrt(totRotMatrix_.xx()*totRotMatrix_.xx()
        + totRotMatrix_.yy()*totRotMatrix_.yy());

    if (sy > SMALL)
    {
        eulerAngles.x() =
            Foam::atan2(totRotMatrix_.zy(),totRotMatrix_.zz());
        eulerAngles.y() = Foam::atan2(-totRotMatrix_.zx(),sy);
        eulerAngles.z() =
            Foam::atan2(totRotMatrix_.yx(),totRotMatrix_.xx());
    }
    else
    {
        eulerAngles.x() =
            Foam::atan2(-totRotMatrix_.yz(),totRotMatrix_.yy());
        eulerAngles.y() = Foam::atan2(-totRotMatrix_.zx(),sy);
        eulerAngles.z() = 0.0;
    }

    geomModel_->bodyRotatePoints(angle,Axis_);
    geomModel_->bodyMovePoints(transIncr);

    // geomModel_->synchronPos();

    // InfoH << iB_Info;
    // InfoH << "-- body " << bodyId_ << " CoM                  : "
    //     << geomModel_->getCoM() << endl;
    // InfoH << "-- body " << bodyId_ << " linear velocity      : "
    //     << Vel_ << endl;
    // InfoH << "-- body " << bodyId_ << " angluar velocity     : "
    //     << omega_ << endl;
    // InfoH << "-- body " << bodyId_ << " axis of rotation     : "
    //     << Axis_ << endl;
    // InfoH << "-- body " << bodyId_ << " total rotation matrix: "
    //     << totRotMatrix_ << endl;
}
//---------------------------------------------------------------------------//
void immersedBody::printBodyInfo()
{
    InfoH << iB_Info;
    InfoH << "-- body " << bodyIdStr_ << " CoM                  : "
        << geomModel_->getCoM() << endl;
    InfoH << "-- body " << bodyIdStr_ << " linear velocity      : "
        << Vel_ << endl;
    InfoH << "-- body " << bodyIdStr_ << " angluar velocity     : "
        << omega_ << endl;
    InfoH << "-- body " << bodyIdStr_ << " axis of rotation     : "
        << Axis_ << endl;
    InfoH << "-- body " << bodyIdStr_ << " total rotation matrix: "
        << totRotMatrix_ << endl;
}
//---------------------------------------------------------------------------//
void immersedBody::updateVectorField
(
    volVectorField& VS,
    word VName,
    volScalarField& body
)
{
    // check dictionary for parameters (only noSlip allowed)
    word BC = word(immersedDict_.subDict(VName).lookup("BC"));

    List<DynamicLabelList> intLists;
    List<DynamicLabelList> surfLists;
    List<DynamicLabelList> haloLists;
    DynamicVectorList refCoMList;

    geomModel_->getReferencedLists(
        intLists,
        surfLists,
        haloLists,
        refCoMList
    );

    if (BC=="noSlip")
    {
        // if STATICBODY set to zero
        if ( bodyOperation_==0)
        {
            forAll (intLists, i)
            {
                DynamicLabelList& intListI = intLists[i];
                forAll (intListI, intCell)
                {
                    label cellI = intListI[intCell];
                    VS[cellI]   = Vel_;
                }
            }

            forAll (surfLists, i)
            {
                DynamicLabelList& surfListI = surfLists[i];
                forAll (surfListI, surfCell)
                {
                    label cellI = surfListI[surfCell];
                    VS[cellI]   = Vel_;
                }
            }
        }
        else
        {
            forAll (intLists, i)
            {
                DynamicLabelList& intListI = intLists[i];
                forAll (intListI, intCell)
                {
                    label cellI = intListI[intCell];
                    vector planarVec =  mesh_.C()[cellI] - refCoMList[i]
                                    - Axis_*(
                                    (mesh_.C()[cellI] - refCoMList[i])
                                    &Axis_);

                    vector VSvalue = (-(planarVec^Axis_)*omega_ + Vel_);
                    VS[cellI] = VSvalue;
                }
            }

            forAll (surfLists, i)
            {
                DynamicLabelList& surfListI = surfLists[i];
                forAll (surfListI, surfCell)
                {
                    label cellI = surfListI[surfCell];

                    vector planarVec =  mesh_.C()[cellI] - refCoMList[i]
                                    - Axis_*(
                                    (mesh_.C()[cellI] - refCoMList[i])
                                    &Axis_);

                    vector VSvalue = (-(planarVec^Axis_)*omega_ + Vel_);
                    VS[cellI] = VSvalue;
                }
            }
        }
    }
}
//---------------------------------------------------------------------------//
// reset body field for this immersed object
vectorField immersedBody::getUatIbPoints()
{
    const List<point>& ibPoints = intpInfo_->getIbPoints();
    vectorField ibPointsVal(ibPoints.size());
    forAll(ibPoints, pointI)
    {
        // vector planarVec =  geomModel_->getLVec(ibPoints[pointI])
        //                     - Axis_*(
        //                     (geomModel_->getLVec(ibPoints[pointI]))&Axis_);

        vector planarVec =  ibPoints[pointI] - geomModel_->getCoM()
                            - Axis_*(
                            (ibPoints[pointI]-geomModel_->getCoM())&Axis_);

        vector VSvalue = (-(planarVec^Axis_)*omega_ + Vel_);
        ibPointsVal[pointI] = VSvalue;
    }

    return ibPointsVal;
}
//---------------------------------------------------------------------------//
// function to reset body
void immersedBody::resetBody
(
    volScalarField& body
)
{
    geomModel_->resetBody(body);
}
//---------------------------------------------------------------------------//
// function to move the body after the contact
void immersedBody::postContactUpdateBodyField
(
    volScalarField& body,
    volScalarField& refineF
)
{
    createImmersedBody(body,refineF,false);
}
//---------------------------------------------------------------------------//
void immersedBody::recreateBodyField
(
    volScalarField& body,
    volScalarField& refineF
)
{
    octreeField_ = Field<label>(mesh_.nCells(), 0);
    geomModel_->getSurfaceCellList()[Pstream::myProcNo()].clear();
    geomModel_->getInternalCellList()[Pstream::myProcNo()].clear();
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

    const List<DynamicLabelList>& surfCells = geomModel_->getSurfaceCellList();
    forAll (surfCells[Pstream::myProcNo()],sCellI)
    {
        label cellI(surfCells[Pstream::myProcNo()][sCellI]);

        scalar dCell(Foam::pow(mesh_.V()[cellI],0.3333));
        scalar CoNumCell(VelMag*mesh_.time().deltaT().value()/dCell);

        CoNumCell+=rotCoNumB/dCell;

        CoNum_      = max(CoNum_,CoNumCell);
        meanCoNum_ += CoNumCell;
        auxCntr    += 1;
    }
    const List<DynamicLabelList>& intCells = geomModel_->getInternalCellList();
    forAll (intCells[Pstream::myProcNo()],iCellI)
    {
        label cellI(intCells[Pstream::myProcNo()][iCellI]);

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

    InfoH << iB_Info << "-- body " << bodyIdStr_
        << " Courant Number mean: " << meanCoNum_
        << " max: " << CoNum_ << endl;

}

//---------------------------------------------------------------------------//
// print out body linear and angular momentum
void immersedBody::printMomentum()
{
    vector L(geomModel_->getI()&(Axis_*omega_));
    vector p(geomModel_->getM()*Vel_);

    InfoH << iB_Info;
    InfoH << "-- body " << bodyIdStr_ << "  linear momentum:" << p
         << " magnitude: " << mag(p) <<endl;
    InfoH << "-- body " << bodyIdStr_ << " angular momentum:" << L
         << " magnitude: " << mag(L) <<endl;
}
//---------------------------------------------------------------------------//
// print out body statistics
void immersedBody::printStats()
{
    vector L(geomModel_->getI()&(Axis_*omega_));
    vector p(geomModel_->getM()*Vel_);

    InfoH << iB_Info << "-- body " << bodyIdStr_ << "  linear momentum:" << p
        << " magnitude: " << mag(p) <<endl;
    InfoH << "-- body " << bodyIdStr_ << " angular momentum:" << L
        << " magnitude: " << mag(L) <<endl;
    InfoH << basic_Info << "-- body " << bodyIdStr_ << " CoM :"
        << geomModel_->getCoM() << endl;
    InfoH << basic_Info << "-- body " << bodyIdStr_ << "  linear velocity:"
        << Vel_ << " magnitude: " << mag(Vel_) <<endl;
    InfoH << "-- body " << bodyIdStr_ << " angular velocity:" << omega_
        << " magnitude: " << mag(omega_) <<endl;
    InfoH << "-- body " << bodyIdStr_ << "    rotation axis:" << Axis_
        << " magnitude: " << mag(Axis_) <<endl;
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
    geomModel_->resetBody(body);
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
    const List<DynamicLabelList>& intCells = geomModel_->getInternalCellList();
    forAll (intCells[Pstream::myProcNo()],iCellI)
    {
        cellI   = intCells[Pstream::myProcNo()][iCellI];
        meanV  += U[cellI]*mesh_.V()[cellI];
        meanC  += curlU[cellI]*mesh_.V()[cellI];
        totVol += mesh_.V()[cellI];
    }
    reduce(meanV, sumOp<vector>());
    reduce(meanC, sumOp<vector>());
    reduce(totVol, sumOp<scalar>());
    if(totVol > 0)
    {
        Vel_ = meanV/(totVol);
        meanC/=(totVol);
    }
    vector Omega(0.5*meanC);
    if(updateTorque_)
    {
        omega_ = mag(Omega);
        if (omega_ < SMALL)
        {
            Axis_ = vector::one;
            if (mesh_.nGeometricD() < 3)
            {
                const vector validDirs =
                    (mesh_.geometricD() + Vector<label>::one)/2;
                Axis_ -= validDirs;
            }
        }
        else
        {
            Axis_ =  Omega/(omega_+SMALL);
            if (mesh_.nGeometricD() < 3)
            {
                const vector validDirs =
                    (mesh_.geometricD() + Vector<label>::one)/2;
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
    InfoH << basic_Info << "-- body " << bodyIdStr_
        << "initial movement variables:" << endl;
    printStats();
}
//---------------------------------------------------------------------------//
void immersedBody::pimpleUpdate
(
    volScalarField& body,
    volVectorField& fPress,
    volVectorField& fVisc,
    const bool applyAddedMass
)
{
    updateCoupling(body, fPress, fVisc, applyAddedMass);
    updateMovement(VelOld_, AxisOld_, omegaOld_);
}
//---------------------------------------------------------------------------//
void immersedBody::pimpleUpdate
(
    volScalarField& body,
    volVectorField& f,
    const bool kinematicForce,
    const bool applyAddedMass
)
{
    updateCoupling(body, f, kinematicForce, applyAddedMass);
    updateMovement(VelOld_, AxisOld_, omegaOld_);
}
//---------------------------------------------------------------------------//
void immersedBody::checkIfInDomain(volScalarField& body)
{
    if(geomModel_->getM0() < SMALL)
    {
        switchActiveOff(body);
        geomModel_->resetBody(body);
    }

    InfoH << iB_Info << "-- body " << bodyIdStr_ << " current M/M0: "
        << geomModel_->getM()/geomModel_->getM0() << endl;
    // if only 1% of the initial particle mass remains in the domain, switch it off
    if (geomModel_->getM()/(geomModel_->getM0()+SMALL) < 1e-2 && case3D)
    {
        switchActiveOff(body);
        geomModel_->resetBody(body);
    }
    else if (!case3D && geomModel_->getNCells() <= 1 && !geomModel_->isCluster())
    {
        switchActiveOff(body);
        geomModel_->resetBody(body);
        InfoH << iB_Info << "-- body " << bodyIdStr_ << " switched off" << endl;
    }
}
//---------------------------------------------------------------------------//
void immersedBody::setRestartSim(vector vel, scalar angVel, vector axisRot, bool setStatic, label timesInContact)
{
    Vel_ = vel;
    omega_ = angVel;
    Axis_ = axisRot;
    ibContactClass_->setTimeStepsInContWStatic(timesInContact);
    InfoH << iB_Info << "-- body " << bodyIdStr_
        << " timeStepsInContWStatic_: "
        << ibContactClass_->getTimeStepsInContWStatic() << endl;
    if(setStatic)
    {
        bodyOperation_ = 0;
        omega_ = 0;
        Vel_ *= 0;
        InfoH << basic_Info << "-- body " << bodyIdStr_ << " set as Static" << endl;
    }
}
//---------------------------------------------------------------------------//
void immersedBody::checkBodyOp()
{
    if(bodyOperation_ != 5 || timesToSetStatic_ == -1)
        return;

    if(!ibContactClass_->checkInContactWithStatic() && ibContactClass_->getTimeStepsInContWStatic() > 0)
    {
        ibContactClass_->setTimeStepsInContWStatic(0);
        return;
    }

    if(ibContactClass_->checkInContactWithStatic())
    {
        ibContactClass_->setTimeStepsInContWStatic(ibContactClass_->getTimeStepsInContWStatic() + 1);
        InfoH << iB_Info << "-- body " << bodyIdStr_
            << " timeStepsInContWStatic_: "
            << ibContactClass_->getTimeStepsInContWStatic() << endl;

        if(ibContactClass_->getTimeStepsInContWStatic() == 1)
        {
            staticContactPost_ = geomModel_->getCoM();
        }
        else
        {
            if(mag(staticContactPost_ - geomModel_->getCoM()) > 0.05 * geomModel_->getDC())
            {
                ibContactClass_->setTimeStepsInContWStatic(0);
                return;
            }
        }

        if(ibContactClass_->getTimeStepsInContWStatic() >= timesToSetStatic_)
        {
            bodyOperation_ = 0;
            omega_ = 0;
            Vel_ *= 0;
            InfoH << basic_Info << "-- body " << bodyIdStr_ << " set as Static" << endl;
        }
    }

    ibContactClass_->inContactWithStatic(false);
}

//---------------------------------------------------------------------------//
void immersedBody::updateRhoF
(
    const scalar rho
)
{    
    rhoF_ = rho;
}
void immersedBody::updateRhoF                                           //variant for VOF
(
    const volScalarField& rho,
    const volScalarField& body
)
{
    scalar fluidMass(0);
    scalar fluidVol(0);

    List<DynamicLabelList> relevantLists;
    geomModel_->getReferencedHaloCellList(relevantLists);
    // geomModel_->getReferencedInternalCellList(relevantLists);
    DynamicVectorList refCoMList;
    geomModel_->getReferencedCoMList(refCoMList);
    
    // Note (MI): in this case, we do not want to take into account the
    //            fluid composition inside the particle (frozen alpha field)
    // - we calculate the density of the surrounding fluid only from
    //   HALO cells
    // - weighting of the cell is done based on the fluid volume fraction
    
    // compute the weighted average of density        
    forAll (relevantLists, i)
    {
        DynamicLabelList& relevantListI = relevantLists[i];
        forAll (relevantListI, rCell)
        {
            label cellI = relevantListI[rCell];

            fluidMass += rho[cellI]*mesh_.V()[cellI]*(1.0 - body[cellI]);
            fluidVol  += mesh_.V()[cellI]*(1.0 - body[cellI]);
            // fluidMass += rho[cellI]*mesh_.V()[cellI];
            // fluidVol  += mesh_.V()[cellI];
        }
    }
    
    reduce(fluidMass, sumOp<scalar>());
    reduce(fluidVol, sumOp<scalar>());
    
    
    if (fluidVol > SMALL)
    {
        rhoF_ = fluidMass/fluidVol;
    }
    else
    {
        rhoF_ = 1.0;
    }
    InfoH << iB_Info << "-- body: " << bodyIdStr_ << ": rhoF = " << rhoF_ << endl;
}
