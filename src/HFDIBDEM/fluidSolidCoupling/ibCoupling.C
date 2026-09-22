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

Description
    immersed-boundary coupling model

Note (MI, 20260922):
    - currently, the implementations for VOF and single-phase
      are slightly different but mostly doubled
    - the difference in implementation comes from uncertainty about
      physics modeling
    - also, the updateCouplingTail() should be made more generic and
      shared between the two implementations

SourceFiles
    ibCoupling.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#include "ibCoupling.H"

#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Constructors ----------------------------------------------------

ibCoupling::ibCoupling
(
    const fvMesh& mesh,
    std::shared_ptr<geomModel>& geomModel,
    interpolationInfo& intpInfo,
    const word& bodyIdStr
)
:
fluidSolidCoupling(mesh, geomModel, intpInfo, bodyIdStr)
{}

// Destructors -----------------------------------------------------

ibCoupling::~ibCoupling()
{}

// Member functions ------------------------------------------------

void ibCoupling::updateCoupling
(
    const fluidContext& ctx,
    solidContext& sCtx
)
{
    if (!ctx.f)
    {
        FatalErrorInFunction
            << "ibCoupling: no force field in the fluid context"
            << abort(FatalError);
    }

    const volVectorField& f = *ctx.f;

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

    const List<point>& ibPoints = intpInfo_.getIbPoints();              //get surface points

    if (ctx.situation == fluidSituation::singlePhase)
    {
        updateCouplingSinglePhase
        (
            ctx,
            f,
            intLists,
            surfLists,
            refCoMList,
            ibPoints,
            FV,
            TA,
            FAdded
        );
    }
    else
    {
        updateCouplingVOF
        (
            ctx,
            sCtx,
            f,
            intLists,
            surfLists,
            haloLists,
            refCoMList,
            ibPoints,
            FV,
            TA,
            FAdded
        );
    }

    updateCouplingTail(ctx, sCtx, FV, TA, FAdded);
}

//---------------------------------------------------------------------------//
void ibCoupling::updateCouplingSinglePhase
(
    const fluidContext& ctx,
    const volVectorField& f,
    List<DynamicLabelList>& intLists,
    List<DynamicLabelList>& surfLists,
    DynamicVectorList& refCoMList,
    const List<point>& ibPoints,
    vector& FV,
    vector& TA,
    vector& FAdded
)
{
    // Note (MI): current idea is that action-reaction force
    //            between solid and fluid happens only on the body
    //            surface, what happens inside is solid-to-solid??
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

    forAll (surfLists, i)
    {
        DynamicLabelList& surfListI = surfLists[i];
        forAll (surfListI, surfCell)
        {
            label cellI = surfListI[surfCell];

            // Note (MI): current idea is that there is a discrepancy
            //            between the force required to fullfill boundary
            //            conditions and the force actually acting on
            //            the fluid
            // scalar fScale = 1.0*body[cellI]+0.5;
            // vector fCell = (1.0 - body[cellI])*f[cellI];
            scalar fScale = 1.0;
            vector fCell =  f[cellI]*mesh_.V()[cellI];
            fCell *= fScale;

            const vector& surfPoint = ibPoints[intpInfo_.findIbPoint(cellI)];

            FV -=  fCell;
            TA -=  (surfPoint - refCoMList[i])^fCell;
            // Note (MI): f.prevIter() is only valid when the solver stored
            //            it; accumulate added-mass force only when needed
            if (ctx.applyAddedMass)
            {
                FAdded -= (f.prevIter()[cellI] - f[cellI])
                    *mesh_.V()[cellI];//under construction
            }
        }
    }
}

//---------------------------------------------------------------------------//
void ibCoupling::updateCouplingVOF
(
    const fluidContext& ctx,
    solidContext& sCtx,
    const volVectorField& f,
    List<DynamicLabelList>& intLists,
    List<DynamicLabelList>& surfLists,
    List<DynamicLabelList>& haloLists,
    DynamicVectorList& refCoMList,
    const List<point>& ibPoints,
    vector& FV,
    vector& TA,
    vector& FAdded
)
{
    if (!ctx.rho || !ctx.body)
    {
        FatalErrorInFunction
            << "ibCoupling: no density or body field in the fluid context"
            << abort(FatalError);
    }

    const volScalarField& rho = *ctx.rho;
    const volScalarField& body = *ctx.body;
    const dimensionedScalar& rhoF = *sCtx.rhoF;

    List<List<DynamicLabelList>>& surfToHaloAddressing
        = geomModel_->getSurfToHaloLabels();

    forAll (intLists, i)
    {
        DynamicLabelList& intListI = intLists[i];
        forAll (intListI, intCell)
        {
            label cellI = intListI[intCell];

            scalar fScale = rhoF.value()/rho[cellI];
            vector fCell =  f[cellI]*mesh_.V()[cellI];
            fCell *= fScale;

            FV -=  fCell;
            TA -=  ((mesh_.C()[cellI] - refCoMList[i])^fCell);
            if (ctx.applyAddedMass)
            {
                FAdded -= (f.prevIter()[cellI] - f[cellI])
                    *mesh_.V()[cellI];
            }
        }
    }

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

            FV -=  fCell;
            TA -=  (ibPoints[intpInfo_.findIbPoint(cellI)] - refCoMList[i])^fCell;
            if (ctx.applyAddedMass)
            {
                FAdded -= (f.prevIter()[cellI] - f[cellI])
                    *mesh_.V()[cellI];//under construction
            }
        }
    }
}

//---------------------------------------------------------------------------//
void ibCoupling::updateCouplingTail
(
    const fluidContext& ctx,
    solidContext& sCtx,
    const vector& FV,
    const vector& TA,
    const vector& FAdded
)
{
    forces& FCoupling = *sCtx.FCoupling;
    const forces& FCouplingOld = *sCtx.FCouplingOld;
    scalar& couplingHistCoef = *sCtx.couplingHistCoef;
    dimensionedScalar& rhoF = *sCtx.rhoF;
    const vector& a = *sCtx.a;

    vector FVl(FV);
    vector TAl(TA);
    vector FAddedl(FAdded);

    if (ctx.kinematicForce)
    {
        FVl *= rhoF.value();
        TAl *= rhoF.value();
        FAddedl *= rhoF.value();
    }

    scalar rhoS = geomModel_->getRhoS().value();
    // FVl /= rhoS;
    // TAl /= rhoS;
    // FAddedl /= rhoS;
    // FVl *= rhoF.value();
    // TAl *= rhoF.value();
    // FAddedl *= rhoF.value();

    // FAddedl = FCouplingOld.F - FVl;

    FAddedl *= rhoF.value()/rhoS;

    FCoupling = couplingHistCoef*forces(FVl, TAl) + (1.0-couplingHistCoef)*FCouplingOld;

    if (ctx.applyAddedMass)
    {
        const scalar m0 = geomModel_->getM0();
        const scalar massSign = ((FVl & FAddedl) < 0.0) ? -1.0 : 1.0;
        // scalar massAdded = min(1.0*m0, mag(FAddedl)/(mag(a_) + SMALL));
        scalar massAdded = mag(FAddedl)/(mag(a) + SMALL);
        massAdded *= massSign;
        InfoH << iB_Info << "-- body " << bodyIdStr_ << " massAdded: " << massAdded
            << " m0: " << m0 << endl;
        InfoH << iB_Info << "-- body " << bodyIdStr_ << " orig coupling force: " << FCoupling.F << " orig coupling torque: " << FCoupling.T << endl;
        const scalar scale = (m0 + massAdded)/m0;
        FCoupling.F *= scale;
        FCoupling.T *= scale;
        InfoH << iB_Info << "-- body " << bodyIdStr_ << " scld coupling force: " << FCoupling.F << " scld coupling torque: " << FCoupling.T << endl;
    }

    couplingHistCoef = max(couplingHistCoef*0.95, 0.5);

    InfoH << iB_Info << "-- body " << bodyIdStr_ << " coupling coefficient: " << couplingHistCoef << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
