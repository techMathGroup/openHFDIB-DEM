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
    abstract base class for fluid-solid coupling models

SourceFiles
    fluidSolidCoupling.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#include "fluidSolidCoupling.H"
#include "ibCoupling.H"

#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Constructors ----------------------------------------------------

fluidSolidCoupling::fluidSolidCoupling
(
    const fvMesh& mesh,
    std::shared_ptr<geomModel>& geomModel,
    interpolationInfo* intpInfo,
    const word& bodyIdStr
)
:
mesh_(mesh),
geomModel_(geomModel),
intpInfo_(intpInfo),
bodyIdStr_(bodyIdStr),
FCoupling_(),
FCouplingOld_(),
couplingHistCoef_(1.0),
rhoF_(1.0)
{}

// Destructors -----------------------------------------------------

fluidSolidCoupling::~fluidSolidCoupling()
{}

// Member functions ------------------------------------------------

void fluidSolidCoupling::updateRhoF
(
    const volScalarField& rho,
    const volScalarField& body
)
{
    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<vector> DynamicVectorList;

    scalar fluidMass(0);
    scalar fluidVol(0);

    List<DynamicLabelList> relevantLists;
    geomModel_->getReferencedHaloCellList(relevantLists);
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

autoPtr<fluidSolidCoupling> fluidSolidCoupling::New
(
    const dictionary& bodyDict,
    const fvMesh& mesh,
    std::shared_ptr<geomModel>& geomModel,
    interpolationInfo* intpInfo,
    const word& bodyIdStr
)
{
    const word modelType
    (
        bodyDict.lookupOrDefault<word>("couplingModel", "IB")
    );

    if (modelType == "IB")
    {
        return autoPtr<fluidSolidCoupling>
        (
            new ibCoupling(mesh, geomModel, intpInfo, bodyIdStr)
        );
    }

    FatalErrorInFunction
        << "Unknown fluidSolidCoupling model type: " << modelType
        << ". Valid model types: IB"
        << abort(FatalError);

    return autoPtr<fluidSolidCoupling>(nullptr);
}

//---------------------------------------------------------------------------//
void fluidSolidCoupling::writeCouplingInfo(dictionary& dict) const
{
    dict.add("FCouplingF", FCoupling_.F);
    dict.add("FCouplingT", FCoupling_.T);
    dict.add("FCouplingOldF", FCouplingOld_.F);
    dict.add("FCouplingOldT", FCouplingOld_.T);
    dict.add("couplingHistCoef", couplingHistCoef_);
    dict.add("rhoF", rhoF_.value());
}

//---------------------------------------------------------------------------//
void fluidSolidCoupling::readCouplingInfo(const dictionary& dict)
{
    FCoupling_.F = dict.lookupOrDefault<vector>("FCouplingF", vector::zero);
    FCoupling_.T = dict.lookupOrDefault<vector>("FCouplingT", vector::zero);
    FCouplingOld_.F = dict.lookupOrDefault<vector>
    (
        "FCouplingOldF",
        vector::zero
    );
    FCouplingOld_.T = dict.lookupOrDefault<vector>
    (
        "FCouplingOldT",
        vector::zero
    );
    couplingHistCoef_ = dict.lookupOrDefault<scalar>
    (
        "couplingHistCoef",
        1.0
    );
    rhoF_ = dict.lookupOrDefault<scalar>("rhoF", 1.0);

    InfoH << iB_Info << "-- body " << bodyIdStr_
        << " restart coupling state: FCoupling = (" << FCoupling_.F << " "
        << FCoupling_.T << "), couplingHistCoef = " << couplingHistCoef_
        << ", rhoF = " << rhoF_.value() << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
