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
    abstract base class for fluid-solid coupling models. Computes the
    fluid-to-solid force on a body from a fluidContext. The coupling
    state (FCoupling_ and history) still lives in immersedBody for now;
    it will move here in a later step.

    Note (MI): forward (S->F) coupling is intentionally NOT part of this
               class yet. Extension point: virtual updateImposedField
               with a default no-op; unresolved models will override
               with point forces.

SourceFiles
    fluidSolidCoupling.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef fluidSolidCoupling_H
#define fluidSolidCoupling_H

#include "couplingContexts.H"

#include "fvMesh.H"
#include "geomModel.H"
#include "interpolationInfo.H"
#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class fluidSolidCoupling Declaration
\*---------------------------------------------------------------------------*/

class fluidSolidCoupling
{
    protected:

    // Protected data ----------------------------------------------

    const fvMesh& mesh_;

    // geomModel reference (aliased to the owner's shared_ptr, survives
    // geomModel replacement, e.g. the periodic-body swap)
    std::shared_ptr<geomModel>& geomModel_;

    // interpolation info of the owning body (ibPoints, findIbPoint);
    // null for bodies created without an interpolation scheme - only
    // the coupling physics needs it, the state below does not
    interpolationInfo* intpInfo_;

    // bodyId string for output
    const word& bodyIdStr_;

    // coupling state (moved from immersedBody, step C)

    forces FCoupling_;
    forces FCouplingOld_;
    scalar couplingHistCoef_;
    dimensionedScalar rhoF_;

    public:

    // Constructors -----------------------------------------------

    fluidSolidCoupling
    (
        const fvMesh& mesh,
        std::shared_ptr<geomModel>& geomModel,
        interpolationInfo* intpInfo,
        const word& bodyIdStr
    );

    // Destructors ------------------------------------------------

    virtual ~fluidSolidCoupling();

    // Member functions -------------------------------------------

    //- compute the coupling force on the body from the fluid context
    virtual void updateCoupling
    (
        const fluidContext& ctx,
        const vector& a
    ) = 0;

    //- coupling force history accessors (movement / print consumers)
    const forces& getFCoupling() const
    {
        return FCoupling_;
    }

    const forces& getFCouplingOld() const
    {
        return FCouplingOld_;
    }

    //- fluid density seen by the body (buoyancy consumers)
    const dimensionedScalar& getRhoF() const
    {
        return rhoF_;
    }

    //- reset the under-relaxation history after a contact event
    void resetCouplingHistory()
    {
        couplingHistCoef_ = 1.0;
    }

    //- carry the current coupling force into the old-iteration slot
    //- (per PIMPLE outer iteration)
    void resetPostPimpleState()
    {
        FCouplingOld_ = FCoupling_;
    }

    //- update the fluid density seen by the body (single-phase
    //- solvers: constant fluid density)
    void updateRhoF(const scalar rho)
    {
        rhoF_ = rho;
    }

    //- update the fluid density seen by the body (VOF solvers:
    //- halo-cell fluid-volume-weighted average)
    void updateRhoF(const volScalarField& rho, const volScalarField& body);

    //- write the coupling state into the restart record
    //- (entries are optional for readers: lookupOrDefault with
    //- constructor defaults keeps old restart files loadable)
    void writeCouplingInfo(dictionary& dict) const;

    //- read the coupling state back from a restart record
    void readCouplingInfo(const dictionary& dict);

    //- select and construct a coupling model from the body dictionary
    static autoPtr<fluidSolidCoupling> New
    (
        const dictionary& bodyDict,
        const fvMesh& mesh,
        std::shared_ptr<geomModel>& geomModel,
        interpolationInfo* intpInfo,
        const word& bodyIdStr
    );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
// ************************************************************************* //
