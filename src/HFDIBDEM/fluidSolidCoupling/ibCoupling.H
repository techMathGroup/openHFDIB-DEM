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
    immersed-boundary coupling model. Fluid-to-solid force assembled from
    the surface-cell force field. State (FCoupling_ and history) stays in
    immersedBody for now - this class computes into it through the
    references handed to updateCoupling.

    The computation switches on ctx.situation FIRST, then on pointer
    nullness. Only the singlePhase situation is implemented; voF follows
    in a later step.

SourceFiles
    ibCoupling.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef ibCoupling_H
#define ibCoupling_H

#include "fluidSolidCoupling.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ibCoupling Declaration
\*---------------------------------------------------------------------------*/

class ibCoupling
:
    public fluidSolidCoupling
{
    private:

    // Private data types -----------------------------------------

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

    // Private member functions -----------------------------------

    // Note (MI): the assemblies should be identical/similar for both
    //            situations, but the physics is not clear yet.
    //- single-phase force assembly (surf cells only)
    void updateCouplingSinglePhase
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
    );

    //- VOF force assembly (int + surf cells, halo density estimate)
    void updateCouplingVOF
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
    );

    //- shared post-processing: scaling, under-relaxation, added mass
    void updateCouplingTail
    (
        const fluidContext& ctx,
        solidContext& sCtx,
        const vector& FV,
        const vector& TA,
        const vector& FAdded
    );

    public:

    // Constructors -----------------------------------------------

    ibCoupling
    (
        const fvMesh& mesh,
        std::shared_ptr<geomModel>& geomModel,
        interpolationInfo& intpInfo,
        const word& bodyIdStr
    );

    // Destructors ------------------------------------------------

    virtual ~ibCoupling();

    // Member functions -------------------------------------------

    void updateCoupling
    (
        const fluidContext& ctx,
        solidContext& sCtx
    );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
// ************************************************************************* //
