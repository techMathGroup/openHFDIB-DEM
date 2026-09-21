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
    contexts passed between the fluid solver and the fluid-solid coupling
    models. One combined force field is always provided - the solvers
    resolve the pressure and viscous contributions themselves.

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef couplingContexts_H
#define couplingContexts_H

#include "volFieldsFwd.H"
#include "dimensionedScalar.H"

#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- which fluid situation the calling solver is in
enum class fluidSituation
{
    singlePhase,                                                        //one fluid
    voF                                                                   //VOF
};

//- data bundle describing the fluid as seen by the coupling model
struct fluidContext
{
    const volScalarField* body;                                         //body field
    const volVectorField* f;                                        //combined force
    const volScalarField* rho;                                        //fluid density
    fluidSituation situation;                                         //solver type
    bool kinematicForce;                                    //multiply output by rhoF
    bool applyAddedMass;                                     //added-mass scaling
};

//- coupling state owned by immersedBody (until the state move). The
//- model computes into it; pointers are refreshed by the caller per call
struct solidContext
{
    forces* FCoupling;                                             //output force
    const forces* FCouplingOld;                             //previous iteration
    scalar* couplingHistCoef;                                 //history coefficient
    dimensionedScalar* rhoF;                                       //fluid density
    const vector* a;                                            //body acceleration
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
// ************************************************************************* //
