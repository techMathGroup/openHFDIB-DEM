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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Constructors ----------------------------------------------------

fluidSolidCoupling::fluidSolidCoupling
(
    const fvMesh& mesh,
    std::shared_ptr<geomModel>& geomModel,
    interpolationInfo& intpInfo,
    const word& bodyIdStr
)
:
mesh_(mesh),
geomModel_(geomModel),
intpInfo_(intpInfo),
bodyIdStr_(bodyIdStr)
{}

// Destructors -----------------------------------------------------

fluidSolidCoupling::~fluidSolidCoupling()
{}

// Member functions ------------------------------------------------

autoPtr<fluidSolidCoupling> fluidSolidCoupling::New
(
    const dictionary& bodyDict,
    const fvMesh& mesh,
    std::shared_ptr<geomModel>& geomModel,
    interpolationInfo& intpInfo,
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
