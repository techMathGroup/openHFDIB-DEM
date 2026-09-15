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

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/
#include "addModelRepeatSamePosition.H"

using namespace Foam;

//---------------------------------------------------------------------------//
addModelRepeatSamePosition::addModelRepeatSamePosition
(
    const dictionary& addModelDict,
    const Foam::fvMesh& mesh,
    std::unique_ptr<geomModel> bodyGeomModel,
    List<labelList>& cellPoints
)
:
addModel(mesh, std::move(bodyGeomModel), cellPoints),
addModelDict_(addModelDict),
addMode_(word(addModelDict_.lookup("addModel"))),
bodyAdded_(false),
coeffsDict_(addModelDict_.subDict(addMode_+"Coeffs")),
useNTimes_(readLabel(coeffsDict_.lookup("useNTimes"))),
timeBetweenUsage_(readScalar(coeffsDict_.lookup("timeBetweenUsage"))),
addedOnTimeLevel_(0)
{}

addModelRepeatSamePosition::~addModelRepeatSamePosition()
{
}


//---------------------------------------------------------------------------//
bool addModelRepeatSamePosition::shouldAddBody(const volScalarField& body)
{
    scalar timeVal(mesh_.time().value());
    scalar deltaTime(mesh_.time().deltaT().value());
    scalar tmFrac(timeVal/timeBetweenUsage_);
    tmFrac -=  floor(tmFrac+deltaTime);

    InfoH << addModel_Info << "-- addModelMessage-- "
        << "Time/(Time beween usage) - floor(Time/Time beween usage): "
        << tmFrac << endl;

    InfoH << "-- addModelMessage-- "
        << "Number of bodies added on this time level: "
        << addedOnTimeLevel_ << endl;

    bool tmLevelOk(tmFrac < deltaTime);

    if (not tmLevelOk){addedOnTimeLevel_ = 0;}

    return (tmLevelOk and useNTimes_ > 0 and addedOnTimeLevel_ == 0);
}

std::shared_ptr<geomModel> addModelRepeatSamePosition::addBody
(
    const volScalarField& body,
    PtrList<immersedBody>& immersedBodies
)
{
    volScalarField helpBodyField_ = body;
    geomModel_->createImmersedBody(
        helpBodyField_,
        octreeField_,
        cellPoints_
    );

    bool canAddBodyI = !isBodyInContact(immersedBodies);

    reduce(canAddBodyI, andOp<bool>());

    bodyAdded_ = canAddBodyI;
    if (bodyAdded_) {useNTimes_--;}
    addedOnTimeLevel_++;

    InfoH << addModel_Info << "-- addModelMessage-- "
        << "will try to use the body " << useNTimes_ << " more times" << endl;

    return geomModel_->getCopy();
}
//---------------------------------------------------------------------------//
void addModelRepeatSamePosition::recreateBoundBox()
{
    // no random bounding box to recompute in this model, but the octreeField_
    // must be re-sized after a mesh change (e.g. refinement) so that the
    // candidate-cell search in addBody does not write out of bounds
    octreeField_ = Field<label>(mesh_.nCells(), 0);
}
