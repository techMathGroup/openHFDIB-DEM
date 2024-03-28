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

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
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

//~ void init()
//~ {
//~ }

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
