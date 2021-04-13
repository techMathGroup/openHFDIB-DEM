/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________
                       | | | ||  ___|  _  \_   _| ___ \     H ybrid
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /     F ictitious
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \     D omain
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /     I mmersed
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/      B oundary
      | |
      |_|
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
    Martin Isoz (2019-*), Martin Šourek (2019-*), 
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "addModelRepeatSamePosition.H"

using namespace Foam;

//---------------------------------------------------------------------------//
addModelRepeatSamePosition::addModelRepeatSamePosition
(
    const dictionary& addModelDict,
    const word        stlName,
    const Foam::dynamicFvMesh& mesh
)
:
addModelDict_(addModelDict),
addMode_(word(addModelDict_.lookup("addModel"))),
stlName_(stlName),
bodyAdded_(false),
mesh_(mesh),
coeffsDict_(addModelDict_.subDict(addMode_+"Coeffs")),
useNTimes_(readLabel(coeffsDict_.lookup("useNTimes"))),
timeBetweenUsage_(readScalar(coeffsDict_.lookup("timeBetweenUsage"))),
addedOnTimeLevel_(0)
{
}
    
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
    
    Info << "-- addModelMessage-- " << "Time/(Time beween usage) - floor(Time/Time beween usage): " 
         << tmFrac << endl;
         
    Info << "-- addModelMessage-- " << "Number of bodies added on this time level: " << addedOnTimeLevel_ << endl;
         
    bool tmLevelOk(tmFrac < deltaTime);
    
    if (not tmLevelOk){addedOnTimeLevel_ = 0;}
             
    return (tmLevelOk and useNTimes_ > 0 and addedOnTimeLevel_ == 0);
}

triSurface addModelRepeatSamePosition::addBody
(
    const   volScalarField& body
)
{
    triSurfaceMesh bodySurfMesh
    (
        IOobject
        (
            stlName_ +".stl",
            "constant",
            "triSurface",
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    bool canAddBodyI(canAddBody(body,bodySurfMesh));
    reduce(canAddBodyI, andOp<bool>());
    
    bodyAdded_ = canAddBodyI;
    if (bodyAdded_) {useNTimes_--;}
    addedOnTimeLevel_++;
    
    Info << "-- addModelMessage-- " << "will try to use the body " << useNTimes_ << " more times" << endl;
    
    triSurface triToRet(bodySurfMesh);
    
    return triToRet;
}

bool addModelRepeatSamePosition::canAddBody
(
    const volScalarField& body,
    const triSurfaceMesh& bodySurfMesh
)
{
    #include "canAddBodySource.H"
}
