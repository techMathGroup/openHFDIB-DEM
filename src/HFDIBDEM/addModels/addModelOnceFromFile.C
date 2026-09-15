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
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "addModelOnceFromFile.H"

#include <memory>



#include "IFstream.H"
#include "OSspecific.H"


using namespace Foam;

//---------------------------------------------------------------------------//
addModelOnceFromFile::addModelOnceFromFile
(
    const dictionary& addModelDict,
    const Foam::fvMesh& mesh,
    const bool startTime0,
    std::unique_ptr<geomModel> bodyGeomModel,
    List<labelList>& cellPoints,
    word& bodyGeom,
    scalar thrSurf
)
:
addModel(mesh, std::move(bodyGeomModel), cellPoints),
addModelDict_(addModelDict),
addMode_(word(addModelDict_.lookup("addModel"))),
coeffsDict_(addModelDict_.subDict(addMode_+"Coeffs")),
bodyAdded_(false),
fileName_("constant/" + (word(coeffsDict_.lookup("fileName")))),
ifStream_(fileName_.toAbsolute()),
bodyGeom_(bodyGeom),
thrSurf_(thrSurf)
{
    if(!ifStream_.opened())
    {
        FatalErrorIn("addModelOnceFromFile::addModelOnceFromFile()")
        << "Cannot open IFstream for file: "
        << fileName_.toAbsolute() << nl << exit(FatalError);
    }

    if(!startTime0)
    {
        ifStream_.setEof();
    }
}

addModelOnceFromFile::~addModelOnceFromFile()
{
}
//---------------------------------------------------------------------------//
void addModelOnceFromFile::addSphere(string& line)
{
    IStringStream stringStream(line);
    vector position(stringStream);

    geomModel_->bodyMovePoints(position - geomModel_->getCoM());
}
//---------------------------------------------------------------------------//
void addModelOnceFromFile::addSTL(string& line)
{
    if(bodyGeom_ == "convex")
    {
        word stlPath("constant/triSurface/" + line);
        geomModel_ = std::unique_ptr<convexBody>
            (new convexBody(mesh_, stlPath, thrSurf_));
    }
    else
    {
        word stlPath("constant/triSurface/" + line);
        geomModel_ = std::unique_ptr<nonConvexBody>
            (new nonConvexBody(mesh_, stlPath, thrSurf_));
    }
}
//---------------------------------------------------------------------------//
std::shared_ptr<geomModel> addModelOnceFromFile::addBody
(
    const volScalarField& body,
    PtrList<immersedBody>& immersedBodies
)
{
    string line;
    ifStream_.getLine(line);
    if(line == "")
    {
        InfoH << addModel_Info << "-- addModelMessage-- "
              << "Skipping empty line at " << ifStream_.lineNumber() - 1 <<endl;
        bodyAdded_ = false;
        return geomModel_->getCopy();
    }

    if (bodyGeom_ == "sphere")
    {
        addSphere(line);
    }
    else
    {
        addSTL(line);
    }

    bodyAdded_ = true;
    return geomModel_->getCopy();
}
