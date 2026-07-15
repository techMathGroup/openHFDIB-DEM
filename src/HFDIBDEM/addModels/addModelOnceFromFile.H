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

Description
    class for addition and removal of immersedBodies

SourceFiles
    addModelOnceFromFile.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef addModelOnceFromFile_H
#define addModelOnceFromFile_H

#include "addModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class addModelOnce Declaration
\*---------------------------------------------------------------------------*/

class addModelOnceFromFile
:
    public addModel
{

    private:
    const dictionary&       addModelDict_;
    word                    addMode_;
    const dictionary&       coeffsDict_;
    bool                    bodyAdded_;

    fileName                fileName_;

    //there might be a problem in OF2406 and similar
    IFstream                ifStream_;


    word                    bodyGeom_;
    scalar                  thrSurf_;

    // Private Member Functions ------------------------------------

    void addSphere(string& line);

    void addSTL(string& line);

    public:

    // Constructors ------------------------------------------------
    addModelOnceFromFile
    (
        const dictionary& addModelDict,
        const Foam::fvMesh& mesh,
        const bool startTime0,
        std::unique_ptr<geomModel> bodyGeomModel,
        List<labelList>& cellPoints,
        word& bodyGeom,
        scalar thrSurf
    );

    // Destructors -------------------------------------------------
    ~addModelOnceFromFile();

    // Public methods ----------------------------------------------

    bool shouldAddBody(const volScalarField& body)
    {
        return not ifStream_.eof();
    }

    std::shared_ptr<geomModel> addBody(
        const volScalarField& body,
        PtrList<immersedBody>& immersedBodies
        );

    bool getBodyAdded()
    {
        return bodyAdded_;
    }
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
