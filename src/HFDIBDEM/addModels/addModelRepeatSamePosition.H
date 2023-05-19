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
    addModelRepeatSamePosition.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef addModelRepeatSamePosition_H
#define addModelRepeatSamePosition_H

#include "addModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class addModelRepeatSamePosition Declaration
\*---------------------------------------------------------------------------*/

class addModelRepeatSamePosition
:
    public addModel
{

    private:
    const dictionary&       addModelDict_;
    word                    addMode_;
    bool                    bodyAdded_;

    const dictionary&       coeffsDict_;
    label                   useNTimes_;
    scalar                  timeBetweenUsage_;
    label                   addedOnTimeLevel_;

    public:

    // Constructors ------------------------------------------------
    addModelRepeatSamePosition
    (
        const dictionary& addModelDict,
        const Foam::fvMesh& mesh,
        std::unique_ptr<geomModel> bodyGeomModel,
        List<labelList>& cellPoints
    );

    // Destructors -------------------------------------------------
    ~addModelRepeatSamePosition();

    // Public methods ----------------------------------------------
    //~ void init();

    bool shouldAddBody(const volScalarField& body);

    std::shared_ptr<geomModel> addBody(
        const volScalarField& body,
        PtrList<immersedBody>& immersedBodies
    );

    bool checkIfAdd()
    {
        return true;
    }

    void translateRandom()
    {
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "Not needed"<< endl;
    }

    void rotateRandom()
    {
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "Not needed"<< endl;
    }

    word getAddMode()
    {
        return addMode_;
    }

    bool getBodyAdded(){return bodyAdded_;}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
