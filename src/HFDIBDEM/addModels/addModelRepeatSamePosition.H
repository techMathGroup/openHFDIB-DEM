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

    void recreateBoundBox();

    bool getBodyAdded(){return bodyAdded_;}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
