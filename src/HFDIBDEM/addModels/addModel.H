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
    addModel.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef addModel_H
#define addModel_H

#include "dictionary.H"
#include "fvCFD.H"
// #include "fvMesh.H"
#include "triSurfaceMesh.H"

#include "geomModel.H"
#include "convexBody.H"
#include "nonConvexBody.H"
#include "sphereBody.H"
#include "periodicBody.H"

#include "contactModel.H"
#include "immersedBody.H"
#include "ibContactVars.H"

#include "virtualMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// class immersedBody;

/*---------------------------------------------------------------------------*\
                          Class addModel Declaration
\*---------------------------------------------------------------------------*/

class addModel
{

    private:

    protected:

    const  fvMesh&       mesh_;
    List<labelList>& cellPoints_;
    Field<label> octreeField_;
    std::unique_ptr<geomModel> geomModel_;

    public:

    // Constructors ------------------------------------------------
    addModel(
        const fvMesh& mesh,
        std::unique_ptr<geomModel> bodyGeomModel,
        List<labelList>& cellPoints
    )
    :
    mesh_(mesh),
    cellPoints_(cellPoints),
    geomModel_(std::move(bodyGeomModel))
    {
        octreeField_ = Field<label>(mesh_.nCells(), 0);
    }

    // Destructors -------------------------------------------------
    virtual ~addModel();

    // Virtual methods to be overriden------------------------------
    // Note:
    //  - all these needs to have the same arugments accross
    //    all the addModels (we need to keep a common interface)

    // -- init
    virtual void init()
    {
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "notIplemented"<<endl;
    }
    // Note: rewrite this function if you need initialization
    //       and place it into constructor

    virtual bool shouldAddBody(const volScalarField& body)
    {
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "notIplemented"<<endl;
        return false;
    }
    // Note: check current conditions in simulation and return
    //       true if a body should be added

    virtual std::shared_ptr<geomModel> addBody(
        const volScalarField& body,
        PtrList<immersedBody>& immersedBodies
        ) = 0;
    // Note: needs to provide a properly scaled/oriented/translated
    //       triSurface that will be used for IB initialization

    // -- getters
    virtual word getAddMode()
    {
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "notIplemented"<<endl;
        return "NO";
    }

    virtual bool getBodyAdded()
    {
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "notIplemented"<<endl;
        return false;
    }
    // Note: this has to return true if addBody ends with SUCCESS,
    //       return false otherwise

    virtual void recreateBoundBox()
    {
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "notIplemented"<<endl;
    }

    bool isBodyInContact(PtrList<immersedBody>& immersedBodies);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
