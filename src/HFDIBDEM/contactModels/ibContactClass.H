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
    class for storage of contact variables

SourceFiles
    ibContactClass.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef ibContactClass_H
#define ibContactClass_H

#include "dictionary.H"
#include "fvCFD.H"

#include "geomModel.H"
#include "ibContactVars.H"
#include "wallContactVars.H"
#include "materialInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ibContactClass Declaration
\*---------------------------------------------------------------------------*/

class ibContactClass
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

    protected:

    std::shared_ptr<geomModel>& geomModel_;

    // is the particle in contact with wall
    bool    isInWallContact_;

    bool inContactWithStatic_;

    label timeStepsInContWStatic_;

    const materialInfo& matInfo_;

    public:

    // Constructors ------------------------------------------------
    ibContactClass
    (
        std::shared_ptr<geomModel>& geomModel,
        const string& material
    );

    ibContactClass
    (
        const ibContactClass& other
    );

    // Destructors -------------------------------------------------
    ~ibContactClass();

    geomModel& getGeomModel()
    {
        return *geomModel_;
    }

    std::shared_ptr<geomModel> getGeomModelPtr()
    {
        return geomModel_;
    }

    const bool& checkWallContact() const
    {
        return isInWallContact_;
    }

    void setWallContact(bool inContact)
    {
        isInWallContact_ = inContact;
    }

    List<DynamicLabelList>& getSurfCells()
    {
        return geomModel_->getSurfaceCellList();
    }

    List<DynamicLabelList>& getIntCells()
    {
        return geomModel_->getInternalCellList();
    }

    const materialInfo& getMatInfo()
    {
        return matInfo_;
    }

    const bool& checkInContactWithStatic() const
    {
        return inContactWithStatic_;
    }

    void inContactWithStatic(bool inContact)
    {
        inContactWithStatic_ = inContact;
    }

    const label& getTimeStepsInContWStatic() const
    {
        return timeStepsInContWStatic_;
    }

    void setTimeStepsInContWStatic(label setter)
    {
        timeStepsInContWStatic_ = setter;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
