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
    contactInfo.C

Contributors
    Martin Isoz (2019-*), Martin Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef contactInfo_H
#define contactInfo_H

#include "dictionary.H"
#include "fvCFD.H"

#include "geomModel.H"
#include "contactVars.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class contactInfo Declaration
\*---------------------------------------------------------------------------*/

class contactInfo
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

    protected:

    geomModel& geomModel_;

    // is the particle in contact with wall
    bool    isInWallContact_;

    List<DynamicLabelList>& surfCells_;
    List<DynamicLabelList>& intCells_;

    bool inContactWithStatic_;

    label timeStepsInContWStatic_;

    // faces in contact with walls
    List<DynamicLabelList>  wallContactFaces_;

    // spring stiffness for normal contact force
    scalar  kN_;

    // spring stiffness for tangential contact force
    scalar  kt_;

    // damping coefficient for normal contact force
    scalar gammaN_;

    // damping coefficient for tangential contact force
    scalar gammat_;

    // friction coefficient
    scalar mu_;

    // adhesion coefficient for normal contact force
    scalar adhN_;

    // adhesion coefficient for normal contact force
    scalar adhEqui_;

    // tangential force from last timestep
    // first label -1 indicates wall other for bodyID
    DynamicList<Tuple2<label,Tuple2<label,vector>>> historyFt_;

    public:

    // Constructors ------------------------------------------------
    contactInfo
    (
        geomModel& geomModel,
        List<DynamicLabelList>& surfCells,
        List<DynamicLabelList>& intCells,
        scalar kN,
        scalar kt,
        scalar gammaN,
        scalar gammat,
        scalar mu,
        scalar adhN,
        scalar adhEqui
    );

    // Destructors -------------------------------------------------
    ~contactInfo();

    geomModel& getGeomModel()
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
        return surfCells_;
    }

    List<DynamicLabelList>& getIntCells()
    {
        return intCells_;
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

    List<DynamicLabelList>& wallContactFaces()
    {
        return wallContactFaces_;
    }

    void distributeWallContactFaces()
    {
        // distribute wallContactfaces over processors
        Pstream::gatherList(wallContactFaces_, 0);
        Pstream::scatter(wallContactFaces_, 0);
    }

    const scalar& getkN()   const
    {
        return kN_;
    }

    const scalar& getkt()   const
    {
        return kt_;
    }

    const scalar& getgammaN()   const
    {
        return gammaN_;
    }

    const scalar& getgammat()   const
    {
        return gammat_;
    }

    const scalar& getmu()   const
    {
        return mu_;
    }

    const scalar& getAdhN()   const
    {
        return adhN_;
    }

    const scalar& getAdhEqui()   const
    {
        return adhEqui_;
    }

    DynamicList<Tuple2<label,Tuple2<label,vector>>>& getHistoryhistoryFt()
    {
        return historyFt_;
    }

    void clearContactInfo()
    {
        isInWallContact_ = false;
        wallContactFaces_.clear();
        wallContactFaces_.setSize(Pstream::nProcs());
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
