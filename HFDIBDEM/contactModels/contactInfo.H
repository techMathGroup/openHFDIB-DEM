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
#include "wallContactVars.H"

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

    List<DynamicLabelList> surfCells_;
    List<DynamicLabelList> intCells_;

    bool inContactWithStatic_;

    label timeStepsInContWStatic_;

    // faces in contact with walls
    List<wallContactVars> wallContactsVars;

    // Young modulus
    scalar  Y_;

    // poisson ration
    scalar  nu_;

    // dissipative coef
    scalar gamma_;

    // friction coefficient
    scalar mu_;

    // adhesion coefficient for normal contact force
    scalar adhN_;

    // tangential force from last timestep
    // first label -1 indicates wall other for bodyID
    DynamicList<Tuple2<label,Tuple2<label,vector>>> historyFt_;

    public:

    // Constructors ------------------------------------------------
    contactInfo
    (
        geomModel& geomModel,
        scalar Y,
        scalar nu,
        scalar gamma,
        scalar mu,
        scalar adhN
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

    const scalar& getY()   const
    {
        return Y_;
    }

    const scalar& getNu()   const
    {
        return nu_;
    }

    const scalar& getGamma()   const
    {
        return gamma_;
    }

    const scalar& getmu()   const
    {
        return mu_;
    }

    const scalar& getAdhN()   const
    {
        return adhN_;
    }

    DynamicList<Tuple2<label,Tuple2<label,vector>>>& getHistoryFt()
    {
        return historyFt_;
    }

    void clearHistoryFt()
    {
        DynamicList<Tuple2<label,Tuple2<label,vector>>> newHistory;

        forAll (historyFt_,Fti)
        {
            if (historyFt_[Fti].second().first() == 1)
            {
                Tuple2<label,vector> help(-1,historyFt_[Fti].second().second());
                Tuple2<label,Tuple2<label,vector>> help2(historyFt_[Fti].first(),help);
                newHistory.append(help2);
            }
        }

        historyFt_ = newHistory;
    }

    void clearContactInfo()
    {
        isInWallContact_ = false;
        wallContactsVars.clear();
    }

    List<wallContactVars>& getWallContactVar()
    {
        return wallContactsVars;
    }

    void clearCellsList()
    {
        surfCells_[Pstream::myProcNo()].clear();
        intCells_[Pstream::myProcNo()].clear();
    }

    bool cellNotInLists(label cell);

    void appendLists(List<DynamicLabelList> surfCells, List<DynamicLabelList> intCells);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
