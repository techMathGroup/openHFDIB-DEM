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
    prtContactInfo.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef prtContactInfo_H
#define prtContactInfo_H

#include "dictionary.H"
#include "fvCFD.H"

#include "geomModel.H"
#include "ibContactVars.H"
#include "ibContactClass.H"
#include "prtSubContactInfo.H"
#include "physicalProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class prtContactInfo Declaration
\*---------------------------------------------------------------------------*/

class prtContactInfo
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

    // bool isFileCreated_ = true;

    protected:

    ibContactClass& cIbContactClass_;
    ibContactVars& cContactVars_;

    ibContactClass& tIbContactClass_;
    ibContactVars& tContactVars_;

    Tuple2<label,label> contactPair_;

    std::vector<std::shared_ptr<prtSubContactInfo>> contactList_;
    std::vector<std::shared_ptr<prtSubContactInfo>> newContactList_;

    physicalProperties physicalProperties_;

    std::shared_ptr<prtSubContactInfo> matchSubContact
    (
        boundBox& bbox,
        physicalProperties& physicalProperties,
        Tuple2<label,label>& contactPair
    );

    void limitBBox(boundBox& bbox);

    public:

    // Constructors ------------------------------------------------
    prtContactInfo
    (
        ibContactClass& cClass,
        ibContactVars& cVars,
        ibContactClass& tClass,
        ibContactVars& tVars
    );

    // Destructors -------------------------------------------------
    ~prtContactInfo();
    //-----------------------------Func-------------------------------//

    const Tuple2<label, label>& getCPair()
    {
        return contactPair_;
    }

    ibContactClass& getcClass()
    {
        return cIbContactClass_;
    }

    ibContactClass& gettClass()
    {
        return tIbContactClass_;
    }

    ibContactVars& getcVars()
    {
        return cContactVars_;
    }

    ibContactVars& gettVars()
    {
        return tContactVars_;
    }

    scalar getMu()
    {
        return physicalProperties_.aMu_;
    }

    void getContacts_Sphere();

    void getContacts_ArbShape
    (
        scalar cellV
    );

    void swapContactLists();

    bool contactResolved();

    std::vector<std::shared_ptr<prtSubContactInfo>>& getPrtSCList()
    {
        return contactList_;
    }

    void syncContactList();

    void registerContactList(DynamicList<prtSubContactInfo*>& contactList);

    void clearData();

    void syncData();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
