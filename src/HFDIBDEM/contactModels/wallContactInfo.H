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
    wallContactInfo.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef wallContactInfo_H
#define wallContactInfo_H

#include "dictionary.H"
#include "fvCFD.H"

#include "geomModel.H"
#include "ibContactVars.H"
#include "ibContactClass.H"
#include "wallContactVars.H"
#include "materialInfo.H"
#include "physicalProperties.H"
#include "forces.H"
#include "spectatorMesh.H"
#include "virtualMeshWall.H"
#include "wallSubContactInfo.H"
#include "wallPlaneInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class wallContactInfo Declaration
\*---------------------------------------------------------------------------*/

class wallContactInfo
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;
    typedef HashSet<vector,Hash<vector>> vectorHashSet;

    ibContactClass& ibContactClass_;
    ibContactVars& ibContactVars_;

    label bodyId_;

    // Mean model parameters
    HashTable<physicalProperties,string,Hash<string>> wallMeanPars_;

    std::vector<std::shared_ptr<wallSubContactInfo>> subCList_;

    scalar reduceM_;    // Reduced mass

    forces outForce_; // Out variable for force

    autoPtr<spectatorMesh> SM_;

    autoPtr<virtualMeshWallInfo> vmwInfo_;

    List<string> contactPatches_;

    public:

    // Constructors ------------------------------------------------
    wallContactInfo
    (
        ibContactClass& cClass,
        ibContactVars& cVars
    );

    // Destructors -------------------------------------------------
    ~wallContactInfo();

    HashTable<physicalProperties,string,Hash<string>>& getWallMeanPars()
    {
        return wallMeanPars_;
    }

    ibContactClass& getcClass()
    {
        return ibContactClass_;
    }

    ibContactVars& getcVars()
    {
        return ibContactVars_;
    }

    label getBodyId()
    {
        return bodyId_;
    }

    boundBox constructBoundBox
    (
        boundBox& bodyBB
    );

    void constructSM();

    bool isInsidePlane
    (
        vector checkedPoint,
        const string& wall,
        const HashTable<List<vector>,string,Hash<string>>& wallPatches = wallPlaneInfo::getWallPlaneInfo()
    );

    bool detectWallContact(
        const HashTable<List<vector>,string,Hash<string>>& wallPatches = wallPlaneInfo::getWallPlaneInfo()
    );

    void findContactAreas();

    List<DynamicList<vector>> detectPossibleSMContact
    (
        DynamicVectorList& contactPoints,
        List<string>& contactPatches
    );

    boundBox getSCBBox
    (
        DynamicVectorList& subContactAreas
    );

    boundBox constructVMBox
    (
        boundBox& baseContactAreaBB,
        string& wallName
    );

    vector getPlanePoint
    (
        vector pointInDomain,
        string wallName
    );

    void checkSMElement
    (
        vector& index,
        List<string> contactPatches
    );

    void checkElement
    (
        vector& index,
        bool& inMesh,
        bool& inBody,
        bool& allInMesh,
        bool& allInBody
    );

    boundBox contactPlaneBBox
    (
        boundBox contactBoundBox,
        string wallName
    );

    boundBox correctSMBBforWall
    (
        boundBox bB,
        string wallName
    );

    void clearOldContact();

    void setNewSubContact(
        const List<Tuple2<point,boundBox>>& contactBBData,
        const List<Tuple2<point,boundBox>>& planeBBData,
        const List<string>& contactPatches,
        const List<Tuple2<point,boundBox>>& isInternal,
        boundBox BB
    );

    std::vector<std::shared_ptr<wallSubContactInfo>>& getWallSCList()
    {
        return subCList_;
    }

    void registerSubContactList(DynamicList<wallSubContactInfo*>& contactList);

    List<string>& getContactPatches()
    {
        return contactPatches_;
    }

    // void addSCToEnd
    // (
    //     std::vector<std::shared_ptr<wallSubContactInfo>>& sCW
    // )
    // {
    //     subCList_.push_back(sCW);
    // }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
