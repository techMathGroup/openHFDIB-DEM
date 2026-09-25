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
    wallSubContactInfo.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef wallSubContactInfo_H
#define wallSubContactInfo_H

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class wallSubContactInfo Declaration
\*---------------------------------------------------------------------------*/

class wallSubContactInfo
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;
    // list of collision pathes 
    List<string> contactPatches_;

    //bool list for internal SM elements
    List<Tuple2<point,boundBox>> internalBBData_;

    // Mean model parameters
    HashTable<physicalProperties,string,Hash<string>> wallMeanPars_;

    boundBox BB_;

    label bodyId_;

    wallContactVars wallCntvar_;

    scalar reduceM_;    // Reduced mass

    forces outForce_; // Out variable for force

    List<autoPtr<virtualMeshWallInfo>> vmWInfoList_; 

    List<autoPtr<virtualMeshWallInfo>> vmPlaneInfoList_;

    bool contactResolved_ = false;

    public:

    // Constructors ------------------------------------------------
    wallSubContactInfo
    (
    List<Tuple2<point,boundBox>> contactBBData,
    List<Tuple2<point,boundBox>> planeBBData,
    List<string> contactPatches,
    List<Tuple2<point,boundBox>> internalBBData,
    HashTable<physicalProperties,string,Hash<string>> wallMeanPars,
    boundBox BB,
    label bodyId
    );

    // Destructors -------------------------------------------------
    ~wallSubContactInfo();

    label getVMContactSize()
    {
        return vmWInfoList_.size();
    }
    
    label getVMPlaneSize()
    {
        return vmPlaneInfoList_.size();
    }

    autoPtr<virtualMeshWallInfo>& getVMContactInfo
    (
        label ID
    );

    autoPtr<virtualMeshWallInfo>& getVMPlaneInfo
    (
        label ID
    );

    wallContactVars& getWallCntVars()
    {
        return wallCntvar_;
    }

    void clearWallCntVars();

    scalar getMu(wallContactVars& wallCntvar)
    {
        return wallCntvar.getMeanCntPar().aMu_;
    }

    vector getLVec
    (
        wallContactVars& wallCntvar,
        ibContactClass ibCClass
    );

    vector getVeli
    (
        wallContactVars& wallCntvar,
        ibContactVars& cVars
    );

    void evalVariables
    (
        wallContactVars& wallCntvar,
        ibContactClass& ibCClass,
        ibContactVars& cVars
    );

    vector getFNe(wallContactVars& wallCntvar);

    vector getFA(wallContactVars& wallCntvar);

    vector getFNd(wallContactVars& wallCntvar);

    vector getFt(wallContactVars& wallCntvar, scalar deltaT);

    forces& getOutForce()
    {
        return outForce_;
    }

    void setResolvedContact(bool resolved)
    {
        contactResolved_ = resolved;
    }

    bool getContactResolved()
    {
        return contactResolved_;
    }

    bool isContactPatch(string patchName)
    {
        return wallMeanPars_.found(patchName);
    }

    HashTable<physicalProperties,string,Hash<string>>& getWallMeanPars()
    {
        return wallMeanPars_;
    }

    List<autoPtr<virtualMeshWallInfo>>& getVMContactInfo()
    {
        return vmWInfoList_;
    }

    List<autoPtr<virtualMeshWallInfo>>& getVMPlaneInfo()
    {
        return vmPlaneInfoList_;
    }

    const List<Tuple2<point,boundBox>>& getInternalElements()
    {
        return internalBBData_;
    }

    const List<string>& getContactPatches()
    {
        return contactPatches_;
    }

    void syncData();

    void syncContactResolve();

    label& getBodyId()
    {
        return bodyId_;
    }

    boundBox& getsWCBB()
    {
        return BB_;
    }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //