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

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "wallSubContactInfo.H"

#include "interAdhesion.H"
#include "wallMatInfo.H"

#include "virtualMeshLevel.H"
#include "wallPlaneInfo.H"
#include "contactModelInfo.H"
using namespace Foam;
//---------------------------------------------------------------------------//
wallSubContactInfo::wallSubContactInfo
(
    List<Tuple2<point,boundBox>> contactBBData,
    List<Tuple2<point,boundBox>> planeBBData,
    List<string> contactPatches,
    List<Tuple2<point,boundBox>> internalBBData,
    HashTable<physicalProperties,string,Hash<string>> wallMeanPars,
    boundBox BB,
    label bodyId
)
:
contactPatches_(contactPatches),
internalBBData_(internalBBData),
wallMeanPars_(wallMeanPars),
BB_(BB),
bodyId_(bodyId)
{
    forAll(contactBBData,cBD)
    {
        vector subVolumeNVector = vector(
            floor((contactBBData[cBD].second().span()[0]/virtualMeshLevel::getCharCellSize())*virtualMeshLevel::getLevelOfDivision()),
            floor((contactBBData[cBD].second().span()[1]/virtualMeshLevel::getCharCellSize())*virtualMeshLevel::getLevelOfDivision()),
            floor((contactBBData[cBD].second().span()[2]/virtualMeshLevel::getCharCellSize())*virtualMeshLevel::getLevelOfDivision())
        );
        if(cmptMin(subVolumeNVector)<SMALL)
        {
            // Pout <<" Trubble with subVolumeNVector "<< subVolumeNVector << endl;
            for(int i=0;i<3;i++)
            {
                if(subVolumeNVector[i] <SMALL)
                {                
                    subVolumeNVector[i] = 1;
                    contactBBData[cBD].second().min()[i] -=virtualMeshLevel::getCharCellSize()/virtualMeshLevel::getLevelOfDivision()*0.5;
                    contactBBData[cBD].second().max()[i] +=virtualMeshLevel::getCharCellSize()/virtualMeshLevel::getLevelOfDivision()*0.5;
                }  
            }
            // Pout <<" Corrected subVolumeNVector "<< subVolumeNVector << endl;
        }
     

        autoPtr<virtualMeshWallInfo> vmWInfo(
            new virtualMeshWallInfo(
                contactBBData[cBD].second(),
                contactBBData[cBD].first(),
                subVolumeNVector,
                virtualMeshLevel::getCharCellSize(),
                pow(virtualMeshLevel::getCharCellSize()/virtualMeshLevel::getLevelOfDivision(),3)
            )  
        );
        vmWInfoList_.append(std::move(vmWInfo));
    }

    forAll(planeBBData,pBD)
    {
        vector subVolumeNVector = vector(
            ceil((planeBBData[pBD].second().span()[0]/virtualMeshLevel::getCharCellSize()))*virtualMeshLevel::getLevelOfDivision(),
            ceil((planeBBData[pBD].second().span()[1]/virtualMeshLevel::getCharCellSize()))*virtualMeshLevel::getLevelOfDivision(),
            ceil((planeBBData[pBD].second().span()[2]/virtualMeshLevel::getCharCellSize()))*virtualMeshLevel::getLevelOfDivision()
        );

        for(int i=0;i<3;i++)
        {
            if(subVolumeNVector[i] == planeBBData[pBD].second().minDim())
            {                
                subVolumeNVector[i] = 1;
                planeBBData[pBD].second().min()[i] -=virtualMeshLevel::getCharCellSize()/virtualMeshLevel::getLevelOfDivision()*0.5;
                planeBBData[pBD].second().max()[i] +=virtualMeshLevel::getCharCellSize()/virtualMeshLevel::getLevelOfDivision()*0.5;
            }  
        }

        autoPtr<virtualMeshWallInfo> vmWInfo(
            new virtualMeshWallInfo(
                planeBBData[pBD].second(),
                planeBBData[pBD].first(),
                subVolumeNVector,
                virtualMeshLevel::getCharCellSize(),
                pow(virtualMeshLevel::getCharCellSize()/virtualMeshLevel::getLevelOfDivision(),3)
            )
        );
        vmPlaneInfoList_.append(std::move(vmWInfo));
    }
}

wallSubContactInfo::~wallSubContactInfo()
{
}
//---------------------------------------------------------------------------//
vector wallSubContactInfo::getLVec(wallContactVars& wallCntvar, ibContactClass ibCClass)
{
    return ibCClass.getGeomModel().getLVec(wallCntvar.contactCenter_);
}
//---------------------------------------------------------------------------//
vector wallSubContactInfo::getVeli(wallContactVars& wallCntvar, ibContactVars& cVars)
{
    return (-((wallCntvar.lVec_-cVars.Axis_
        *((wallCntvar.lVec_) & cVars.Axis_))
        ^ cVars.Axis_)*cVars.omega_+ cVars.Vel_);
}
//---------------------------------------------------------------------------//
void wallSubContactInfo::evalVariables(
    wallContactVars& wallCntvar,
    ibContactClass& ibCClass,
    ibContactVars& cVars
)
{
    reduceM_ =
    (
        ibCClass.getGeomModel().getM0()
        *ibCClass.getGeomModel().getM0()
        /(ibCClass.getGeomModel().getM0()
        +ibCClass.getGeomModel().getM0())
    );

    wallCntvar.lVec_ = getLVec(wallCntvar,ibCClass);
    // wallCntvar.lVec_ = wallCntvar.contactCenter_ - ibCClass.getGeomModel().getCoM();
    wallCntvar.Veli_ = getVeli(wallCntvar, cVars);

    wallCntvar.Vn_ = -(wallCntvar.Veli_ - vector::zero) & wallCntvar.contactNormal_;
    wallCntvar.Lc_ = (contactModelInfo::getLcCoeff())*mag(wallCntvar.lVec_)*mag(wallCntvar.lVec_)/(mag(wallCntvar.lVec_) + mag(wallCntvar.lVec_));
    

    wallCntvar.curAdhN_ = min
    (
        wallCntvar.getMeanCntPar().maxAdhN_,
        max(wallCntvar.curAdhN_, wallCntvar.getMeanCntPar().aY_
            *wallCntvar.contactVolume_
            /(sqr(wallCntvar.Lc_)*8*Foam::constant::mathematical::pi))
    );
}
//---------------------------------------------------------------------------//
vector wallSubContactInfo::getFNe(wallContactVars& wallCntvar)
{
    return (wallCntvar.getMeanCntPar().aY_*wallCntvar.contactVolume_
        /(wallCntvar.Lc_+SMALL))*wallCntvar.contactNormal_;
}
//---------------------------------------------------------------------------//
vector wallSubContactInfo::getFA(wallContactVars& wallCntvar)
{
    return ((sqrt(8*Foam::constant::mathematical::pi
        *wallCntvar.getMeanCntPar().aY_
        *wallCntvar.curAdhN_*wallCntvar.contactVolume_))
        *wallCntvar.contactNormal_);
}
//---------------------------------------------------------------------------//
vector wallSubContactInfo::getFNd(wallContactVars& wallCntvar)
{
    physicalProperties& meanCntPar(wallCntvar.getMeanCntPar());
    return ((meanCntPar.reduceBeta_*sqrt(meanCntPar.aY_
            *reduceM_*wallCntvar.contactArea_/(wallCntvar.Lc_+SMALL))*
            wallCntvar.Vn_)*wallCntvar.contactNormal_);

}
//---------------------------------------------------------------------------//
vector wallSubContactInfo::getFt(wallContactVars& wallCntvar, scalar deltaT)
{
    physicalProperties& meanCntPar(wallCntvar.getMeanCntPar());
    // project last Ft into a new direction
    vector FtLastP(wallCntvar.FtPrev_
        - (wallCntvar.FtPrev_ & wallCntvar.contactNormal_)
        *wallCntvar.contactNormal_);
    // scale projected Ft to have same magnitude as FtLast
    vector FtLastS(mag(wallCntvar.FtPrev_) * (FtLastP/(mag(FtLastP)+SMALL)));
    // compute relative tangential velocity
    // vector cVeliNorm = wallCntvar.Veli_
        // - ((wallCntvar.Veli_ & wallCntvar.contactNormal_)
        // *wallCntvar.contactNormal_);
    vector cVeliNorm = wallCntvar.Veli_*(wallCntvar.Veli_&wallCntvar.contactNormal_);

    vector Vt(wallCntvar.Veli_-(cVeliNorm - vector::zero));
    // compute tangential force
    if(contactModelInfo::getUseMindlinRotationalModel())
    {
        
        scalar kT = 200*8*meanCntPar.aG_*(wallCntvar.contactArea_/(wallCntvar.Lc_+SMALL));
        vector deltaFt(kT*Vt*deltaT + 2*meanCntPar.reduceBeta_*sqrt(kT*reduceM_)*Vt);
        wallCntvar.FtPrev_ = - FtLastS - deltaFt;
    }

    if(contactModelInfo::getUseChenRotationalModel())
    {
   
        vector Ftdi(meanCntPar.reduceBeta_*sqrt(meanCntPar.aG_*reduceM_*wallCntvar.Lc_)*Vt);
        Ftdi += meanCntPar.aG_*wallCntvar.Lc_*Vt*deltaT;
        wallCntvar.FtPrev_ = - FtLastS - Ftdi;
    }

    return wallCntvar.FtPrev_;
}
//---------------------------------------------------------------------------//
void wallSubContactInfo::syncData()
{
    reduce(outForce_.F, sumOp<vector>());
    reduce(outForce_.T, sumOp<vector>());
}
//---------------------------------------------------------------------------//
void wallSubContactInfo::syncContactResolve()
{
    reduce(contactResolved_,orOp<bool>());
}
//---------------------------------------------------------------------------//
autoPtr<virtualMeshWallInfo>& wallSubContactInfo::getVMContactInfo
(
    label ID
)
{
    return vmWInfoList_[ID];
}
//---------------------------------------------------------------------------//
autoPtr<virtualMeshWallInfo>& wallSubContactInfo::getVMPlaneInfo
(
    label ID
)
{
    return vmPlaneInfoList_[ID];
}

//---------------------------------------------------------------------------//
