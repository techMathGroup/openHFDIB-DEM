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

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/
#include "prtSubContactInfo.H"
#include "contactModelInfo.H"

using namespace Foam;

//---------------------------------------------------------------------------//
prtSubContactInfo::prtSubContactInfo
(
    const Tuple2<label,label>& contactPair,
    const physicalProperties& physicalProperties
)
:
contactPair_(contactPair),
physicalProperties_(physicalProperties)
{}

prtSubContactInfo::~prtSubContactInfo()
{}
//---------------------------------------------------------------------------//
vector prtSubContactInfo::getVeli(ibContactVars& cVars, vector& lVec)
{
    return (-((lVec-cVars.Axis_*((lVec) & cVars.Axis_))
        ^ cVars.Axis_)*cVars.omega_+ cVars.Vel_);
}
//---------------------------------------------------------------------------//
void prtSubContactInfo::evalVariables(
    ibContactClass& cIb,
    ibContactClass& tIb,
    ibContactVars& cVars,
    ibContactVars& tVars
)
{
    cLVec_ = cIb.getGeomModel().getLVec(prtCntVars_.contactCenter_);
    // cLVec_ = prtCntVars_.contactCenter_ - cIb.getGeomModel().getCoM();
    tLVec_ = tIb.getGeomModel().getLVec(prtCntVars_.contactCenter_);
    // tLVec_ = prtCntVars_.contactCenter_ - tIb.getGeomModel().getCoM();

    cVeli_ = getVeli(cVars, cLVec_);
    tVeli_ = getVeli(tVars, tLVec_);

    Vn_ = -(cVeli_ - tVeli_) & prtCntVars_.contactNormal_;
    Lc_ = (contactModelInfo::getLcCoeff())*mag(cLVec_)*mag(tLVec_)/(mag(cLVec_) + mag(tLVec_));

    physicalProperties_.curAdhN_ = min
    (
        physicalProperties_.maxAdhN_,
        max(physicalProperties_.curAdhN_, physicalProperties_.aY_*prtCntVars_.contactVolume_/(sqr(Lc_)*8
            *Foam::constant::mathematical::pi))
    );
}
//---------------------------------------------------------------------------//
vector prtSubContactInfo::getFNe()
{
    return (physicalProperties_.aY_*prtCntVars_.contactVolume_/(Lc_+SMALL))
        *prtCntVars_.contactNormal_;
}
//---------------------------------------------------------------------------//
vector prtSubContactInfo::getFA()
{
    return ((sqrt(8*Foam::constant::mathematical::pi*physicalProperties_.aY_
        *physicalProperties_.curAdhN_*prtCntVars_.contactVolume_))
        *prtCntVars_.contactNormal_);
}
//---------------------------------------------------------------------------//
vector prtSubContactInfo::getFNd()
{
    return (physicalProperties_.reduceBeta_*sqrt(physicalProperties_.aY_
        *physicalProperties_.reduceM_*prtCntVars_.contactArea_/(Lc_+SMALL))*
        Vn_)*prtCntVars_.contactNormal_;

}
//---------------------------------------------------------------------------//
vector prtSubContactInfo::getFt(scalar deltaT)
{
    // compute relative tangential velocity
    vector FtLastP(FtPrev_ - (FtPrev_ & prtCntVars_.contactNormal_)
        *prtCntVars_.contactNormal_);

    // scale projected Ft to have same magnitude as FtLast


    vector FtLastS(mag(FtPrev_) * (FtLastP/(mag(FtLastP)+SMALL)));

    // compute relative tangential velocity
    vector relVeli(cVeli_ - tVeli_);
    vector veliNomr((relVeli)*(relVeli & prtCntVars_.contactNormal_));
    vector Vt(relVeli-veliNomr);
    // compute tangential force
        //NewDefinition
    if(contactModelInfo::getUseMindlinRotationalModel())
    {

        // tangTune(200) was an empirical tuning factor (cc92e39; Ondra),
        // not from mindlin theory - see demTimeStepInfo::tangTune_
        scalar tangTune(1.0);
        scalar kT = tangTune*8*physicalProperties_.aG_*(prtCntVars_.contactArea_/(Lc_+SMALL));
        vector deltaFt(kT*Vt*deltaT + 2*physicalProperties_.reduceBeta_*sqrt(kT*physicalProperties_.reduceM_)*Vt);
        FtPrev_ = - FtLastS - deltaFt;
    }

    if(contactModelInfo::getUseChenRotationalModel())
    {
        
        vector Ftdi(- physicalProperties_.reduceBeta_*sqrt(physicalProperties_.aG_*physicalProperties_.reduceM_*Lc_)*Vt);
        Ftdi += physicalProperties_.aG_*Lc_*Vt*deltaT;
        FtPrev_ = - FtLastS- Ftdi;
    }

    
    return FtPrev_;
}
//---------------------------------------------------------------------------//
void prtSubContactInfo::setVMInfo
(
    boundBox& bBox,
    scalar subVolumeV,
    scalar emptyScale
)
{
    if (!vmInfo_)
    {
        vmInfo_ = std::make_shared<virtualMeshInfo>(bBox, subVolumeV, emptyScale);
        return;
    }

    vmInfo_->sV = subVolume(bBox);
    vmInfo_->subVolumeV = subVolumeV;
    vmInfo_->emptyScale = emptyScale;
}
//---------------------------------------------------------------------------//
void prtSubContactInfo::setVMInfo(const virtualMeshInfo& vmInfo)
{
    if (!vmInfo_)
    {
        vmInfo_ = std::make_shared<virtualMeshInfo>(vmInfo);
        return;
    }

    vmInfo_->sV = vmInfo.sV;
    vmInfo_->subVolumeV = vmInfo.subVolumeV;
    vmInfo_->emptyScale = vmInfo.emptyScale;
    // vmInfo_->startingPoint = std::move(vmInfo.startingPoint);
    vmInfo_->startingPoint.reset(new point(*vmInfo.startingPoint));
}
//---------------------------------------------------------------------------//
std::shared_ptr<virtualMeshInfo>& prtSubContactInfo::getVMInfo()
{
    return vmInfo_;
}
//---------------------------------------------------------------------------//
void prtSubContactInfo::syncData()
{
    reduce(outForce_.first().F, sumOp<vector>());
    reduce(outForce_.first().T, sumOp<vector>());
    reduce(outForce_.second().F, sumOp<vector>());
    reduce(outForce_.second().T, sumOp<vector>());


    if (vmInfo_)
    {
        point reducePoint = vector::zero;

        if (contactResolved_)
        {
            reducePoint = vmInfo_->getStartingPoint();
        }

        reduce(reducePoint, sumOp<vector>());
        vmInfo_->startingPoint.reset(new point(reducePoint));
    }
}
//---------------------------------------------------------------------------//


// ************************************************************************* //
