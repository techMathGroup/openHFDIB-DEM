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
#include "prtSubContactInfo.H"

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
vector prtSubContactInfo::getLVec(const point& CoM)
{
    return prtCntVars_.contactCenter_ - CoM;
}
//---------------------------------------------------------------------------//
vector prtSubContactInfo::getVeli(ibContactVars& cVars, vector& lVec)
{
    return (-((lVec-cVars.Axis_*((lVec) & cVars.Axis_))
        ^ cVars.Axis_)*cVars.omega_+ cVars.Vel_);
}
//---------------------------------------------------------------------------//
void prtSubContactInfo::evalVariables(
    const point& cCoM,
    const point& tCoM,
    ibContactVars& cVars,
    ibContactVars& tVars
)
{
    cLVec_ = getLVec(cCoM);
    tLVec_ = getLVec(tCoM);

    cVeli_ = getVeli(cVars, cLVec_);
    tVeli_ = getVeli(tVars, tLVec_);

    Vn_ = -(cVeli_ - tVeli_) & prtCntVars_.contactNormal_;
    Lc_ = 4*mag(cLVec_)*mag(tLVec_)/(mag(cLVec_) + mag(tLVec_));

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
    return (physicalProperties_.aGammaN_*sqrt(physicalProperties_.aY_
        *physicalProperties_.reduceM_/pow(Lc_+SMALL,3))*
        (prtCntVars_.contactArea_ * Vn_))*prtCntVars_.contactNormal_;
}
//---------------------------------------------------------------------------//
vector prtSubContactInfo::getFt(scalar deltaT)
{
    // compute relative tangential velocity
    vector cVeliNorm = cVeli_ - ((cVeli_ & prtCntVars_.contactNormal_)
        *prtCntVars_.contactNormal_);

    vector tVeliNorm = tVeli_ - ((tVeli_ & prtCntVars_.contactNormal_)
        *prtCntVars_.contactNormal_);

    vector Vt(cVeliNorm - tVeliNorm);
    // compute tangential force
    vector Ftdi(- physicalProperties_.aGammat_
        *sqrt(physicalProperties_.aG_*physicalProperties_.reduceM_*Lc_)*Vt);
    Ftdi -= physicalProperties_.aG_*Lc_*Vt*deltaT;
    return Ftdi;
}
//---------------------------------------------------------------------------//
void prtSubContactInfo::setVMInfo(boundBox& bBox, scalar subVolumeV)
{
    if (!vmInfo_)
    {
        vmInfo_ = std::make_shared<virtualMeshInfo>(bBox, subVolumeV);
        return;
    }

    vmInfo_->sV = subVolume(bBox);
    vmInfo_->subVolumeV = subVolumeV;
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
    vmInfo_->startingPoint = vmInfo.startingPoint;
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
