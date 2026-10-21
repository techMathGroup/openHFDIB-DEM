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

    // psp: geometric floor on Lc. As a thin body buries itself in the
    // target, the contact center converges onto its CoM and the raw
    // harmonic Lc collapses to ~0, so FNe = aY*VC/Lc diverges although
    // VC grows only marginally. Floor Lc by the smaller body's
    // half-thickness (smallest bounding-box extent) so deep-penetration
    // stiffness stays bounded: force then tracks the overlap volume.
    {
        const boundBox cBBox(cIb.getGeomModel().getBounds());
        const boundBox tBBox(tIb.getGeomModel().getBounds());
        const vector cSpan(cBBox.span());
        const vector tSpan(tBBox.span());
        const scalar cMinSpan(min(min(cSpan.x(), cSpan.y()), cSpan.z()));
        const scalar tMinSpan(min(min(tSpan.x(), tSpan.y()), tSpan.z()));
        const scalar LcFloor(contactModelInfo::getLcCoeff()*min(cMinSpan, tMinSpan));

        if (Lc_ < LcFloor)
        {
            Lc_ = LcFloor;
        }
    }

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
vector prtSubContactInfo::getFNSl()
{
    const surfaceLoadingLaw::result r
    (
        surfaceLoadingLaw::evaluate
        (
            physicalProperties_.sl_,
            Vn_,
            prtCntVars_.contactArea_,
            prtCntVars_.contactVolume_,
            Lc_
        )
    );

    prtCntVars_.surfaceLoading_ = true;
    prtCntVars_.slPMean_ = r.pMean;
    prtCntVars_.slPStr_ = r.pStr;
    prtCntVars_.slPEdge_ = r.pEdge;
    prtCntVars_.slPHug_ = r.pHug;
    prtCntVars_.slBlendF_ = r.blendF;
    prtCntVars_.slZEff_ = r.zEff;
    prtCntVars_.slW_ = r.w;
    prtCntVars_.slVn_ = Vn_;

    return r.F*prtCntVars_.contactNormal_;
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
vector prtSubContactInfo::getFt(scalar deltaT, scalar FtCeil)
{
    return contactModelInfo::tangentialSpringForce
    (
        FtPrev_,
        cVeli_ - tVeli_,
        prtCntVars_.contactNormal_,
        prtCntVars_.contactArea_,
        Lc_,
        physicalProperties_.aG_,
        physicalProperties_.reduceBeta_,
        physicalProperties_.reduceM_,
        deltaT,
        FtCeil
    );
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
