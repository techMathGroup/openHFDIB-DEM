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
#include "prtContactInfo.H"

using namespace Foam;

//---------------------------------------------------------------------------//
prtContactInfo::prtContactInfo
(
    ibContactClass& cClass,
    ibContactVars& cVars,
    ibContactClass& tClass,
    ibContactVars& tVars,
    HashTable<scalar,string,Hash<string>>& matInterAdh
)
:
cIbContactClass_(cClass),
cContactVars_(cVars),
tIbContactClass_(tClass),
tContactVars_(tVars),
FtPrev_(vector::zero)
{
    contactPair_.first() = cVars.bodyId_;
    contactPair_.second() = tVars.bodyId_;
    materialInfo& cMatInfo(cClass.getMatInfo());
    materialInfo& tMatInfo(tClass.getMatInfo());

    scalar adhPot = 0;
    string adhPotKey;

    if(cMatInfo.getMaterial() < tMatInfo.getMaterial())
    {
        adhPotKey = cMatInfo.getMaterial();
        adhPotKey += "-";
        adhPotKey += tMatInfo.getMaterial();
    }
    else
    {
        adhPotKey = tMatInfo.getMaterial();
        adhPotKey += "-";
        adhPotKey += cMatInfo.getMaterial();
    }
    if(matInterAdh.found(adhPotKey))
    {
        adhPot = matInterAdh[adhPotKey];
    }

    // compute mean model parameters
    aY_ = 1/((1 - sqr(cMatInfo.getNu()))/cMatInfo.getY()
        + (1 - sqr(tMatInfo.getNu()))/tMatInfo.getY());
    aG_ = 1/(2*(2 - cMatInfo.getNu())*(1 + cMatInfo.getNu())/cMatInfo.getY()
        + 2*(2 - tMatInfo.getNu())*(1 + tMatInfo.getNu())/tMatInfo.getY());
    aGammaN_ = aY_*(cMatInfo.getGamma()*tMatInfo.getGamma())
        /(cMatInfo.getGamma()+tMatInfo.getGamma()+SMALL);
    aGammat_ = aG_*(cMatInfo.getGamma()*tMatInfo.getGamma())
        /(cMatInfo.getGamma()+tMatInfo.getGamma()+SMALL);
    aMu_ = (cMatInfo.getMu()+tMatInfo.getMu())/2;
    maxAdhN_ = cMatInfo.getAdhN() + tMatInfo.getAdhN() - 2*adhPot;
    curAdhN_ = 0;

    reduceM_ =
    (
        cIbContactClass_.getGeomModel().getM0()
        *tIbContactClass_.getGeomModel().getM0()
        /(cIbContactClass_.getGeomModel().getM0()
        +tIbContactClass_.getGeomModel().getM0())
    );
}

prtContactInfo::~prtContactInfo()
{
}
//---------------------------------------------------------------------------//
vector prtContactInfo::getLVec(ibContactClass& contactClass)
{
    return contactClass.getGeomModel().getLVec(prtCntVars_.contactCenter_);
}
//---------------------------------------------------------------------------//
vector prtContactInfo::getVeli(ibContactVars& cVars, vector& lVec)
{
    return (-((lVec-cVars.Axis_*((lVec) & cVars.Axis_))
        ^ cVars.Axis_)*cVars.omega_+ cVars.Vel_);
}
//---------------------------------------------------------------------------//
void prtContactInfo::evalVariables()
{
    cLVec_ = getLVec(cIbContactClass_);
    InfoH << DEM_Info << "-- Particle-particle contact cLVec_ " << cLVec_ << endl;
    tLVec_ = getLVec(tIbContactClass_);
    InfoH << DEM_Info << "-- Particle-particle contact tLVec_ " << tLVec_ << endl;

    cVeli_ = getVeli(cContactVars_, cLVec_);
    InfoH << DEM_Info << "-- Particle-particle contact cVeli_ " << cVeli_ << endl;
    tVeli_ = getVeli(tContactVars_, tLVec_);
    InfoH << DEM_Info << "-- Particle-particle contact tVeli_ " << tVeli_ << endl;

    Vn_ = -(cVeli_ - tVeli_) & prtCntVars_.contactNormal_;
    InfoH << DEM_Info << "-- Particle-particle contact Vn_ " << Vn_ << endl;
    Lc_ = 4*mag(cLVec_)*mag(tLVec_)/(mag(cLVec_) + mag(tLVec_));
    InfoH << DEM_Info << "-- Particle-particle contact Lc_ " << Lc_ << endl;

    curAdhN_ = min
    (
        maxAdhN_,
        max(curAdhN_, aY_*prtCntVars_.contactVolume_/(sqr(Lc_)*8
            *Foam::constant::mathematical::pi))
    );
}
//---------------------------------------------------------------------------//
vector prtContactInfo::getFNe()
{
    return (aY_*prtCntVars_.contactVolume_/(Lc_+SMALL))
        *prtCntVars_.contactNormal_;
}
//---------------------------------------------------------------------------//
vector prtContactInfo::getFA()
{
    return ((sqrt(8*Foam::constant::mathematical::pi*aY_
        *curAdhN_*prtCntVars_.contactVolume_))
        *prtCntVars_.contactNormal_);
}
//---------------------------------------------------------------------------//
vector prtContactInfo::getFNd()
{
    return (aGammaN_*sqrt(aY_*reduceM_/pow(Lc_+SMALL,3))*
        (prtCntVars_.contactArea_ * Vn_))*prtCntVars_.contactNormal_;
}
//---------------------------------------------------------------------------//
vector prtContactInfo::getFt(scalar deltaT)
{
    // project last Ft into a new direction
    vector FtLastP(FtPrev_ - (FtPrev_ & prtCntVars_.contactNormal_)
        *prtCntVars_.contactNormal_);
    // scale projected Ft to have same magnitude as FtLast
    vector FtLastS(mag(FtPrev_) * (FtLastP/(mag(FtLastP)+SMALL)));
    // compute relative tangential velocity
    vector cVeliNorm = cVeli_ - ((cVeli_ & prtCntVars_.contactNormal_)
        *prtCntVars_.contactNormal_);

    vector tVeliNorm = tVeli_ - ((tVeli_ & prtCntVars_.contactNormal_)
        *prtCntVars_.contactNormal_);

    vector Vt(cVeliNorm - tVeliNorm);
    // compute tangential force
    vector Ftdi(- aGammat_*sqrt(aG_*reduceM_*Lc_)*Vt);
    FtPrev_ = FtLastS - aG_*Lc_*Vt*deltaT + Ftdi;
    return FtPrev_;
}
//---------------------------------------------------------------------------//

// ************************************************************************* //
