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
#include "wallContactInfo.H"

using namespace Foam;

//---------------------------------------------------------------------------//
wallContactInfo::wallContactInfo
(
    ibContactClass& cClass,
    ibContactVars& cVars,
    HashTable<materialInfo,string,Hash<string>>& wInfos,
    HashTable<scalar,string,Hash<string>> matInterAdh
)
:
ibContactClass_(cClass),
ibContactVars_(cVars)
{
    bodyId_ = cVars.bodyId_;
    materialInfo& cMatInfo(cClass.getMatInfo());

    List<string> cntPatches(wInfos.toc());
    forAll(cntPatches, patchI)
    {
        materialInfo& wInfo(wInfos[cntPatches[patchI]]);

        scalar adhPot = 0;
        string adhPotKey;

        if(cMatInfo.getMaterial() < wInfo.getMaterial())
        {
            adhPotKey = cMatInfo.getMaterial();
            adhPotKey += "-";
            adhPotKey += wInfo.getMaterial();
        }
        else
        {
            adhPotKey = wInfo.getMaterial();
            adhPotKey += "-";
            adhPotKey += cMatInfo.getMaterial();
        }
        if(matInterAdh.found(adhPotKey))
        {
            adhPot = matInterAdh[adhPotKey];
        }

        // compute mean model parameters
        scalar aY = 1/((1 - sqr(cMatInfo.getNu()))/cMatInfo.getY()
            + (1 - sqr(wInfo.getNu()))/wInfo.getY());
        scalar aG = 1/(2*(2 - cMatInfo.getNu())*(1 + cMatInfo.getNu())
            /cMatInfo.getY() + 2*(2 - wInfo.getNu())
            *(1 + wInfo.getNu())/wInfo.getY());
        scalar aGammaN = aY*(cMatInfo.getGamma()*wInfo.getGamma())
            /(cMatInfo.getGamma()+wInfo.getGamma()+SMALL);
        scalar aGammat = aG*(cMatInfo.getGamma()*wInfo.getGamma())
            /(cMatInfo.getGamma()+wInfo.getGamma()+SMALL);
        scalar aMu = (cMatInfo.getMu()+wInfo.getMu())/2;
        scalar maxAdhN = cMatInfo.getAdhN() + wInfo.getAdhN() - 2*adhPot;

        wallMeanPars_.insert(
            cntPatches[patchI],
            meanContactPar(aY, aG, aGammaN, aGammat, aMu, maxAdhN)
        );
    }

    reduceM_ = 0;
}

wallContactInfo::~wallContactInfo()
{
}
//---------------------------------------------------------------------------//
wallContactVars& wallContactInfo::getWallCntVars
(
    const vector contactCenter,
    const vector contactNormal,
    const scalar deltaT
)
{
    forAll(wallCntHashes_, hashI)
    {
        wallContactVars& wallCntVarsI(wallcVarsTable_[wallCntHashes_[hashI]]);
        if(
            mag(wallCntVarsI.contactCenter_ - contactCenter)
            <
            2 * mag(wallCntVarsI.Veli_) * deltaT
            &&
            mag(wallCntVarsI.contactNormal_ - contactNormal) < 1e-4
        )
        {
            curUsedHashes_.append(wallCntHashes_[hashI]);
            return wallCntVarsI;
        }
    }

    label newHash = wallCntHashes_.size();

    wallcVarsTable_.insert(newHash, wallContactVars());
    wallCntHashes_.append(newHash);
    curUsedHashes_.append(wallCntHashes_.last());
    return wallcVarsTable_[newHash];
}
//---------------------------------------------------------------------------//
void wallContactInfo::clearWallCntVars()
{
    forAll(wallCntHashes_, hashI)
    {
        bool used = false;
        forAll(curUsedHashes_, hashII)
        {
            if(wallCntHashes_[hashI] == curUsedHashes_[hashII])
            {
                used = true;
                break;
            }
        }

        if(!used)
        {
            if(wallcVarsTable_.found(wallCntHashes_[hashI]))
            {
                wallcVarsTable_.erase(wallCntHashes_[hashI]);
            }
        }
    }

    wallCntHashes_ = curUsedHashes_;
    curUsedHashes_.clear();
}
//---------------------------------------------------------------------------//
vector wallContactInfo::getLVec(wallContactVars& wallCntvar)
{
    return ibContactClass_.getGeomModel().getLVec(wallCntvar.contactCenter_);
}
//---------------------------------------------------------------------------//
vector wallContactInfo::getVeli(wallContactVars& wallCntvar)
{
    return (-((wallCntvar.lVec_-ibContactVars_.Axis_
        *((wallCntvar.lVec_) & ibContactVars_.Axis_))
        ^ ibContactVars_.Axis_)*ibContactVars_.omega_+ ibContactVars_.Vel_);
}
//---------------------------------------------------------------------------//
void wallContactInfo::evalVariables(wallContactVars& wallCntvar)
{
    reduceM_ =
    (
        ibContactClass_.getGeomModel().getM0()
        *ibContactClass_.getGeomModel().getM0()
        /(ibContactClass_.getGeomModel().getM0()
        +ibContactClass_.getGeomModel().getM0())
    );

    wallCntvar.lVec_ = getLVec(wallCntvar);

    wallCntvar.Veli_ = getVeli(wallCntvar);

    wallCntvar.Vn_ = -(wallCntvar.Veli_ - vector::zero) & wallCntvar.contactNormal_;
    wallCntvar.Lc_ = 4*mag(wallCntvar.lVec_)*mag(wallCntvar.lVec_)
        /(mag(wallCntvar.lVec_) + mag(wallCntvar.lVec_));

    wallCntvar.curAdhN_ = min
    (
        wallCntvar.getMeanCntPar().maxAdhN_,
        max(wallCntvar.curAdhN_, wallCntvar.getMeanCntPar().aY_
            *wallCntvar.contactVolume_
            /(sqr(wallCntvar.Lc_)*8*Foam::constant::mathematical::pi))
    );
}
//---------------------------------------------------------------------------//
vector wallContactInfo::getFNe(wallContactVars& wallCntvar)
{
    return (wallCntvar.getMeanCntPar().aY_*wallCntvar.contactVolume_
        /(wallCntvar.Lc_+SMALL))*wallCntvar.contactNormal_;
}
//---------------------------------------------------------------------------//
vector wallContactInfo::getFA(wallContactVars& wallCntvar)
{
    return ((sqrt(8*Foam::constant::mathematical::pi
        *wallCntvar.getMeanCntPar().aY_
        *wallCntvar.curAdhN_*wallCntvar.contactVolume_))
        *wallCntvar.contactNormal_);
}
//---------------------------------------------------------------------------//
vector wallContactInfo::getFNd(wallContactVars& wallCntvar)
{
    meanContactPar& meanCntPar(wallCntvar.getMeanCntPar());
    return (meanCntPar.aGammaN_*sqrt(meanCntPar.aY_*reduceM_
        /pow(wallCntvar.Lc_+SMALL,3))*(wallCntvar.contactArea_*wallCntvar.Vn_))
        *wallCntvar.contactNormal_;
}
//---------------------------------------------------------------------------//
vector wallContactInfo::getFt(wallContactVars& wallCntvar, scalar deltaT)
{
    meanContactPar& meanCntPar(wallCntvar.getMeanCntPar());
    // project last Ft into a new direction
    vector FtLastP(wallCntvar.FtPrev_
        - (wallCntvar.FtPrev_ & wallCntvar.contactNormal_)
        *wallCntvar.contactNormal_);
    // scale projected Ft to have same magnitude as FtLast
    vector FtLastS(mag(wallCntvar.FtPrev_) * (FtLastP/(mag(FtLastP)+SMALL)));
    // compute relative tangential velocity
    vector cVeliNorm = wallCntvar.Veli_
        - ((wallCntvar.Veli_ & wallCntvar.contactNormal_)
        *wallCntvar.contactNormal_);

    vector Vt(cVeliNorm - vector::zero);
    // compute tangential force
    vector Ftdi(- meanCntPar.aGammat_*sqrt(meanCntPar.aG_*reduceM_*wallCntvar.Lc_)*Vt);
    wallCntvar.FtPrev_ = FtLastS - meanCntPar.aG_*wallCntvar.Lc_*Vt*deltaT + Ftdi;
    return wallCntvar.FtPrev_;
}
//---------------------------------------------------------------------------//

// ************************************************************************* //
