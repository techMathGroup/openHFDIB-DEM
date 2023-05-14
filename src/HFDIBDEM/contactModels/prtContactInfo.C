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

#include "interAdhesion.H"
#include "virtualMeshLevel.H"

using namespace Foam;

//---------------------------------------------------------------------------//
prtContactInfo::prtContactInfo
(
    ibContactClass& cClass,
    ibContactVars& cVars,
    ibContactClass& tClass,
    ibContactVars& tVars
)
:
cIbContactClass_(cClass),
cContactVars_(cVars),
tIbContactClass_(tClass),
tContactVars_(tVars)
{
    contactPair_.first() = cVars.bodyId_;
    contactPair_.second() = tVars.bodyId_;
    const materialInfo& cMatInfo(cClass.getMatInfo());
    const materialInfo& tMatInfo(tClass.getMatInfo());

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

    if(interAdhesion::getInterAdhesion().found(adhPotKey))
    {
        adhPot = interAdhesion::getInterAdhesion()[adhPotKey];
    }

    // compute mean model parameters
    physicalProperties_.aY_ = 1/((1 - sqr(cMatInfo.getNu()))/cMatInfo.getY()
        + (1 - sqr(tMatInfo.getNu()))/tMatInfo.getY());
    physicalProperties_.aG_ = 1/(2*(2 - cMatInfo.getNu())*(1 + cMatInfo.getNu())/cMatInfo.getY()
        + 2*(2 - tMatInfo.getNu())*(1 + tMatInfo.getNu())/tMatInfo.getY());
    physicalProperties_.aGammaN_ = (cMatInfo.getGamma()*tMatInfo.getGamma())
        /(cMatInfo.getGamma()+tMatInfo.getGamma()+SMALL);
    physicalProperties_.aGammat_ = (cMatInfo.getGamma()*tMatInfo.getGamma())
        /(cMatInfo.getGamma()+tMatInfo.getGamma()+SMALL);
    physicalProperties_.aMu_ = (cMatInfo.getMu()+tMatInfo.getMu())/2;
    physicalProperties_.maxAdhN_ = cMatInfo.getAdhN() + tMatInfo.getAdhN() - 2*adhPot;
    physicalProperties_.curAdhN_ = 0;
    physicalProperties_.reduceM_ =
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
std::shared_ptr<prtSubContactInfo> prtContactInfo::matchSubContact
(
    boundBox& bbox,
    physicalProperties& physicalProperties,
    Tuple2<label,label>& contactPair
)
{
    for(auto sC : contactList_)
    {
        if (!sC->getVMInfo())
        {
            continue;
        }

        if (bbox.contains(sC->getVMInfo()->getStartingPoint()))
        {
            return sC;
        }
    }

    return std::make_shared<prtSubContactInfo>
        (contactPair, physicalProperties);
}
//---------------------------------------------------------------------------//
void prtContactInfo::limitBBox(boundBox& bbox)
{
    const boundBox cBBox(cIbContactClass_.getGeomModel().getBounds());
    const boundBox tBBox(tIbContactClass_.getGeomModel().getBounds());
    for (label coord = 0; coord < 3; coord++)
    {
        if (bbox.min()[coord] < cBBox.min()[coord])
        {
            bbox.min()[coord] = cBBox.min()[coord];
        }

        if (bbox.max()[coord] > cBBox.max()[coord])
        {
            bbox.max()[coord] = cBBox.max()[coord];
        }

        if (bbox.min()[coord] < tBBox.min()[coord])
        {
            bbox.min()[coord] = tBBox.min()[coord];
        }

        if (bbox.max()[coord] > tBBox.max()[coord])
        {
            bbox.max()[coord] = tBBox.max()[coord];
        }
    }
}
//---------------------------------------------------------------------------//
void prtContactInfo::getContacts_Sphere()
{
    if
    (
        mag(cIbContactClass_.getGeomModel().getCoM()-tIbContactClass_.getGeomModel().getCoM())
        >=
        ((cIbContactClass_.getGeomModel().getDC() / 2) + (tIbContactClass_.getGeomModel().getDC() / 2))
    )
    {
        return;
    }

    newContactList_.emplace_back(std::make_shared<prtSubContactInfo>
        (contactPair_, physicalProperties_)
    );
}
//---------------------------------------------------------------------------//
void prtContactInfo::getContacts_ArbShape
(
    scalar cellV
)
{
    boundBox subCbBox(
        cIbContactClass_.getGeomModel().getBounds().min(),
        cIbContactClass_.getGeomModel().getBounds().max()
    );

    if (!subCbBox.overlaps(tIbContactClass_.getGeomModel().getBounds()))
    {
        return;
    }

    limitBBox(subCbBox);

    scalar charCellSize = pow(cellV,0.333333);
    scalar subVolumeLength = charCellSize/virtualMeshLevel::getLevelOfDivision();
    scalar subVolumeV = pow(subVolumeLength,3);

    newContactList_.emplace_back(matchSubContact(subCbBox, physicalProperties_, contactPair_));
    newContactList_.back()->setVMInfo(subCbBox, subVolumeV);
    return;
}
//---------------------------------------------------------------------------//
void prtContactInfo::swapContactLists()
{
    contactList_.swap(newContactList_);

}
//---------------------------------------------------------------------------//
bool prtContactInfo::contactResolved()
{
    for(auto sC : contactList_)
    {
        if (sC->contactResolved())
        {
            return true;
        }
    }
    return false;
}
//---------------------------------------------------------------------------//
void prtContactInfo::syncContactList()
{
    std::vector<std::shared_ptr<prtSubContactInfo>> syncedContactList;
    // Sync the sub-contact list going through processors
    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        label numOfCnts = 0;
        if (procI == Pstream::myProcNo())
        {
            numOfCnts = contactList_.size();
        }
        reduce(numOfCnts, sumOp<label>());

        for (label i = 0; i < numOfCnts; i++)
        {
            bool vmInfoValid = false;
            if (procI == Pstream::myProcNo())
            {
                vmInfoValid = contactList_[i]->getVMInfo() ? true : false;
            }
            reduce(vmInfoValid, orOp<bool>());

            if (vmInfoValid)
            {
                syncedContactList.emplace_back(std::make_shared<prtSubContactInfo>
                    (contactPair_, physicalProperties_)
                );

                virtualMeshInfo vmInfoToSync;
                if (procI == Pstream::myProcNo())
                {
                    vmInfoToSync = virtualMeshInfo(*(contactList_[i]->getVMInfo()));
                }
                reduce(vmInfoToSync.sV.min(), sumOp<vector>());
                reduce(vmInfoToSync.sV.max(), sumOp<vector>());
                reduce(vmInfoToSync.subVolumeV, sumOp<scalar>());

                point startPointToReduce = vmInfoToSync.getStartingPoint();
                if (procI != Pstream::myProcNo())
                {
                    startPointToReduce = vector::zero;
                }
                reduce(startPointToReduce, sumOp<vector>());

                syncedContactList.back()->setVMInfo(vmInfoToSync);
            }
            else if (procI == 0)
            {
                syncedContactList.emplace_back(std::make_shared<prtSubContactInfo>
                    (contactPair_, physicalProperties_)
                );
            }
        }
    }

    contactList_.swap(syncedContactList);
}
//---------------------------------------------------------------------------//
void prtContactInfo::registerContactList(DynamicList<prtSubContactInfo*>& contactList)
{
    for(auto sC : contactList_)
    {
        contactList.append(sC.get());
    }
}
//---------------------------------------------------------------------------//
void prtContactInfo::clearData()
{
    newContactList_.clear();
    for(auto sC : contactList_)
    {
        sC->clearOutForces();
        sC->setResolvedContact(false);
    }
}
//---------------------------------------------------------------------------//
void prtContactInfo::syncData()
{
    for(auto sC : contactList_)
    {
        sC->syncData();
    }
}
//---------------------------------------------------------------------------//

// ************************************************************************* //
