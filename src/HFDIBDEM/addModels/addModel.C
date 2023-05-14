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
#include "addModel.H"

using namespace Foam;

//---------------------------------------------------------------------------//
addModel::~addModel()
{
}
//---------------------------------------------------------------------------//
bool addModel::isBodyInContact(PtrList<immersedBody>& immersedBodies)
{
    std::shared_ptr<geomModel> gModel = geomModel_->getCopy();
    ibContactClass cIbClass(
        gModel,
        "None"
    );

    vector velAxi (vector::zero);
    scalar helpScalar(0.0);

    ibContactVars cIbVars(
        immersedBodies.size(),
        velAxi,
        helpScalar,
        velAxi,
        helpScalar,
        helpScalar,
        geomModel_->getRhoS()
    );

    wallContactInfo cIBWallCntI(
            cIbClass,
            cIbVars
    );

    bool inContact = contactModel::detectWallContact(
        mesh_,
        cIbClass,
        cIBWallCntI
    );

    if(!inContact)
    {
        forAll(immersedBodies, ibI)
        {
            bool bBoxContact = true;
            boundBox ibIbBox = immersedBodies[ibI].getGeomModel().getBounds();
            boundBox cbBox = geomModel_->getBounds();

            forAll(geometricD,dir)
            {
                if(geometricD[dir] == 1)
                {
                    if(!(ibIbBox.max()[dir] >= cbBox.min()[dir]
                        && ibIbBox.min()[dir] <= cbBox.max()[dir]))
                    {
                        bBoxContact = false;
                        break;
                    }
                }
            }

            if(!bBoxContact)
            {
                continue;
            }

            prtContactInfo prtCInfo (
                cIbClass,
                cIbVars,
                immersedBodies[ibI].getibContactClass(),
                immersedBodies[ibI].getContactVars()
            );

            contactModel::getContacts(
                mesh_,
                prtCInfo
            );

            DynamicList<prtSubContactInfo*> sCList;
            prtCInfo.registerContactList(sCList);
            forAll(sCList,sC)
            {
                prtSubContactInfo* prtSCInfo = sCList[sC];

                if(contactModel::detectPrtPrtContact(
                    mesh_,
                    cIbClass,
                    immersedBodies[ibI].getibContactClass(),
                    *prtSCInfo
                ))
                {
                    inContact = true;
                    break;
                }

            }

        }
    }

    return inContact;
}
