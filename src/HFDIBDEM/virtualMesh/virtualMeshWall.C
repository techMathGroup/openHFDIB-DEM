/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| "_ \ / _ \ "_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
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
#include "virtualMeshWall.H"

#include "virtualMeshLevel.H"

using namespace Foam;

//---------------------------------------------------------------------------//
virtualMeshWall::virtualMeshWall
(
    virtualMeshWallInfo& vMeshWallInfo,
    geomModel& cGeomModel
)
:
cGeomModel_(cGeomModel),
vMeshWallInfo_(vMeshWallInfo),
bbMatrix_(vMeshWallInfo.subVolumeNVector,
    vMeshWallInfo.bBox,
    vMeshWallInfo.charCellSize,
    vMeshWallInfo.subVolumeV)
{}

virtualMeshWall::~virtualMeshWall()
{
}
//---------------------------------------------------------------------------//
bool virtualMeshWall::detectFirstContactPoint()
{
    autoPtr<DynamicVectorList> nextToCheck(
        new DynamicVectorList);

    autoPtr<DynamicVectorList> auxToCheck(
        new DynamicVectorList);

    nextToCheck->append(bbMatrix_.getSVIndexForPoint_Wall(vMeshWallInfo_.getStartingPoint()));
    nextToCheck->append(bbMatrix_.cornerNeighbourSubVolumes(nextToCheck()[0]));
    // InfoH << DEM_Info << " -- VM firstSV : " << nextToCheck()[0] << " point " << bbMatrix_[nextToCheck()[0]].center << endl;
    while (nextToCheck->size() > 0)
    {
        auxToCheck->clear();
        forAll (nextToCheck(),sV)
        {
            subVolumeProperties& cSubVolume = bbMatrix_[nextToCheck()[sV]];
            if (!cSubVolume.toCheck)
            {
                continue;
            }
            checkSubVolume(cSubVolume);

            if (cSubVolume.isCBody)
            {
                vMeshWallInfo_.startingPoint = cSubVolume.center;
                resetSubVolume(cSubVolume);

                return true;

            }
            auxToCheck().append(bbMatrix_.cornerNeighbourSubVolumes(nextToCheck()[sV]));
        }
        autoPtr<DynamicVectorList> helpPtr(nextToCheck.ptr()); // removing const
        nextToCheck.reset(auxToCheck.ptr()); //set -> reset
        auxToCheck = std::move(helpPtr); // adding std::move
    }
    return false;
}
//---------------------------------------------------------------------------////---------------------------------------------------------------------------//
bool virtualMeshWall::detectFirstFaceContactPoint()
{
    autoPtr<DynamicVectorList> nextToCheck(
        new DynamicVectorList);

    autoPtr<DynamicVectorList> auxToCheck(
        new DynamicVectorList);

    nextToCheck->append(bbMatrix_.getSVIndexForPoint_Wall(vMeshWallInfo_.getStartingPoint()));
    nextToCheck->append(bbMatrix_.faceNeighbourSubVolumes(nextToCheck()[0]));
    // InfoH << DEM_Info << " -- VM firstSV : " << nextToCheck()[0] << " point " << bbMatrix_[nextToCheck()[0]].center << endl;
    while (nextToCheck->size() > 0)
    {
        auxToCheck->clear();
        forAll (nextToCheck(),sV)
        {
            subVolumeProperties& cSubVolume = bbMatrix_[nextToCheck()[sV]];
            if (!cSubVolume.toCheck)
            {
                continue;
            }
            checkSubVolume(cSubVolume);

            if (cSubVolume.isCBody)
            {
                vMeshWallInfo_.startingPoint = cSubVolume.center;
                resetSubVolume(cSubVolume);

                return true;

            }
            auxToCheck().append(bbMatrix_.faceNeighbourSubVolumes(nextToCheck()[sV]));
        }
        autoPtr<DynamicVectorList> helpPtr(nextToCheck.ptr());
        nextToCheck.reset(auxToCheck.ptr());
        auxToCheck = std::move(helpPtr);
    }
    return false;
}
//---------------------------------------------------------------------------//
scalar virtualMeshWall::evaluateContact()
{
    label volumeCount = 0;
    contactCenter_ = vector::zero;
    autoPtr<DynamicVectorList> nextToCheck(
        new DynamicVectorList);
    autoPtr<DynamicVectorList> auxToCheck(
        new DynamicVectorList);
    nextToCheck->append(bbMatrix_.getSVIndexForPoint_Wall(vMeshWallInfo_.getStartingPoint()));
    label iterCount(0);
    while (nextToCheck().size() > 0)
    {
        auxToCheck().clear();

        forAll (nextToCheck(),sV)
        {
            subVolumeProperties& cSubVolume = bbMatrix_[nextToCheck()[sV]];
            iterCount++;
            if (!cSubVolume.toCheck)
            {
                continue;
            }

            checkSubVolume(cSubVolume);
            if (cSubVolume.isCBody)
            {
                volumeCount++;
                contactCenter_ += cSubVolume.center;
                auxToCheck->append(bbMatrix_.faceNeighbourSubVolumes(nextToCheck()[sV]));
            }
        }
        autoPtr<DynamicVectorList> helpPtr(nextToCheck.ptr());
        nextToCheck.reset(auxToCheck.ptr());
        auxToCheck = std::move(helpPtr);
    }
    if (volumeCount > 0)
    {
        contactCenter_ /= volumeCount;
    }

    return volumeCount*bbMatrix_.getSubVolumeV();
}
//---------------------------------------------------------------------------//
void virtualMeshWall::checkSubVolume(subVolumeProperties& subVolume)
{
    if (subVolume.toCheck)
    {
        subVolume.isCBody = cGeomModel_.pointInside(subVolume.center);
        subVolume.toCheck = false;
    }
}
//---------------------------------------------------------------------------//
void virtualMeshWall::resetSubVolume(subVolumeProperties& subVolume)
{
    subVolume.toCheck = true;
    subVolume.isCBody = false;
    subVolume.isTBody = false;
    subVolume.isOnEdge = false;   
}
//---------------------------------------------------------------------------//
label virtualMeshWall::getInternalSV()
{
    autoPtr<DynamicVectorList> nextToCheck(
        new DynamicVectorList);
    autoPtr<DynamicVectorList> auxToCheck(
        new DynamicVectorList);

    nextToCheck->append(bbMatrix_.getSVIndexForPoint(vMeshWallInfo_.getStartingPoint()));
    label innerSVCount(0);
    vectorHashSet octreeSvSet;
    
    while (nextToCheck->size() > 0)
    {
        auxToCheck().clear();
        forAll (nextToCheck(),sV)
        {   
            if (!octreeSvSet.found(nextToCheck()[sV]))
            {
                octreeSvSet.insert(nextToCheck()[sV]);
                subVolumeProperties& cSubVolume = bbMatrix_[nextToCheck()[sV]];
                if (cSubVolume.isCBody)
                {
                    bool isNotOnEdge(true);
                    List<vector> neighbourSubVolumes = bbMatrix_.faceNeighbourSubVolumes(nextToCheck()[sV]);
                    neighbourSubVolumes.append(bbMatrix_.edgeNeighbourSubVolumes(nextToCheck()[sV]));
                    neighbourSubVolumes.append(bbMatrix_.cornerNeighbourSubVolumes(nextToCheck()[sV]));
                    forAll(neighbourSubVolumes,nSV)
                    {
                        //~ isNotOnEdge *= bbMatrix_[neighbourSubVolumes[nSV]].isCBody;//OF.com issues warning
                        isNotOnEdge &= bbMatrix_[neighbourSubVolumes[nSV]].isCBody;
                    }
                    if(!isNotOnEdge)
                    {
                        cSubVolume.isOnEdge = true;
                    }
                    auxToCheck().append(bbMatrix_.faceNeighbourSubVolumes(nextToCheck()[sV]));
                }
                if(!cSubVolume.isOnEdge)
                {
                    innerSVCount++;
                }
            }
        }

        autoPtr<DynamicVectorList> helpPtr(nextToCheck.ptr());
        nextToCheck.reset(auxToCheck.ptr());
        auxToCheck = std::move(helpPtr);
    }
    return innerSVCount;
}
// ************************************************************************* //
