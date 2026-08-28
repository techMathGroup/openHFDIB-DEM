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
namespace
{
// Upper bound on the number of sub-volumes a single flood-fill over the
// wall virtual mesh may visit. With the push-guard each sub-volume is
// visited at most once, so a full scan of the matrix is a guaranteed
// terminating bound that never truncates a legitimate search. The bound is
// additionally limited by virtualMeshLevel::maxSubVolumes_ so that an
// accidentally huge virtual mesh cannot monopolise a time step. Computed in
// double to avoid 32-bit label overflow on large matrices.
label maxVSIter(const vector& matrixSize)
{
    scalar nSV =
        matrixSize.x()*matrixSize.y()*matrixSize.z();

    nSV = min(nSV, max(virtualMeshLevel::getMaxSubVolumes(), scalar(1)));

    return label(min(nSV, scalar(labelMax)));
}
}
//---------------------------------------------------------------------------//
void virtualMeshWall::checkAndAppend
(
    const vector& svI,
    DynamicVectorList& auxToCheck
)
{
    vector svIM(svI);
    List<vector> nbrSVI(bbMatrix_.cornerNeighbourSubVolumes(svIM));

    forAll (nbrSVI, nI)
    {
        if (bbMatrix_[nbrSVI[nI]].toCheck)
        {
            auxToCheck.append(nbrSVI[nI]);
        }
    }
}
//---------------------------------------------------------------------------//
void virtualMeshWall::checkAndAppendFace
(
    const vector& svI,
    DynamicVectorList& auxToCheck
)
{
    vector svIM(svI);
    List<vector> nbrSVI(bbMatrix_.faceNeighbourSubVolumes(svIM));

    forAll (nbrSVI, nI)
    {
        if (bbMatrix_[nbrSVI[nI]].toCheck)
        {
            auxToCheck.append(nbrSVI[nI]);
        }
    }
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

    const label iterMax(maxVSIter(bbMatrix_.getMatrixSize()));
    label iterCount(0);

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
            iterCount++;
            checkSubVolume(cSubVolume);

            if (cSubVolume.isCBody)
            {
                vMeshWallInfo_.startingPoint = cSubVolume.center;
                resetSubVolume(cSubVolume);

                return true;

            }
            checkAndAppend(nextToCheck()[sV], auxToCheck());
        }
        autoPtr<DynamicVectorList> helpPtr(nextToCheck.ptr()); // removing const
        nextToCheck.reset(auxToCheck.ptr()); //set -> reset
        auxToCheck = std::move(helpPtr); // adding std::move
        // Only a genuine truncation: cap reached with unvisited sub-volumes
        // still on the frontier. If the frontier is empty the scan completed
        // normally (possibly consuming the whole matrix) and no warning fits.
        if (iterCount >= iterMax && nextToCheck->size() > 0)
        {
            WarningInFunction
                << "virtualMeshWall::detectFirstContactPoint: flood-fill "
                << "visit cap " << iterMax << " reached without finding a "
                << "contact point — no contact is reported this check. "
                << "Virtual mesh bBox: " << bbMatrix_.getBBox()
                << ", matrixSize: " << bbMatrix_.getMatrixSize()
                << ", startingPoint: " << vMeshWallInfo_.getStartingPoint()
                << ". Consider lowering virtualMesh level, increasing "
                << "virtualMesh charCellSize, or raising maxSubVolumes."
                << endl;
            return false;
        }
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

    const label iterMax(maxVSIter(bbMatrix_.getMatrixSize()));
    label iterCount(0);

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
            iterCount++;
            checkSubVolume(cSubVolume);

            if (cSubVolume.isCBody)
            {
                vMeshWallInfo_.startingPoint = cSubVolume.center;
                resetSubVolume(cSubVolume);

                return true;

            }
            checkAndAppendFace(nextToCheck()[sV], auxToCheck());
        }
        autoPtr<DynamicVectorList> helpPtr(nextToCheck.ptr());
        nextToCheck.reset(auxToCheck.ptr());
        auxToCheck = std::move(helpPtr);
        // Only a genuine truncation: cap reached with unvisited sub-volumes
        // still on the frontier (see detectFirstContactPoint).
        if (iterCount >= iterMax && nextToCheck->size() > 0)
        {
            WarningInFunction
                << "virtualMeshWall::detectFirstFaceContactPoint: flood-fill "
                << "visit cap " << iterMax << " reached without finding a "
                << "contact point — no contact is reported this check. "
                << "Virtual mesh bBox: " << bbMatrix_.getBBox()
                << ", matrixSize: " << bbMatrix_.getMatrixSize()
                << ", startingPoint: " << vMeshWallInfo_.getStartingPoint()
                << ". Consider lowering virtualMesh level, increasing "
                << "virtualMesh charCellSize, or raising maxSubVolumes."
                << endl;
            return false;
        }
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

    const label iterMax(maxVSIter(bbMatrix_.getMatrixSize()));

    while (nextToCheck().size() > 0)
    {
        auxToCheck().clear();

        forAll (nextToCheck(),sV)
        {
            subVolumeProperties& cSubVolume = bbMatrix_[nextToCheck()[sV]];
            if (!cSubVolume.toCheck)
            {
                continue;
            }
            iterCount++;

            checkSubVolume(cSubVolume);
            if (cSubVolume.isCBody)
            {
                volumeCount++;
                contactCenter_ += cSubVolume.center;
                checkAndAppendFace(nextToCheck()[sV], auxToCheck());
            }
        }
        autoPtr<DynamicVectorList> helpPtr(nextToCheck.ptr());
        nextToCheck.reset(auxToCheck.ptr());
        auxToCheck = std::move(helpPtr);
        // Only a genuine truncation: cap reached with unvisited sub-volumes
        // still on the frontier (see detectFirstContactPoint).
        if (iterCount >= iterMax && nextToCheck().size() > 0)
        {
            WarningInFunction
                << "virtualMeshWall::evaluateContact: flood-fill visit cap "
                << iterMax << " reached — contact volume is truncated to "
                << volumeCount << " sub-volumes and the reported contact "
                << "force may be underestimated. "
                << "Virtual mesh bBox: " << bbMatrix_.getBBox()
                << ", matrixSize: " << bbMatrix_.getMatrixSize()
                << ", startingPoint: " << vMeshWallInfo_.getStartingPoint()
                << ". Consider lowering virtualMesh level, increasing "
                << "virtualMesh charCellSize, or raising maxSubVolumes."
                << endl;
            break;
        }
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
