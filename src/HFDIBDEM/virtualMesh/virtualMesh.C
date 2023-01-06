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
#include "virtualMesh.H"

#include "virtualMeshLevel.H"

using namespace Foam;

//---------------------------------------------------------------------------//
virtualMesh::virtualMesh
(
    virtualMeshInfo& vMeshInfo,
    geomModel& cGeomModel,
    geomModel& tGeomModel
)
:
cGeomModel_(cGeomModel),
tGeomModel_(tGeomModel),
vMeshInfo_(vMeshInfo),
bbMatrix_(vMeshInfo_.subVolumeNVector,
    vMeshInfo_.bBox,
    vMeshInfo_.charCellSize,
    vMeshInfo_.subVolumeV)
{}

virtualMesh::~virtualMesh()
{
}
//---------------------------------------------------------------------------//
bool virtualMesh::detectFirstContactPoint()
{
    autoPtr<DynamicVectorList> nextToCheck(
        new DynamicVectorList);

    autoPtr<DynamicVectorList> auxToCheck(
        new DynamicVectorList);

    nextToCheck->append(bbMatrix_.getSVIndexForPoint(vMeshInfo_.getStartingPoint()));

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

            if (cSubVolume.isContact())
            {
                vMeshInfo_.startingPoint.reset(new point(cSubVolume.center));
                return true;
            }

            auxToCheck().append(bbMatrix_.cornerNeighbourSubVolumes(nextToCheck()[sV]));
            cSubVolume.toCheck = false;
        }
        const autoPtr<DynamicVectorList> helpPtr(nextToCheck.ptr());
        nextToCheck.set(auxToCheck.ptr());
        auxToCheck = helpPtr;
    }
    return false;
}
//---------------------------------------------------------------------------//
scalar virtualMesh::evaluateContact()
{
    label volumeCount(0);
    contactCenter_ = vector::zero;

    autoPtr<DynamicVectorList> nextToCheck(
        new DynamicVectorList);

    autoPtr<DynamicVectorList> auxToCheck(
            new DynamicVectorList);

    nextToCheck->append(bbMatrix_.getSVIndexForPoint(vMeshInfo_.getStartingPoint()));
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

            checkSubVolume(cSubVolume);

            if (cSubVolume.isContact())
            {
                volumeCount++;
                contactCenter_ += cSubVolume.center;
                auxToCheck->append(bbMatrix_.faceNeighbourSubVolumes(nextToCheck()[sV]));
            }
        }
        const autoPtr<DynamicVectorList> helpPtr(nextToCheck.ptr());
        nextToCheck.set(auxToCheck.ptr());
        auxToCheck = helpPtr;
    }

    if (volumeCount > 0)
    {
        contactCenter_ /= volumeCount;
    }

    return volumeCount*bbMatrix_.getSubVolumeV();
}
//---------------------------------------------------------------------------//
void virtualMesh::checkSubVolume(subVolumeProperties& subVolume)
{
    if (subVolume.toCheck)
    {
        subVolume.isCBody = cGeomModel_.pointInside(subVolume.center);
        subVolume.isTBody = tGeomModel_.pointInside(subVolume.center);
        subVolume.toCheck = false;
    }
}
//---------------------------------------------------------------------------//
void virtualMesh::identifySurfaceSubVolumes()
{
    autoPtr<DynamicVectorList> nextToCheck(
        new DynamicVectorList);
    autoPtr<DynamicVectorList> auxToCheck(
        new DynamicVectorList);

    nextToCheck->append(bbMatrix_.getSVIndexForPoint(vMeshInfo_.getStartingPoint()));

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
                if (cSubVolume.isContact())
                {
                    List<vector> neighbourSubVolumes = bbMatrix_.faceNeighbourSubVolumes(nextToCheck()[sV]);
                    if (std::find_if(neighbourSubVolumes.begin(), neighbourSubVolumes.end(),
                        [this](vector const& index)
                        {
                            return !bbMatrix_[index].isContact();

                        }) != neighbourSubVolumes.end())
                    {
                        cSubVolume.isOnEdge = true;
                        if(detectEdgeSubVolumes(nextToCheck()[sV]))
                        {
                            edgeSubVolumesPoints_.append(cSubVolume.center-shiftEdgeSubVolume(nextToCheck()[sV]));
                        }
                    }
                    auxToCheck().append(bbMatrix_.faceNeighbourSubVolumes(nextToCheck()[sV]));
                }
            }
        }
        const autoPtr<DynamicVectorList> helpPtr(nextToCheck.ptr());
        nextToCheck.set(auxToCheck.ptr());
        auxToCheck = helpPtr;
    }
}
//---------------------------------------------------------------------------//
vector virtualMesh::shiftEdgeSubVolume
(
    vector& subVolumeIndex
)
{
    vector subVolumeCenter = bbMatrix_.getPointInMesh(subVolumeIndex);

    vector shiftVector(vector::zero);
    vector incrementVector(vector::zero);

    DynamicVectorList neigboursList;

    neigboursList.append(bbMatrix_.faceNeighbourSubVolumes(subVolumeIndex));
    neigboursList.append(bbMatrix_.edgeNeighbourSubVolumes(subVolumeIndex));
    neigboursList.append(bbMatrix_.cornerNeighbourSubVolumes(subVolumeIndex));

    forAll(neigboursList,nBIter)
    {
        vector subVolumeVector = (subVolumeCenter - bbMatrix_.getPointInMesh(neigboursList[nBIter]));
        subVolumeVector /= mag(subVolumeVector);
        subVolumeVector *= bbMatrix_[neigboursList[nBIter]].isContact();
        incrementVector += subVolumeVector;
    }
    if (mag(incrementVector)>SMALL)
    {
        incrementVector /= mag(incrementVector);
    }
    shiftVector = incrementVector;

    shiftVector *=(bbMatrix_.getCharCellSize()/virtualMeshLevel::getLevelOfDivision())*-0.5;

    return shiftVector;
}
// //---------------------------------------------------------------------------//
bool virtualMesh::detectEdgeSubVolumes
(
    vector& subVolumeIndex
)
{
    bool isOnBorderONE(false);
    bool isOnBorderTWO(false);

    List<vector> neigboursList = bbMatrix_.faceNeighbourSubVolumes(subVolumeIndex);

    forAll(neigboursList,nBIter)
    {
        if(!bbMatrix_[neigboursList[nBIter]].isCBody && bbMatrix_[neigboursList[nBIter]].isTBody
            && !bbMatrix_[neigboursList[nBIter]].isContact())
        {
            isOnBorderONE = true;
        }
        else if(bbMatrix_[neigboursList[nBIter]].isCBody && !bbMatrix_[neigboursList[nBIter]].isTBody
            && !bbMatrix_[neigboursList[nBIter]].isContact())
        {
            isOnBorderTWO = true;
        }

        if(isOnBorderONE && isOnBorderTWO)
        {
            return true;
        }

    }

    if(isOnBorderONE || isOnBorderTWO)
    {
        neigboursList = bbMatrix_.edgeNeighbourSubVolumes(subVolumeIndex);
        forAll(neigboursList,nBIter)
        {
            if(bbMatrix_[neigboursList[nBIter]].toCheck)
            {
                bbMatrix_[neigboursList[nBIter]].toCheck = false;
                vector pointToCheck = bbMatrix_.getPointInMesh(neigboursList[nBIter]);

                bbMatrix_[neigboursList[nBIter]].isCBody = cGeomModel_.pointInside(pointToCheck);
                bbMatrix_[neigboursList[nBIter]].isTBody = tGeomModel_.pointInside(pointToCheck);
            }

            if(!bbMatrix_[neigboursList[nBIter]].isCBody && bbMatrix_[neigboursList[nBIter]].isTBody
                && !bbMatrix_[neigboursList[nBIter]].isContact())
            {
                isOnBorderONE = true;
            }
            else if(bbMatrix_[neigboursList[nBIter]].isCBody && !bbMatrix_[neigboursList[nBIter]].isTBody
                && !bbMatrix_[neigboursList[nBIter]].isContact())
            {
                isOnBorderTWO = true;
            }
            if(isOnBorderONE && isOnBorderTWO)
            {
                return true;
            }
        }
    }
    return false;
}
//---------------------------------------------------------------------------//
Tuple2<scalar,vector> virtualMesh::get3DcontactNormalAndSurface()
{
// This function is taken from prtContact and just adjustated for higher accuracy
    scalar area(0.0);
    vector normalVec(vector::zero);
    scalar tDC(tGeomModel_.getDC());
    vector normalVector = vector::zero;
    point closestPoint = vector::zero;
    tGeomModel_.getClosestPointAndNormal
    (
        contactCenter_,
        vector::one * tDC,
        closestPoint,
        normalVector
    );

    if (edgeSubVolumesPoints_.size() >= 3)
    {
        bool normOk(false);
        vector center(vector::zero);

        forAll (edgeSubVolumesPoints_,cell)
        {
            center += edgeSubVolumesPoints_[cell];
        }
        center /= edgeSubVolumesPoints_.size();

        scalar xx(0);
        scalar xy(0);
        scalar xz(0);
        scalar yy(0);
        scalar yz(0);
        scalar zz(0);

        forAll (edgeSubVolumesPoints_,cell)
        {
            vector subPoint(edgeSubVolumesPoints_[cell] - center);
            if(subPoint != vector::zero)
                subPoint = subPoint/mag(subPoint);
            xx += subPoint[0] * subPoint[0];
            xy += subPoint[0] * subPoint[1];
            xz += subPoint[0] * subPoint[2];
            yy += subPoint[1] * subPoint[1];
            yz += subPoint[1] * subPoint[2];
            zz += subPoint[2] * subPoint[2];
        }

        xx /= edgeSubVolumesPoints_.size();
        xy /= edgeSubVolumesPoints_.size();
        xz /= edgeSubVolumesPoints_.size();
        yy /= edgeSubVolumesPoints_.size();
        yz /= edgeSubVolumesPoints_.size();
        zz /= edgeSubVolumesPoints_.size();

        vector weightedDir(vector::zero);

        scalar detX(yy*zz-yz*yz);
        vector axisDirX(detX,xz*yz-xy*zz,xy*yz-xz*yy);
        scalar weightX(detX*detX);
        if((weightedDir & axisDirX) < 0.0)
            weightX = -weightX;
        weightedDir += axisDirX * weightX;

        scalar detY(xx*zz-xz*xz);
        vector axisDirY(xz*yz-xy*zz,detY,xy*xz-yz*xx);
        scalar weightY(detY*detY);
        if((weightedDir & axisDirY) < 0.0)
            weightY = -weightY;
        weightedDir += axisDirY * weightY;

        scalar detZ(xx*yy-xy*xy);
        vector axisDirZ(xy*yz-xz*yy,xy*xz-yz*xx,detZ);
        scalar weightZ(detZ*detZ);
        if((weightedDir & axisDirZ) < 0.0)
            weightZ = -weightZ;
        weightedDir += axisDirZ * weightZ;

        if(mag(weightedDir) > SMALL)
        {
            normOk = true;
            normalVec = weightedDir/mag(weightedDir);
        }
        if (!normOk || mag(normalVec) < 1)
            normalVec = normalVector;

        // create best fitting plane
        plane bestFitPlane(contactCenter_, normalVec);
        normalVec = bestFitPlane.normal();
        DynamicPointList commCellsPosInPlane;
        forAll (edgeSubVolumesPoints_,cell)
        {
            commCellsPosInPlane.append(bestFitPlane.nearestPoint(edgeSubVolumesPoints_[cell]));
		}

        vector q1(1.0, 0.0, 0.0);
        vector q2(0.0, 1.0, 0.0);
        if (abs(q1 & bestFitPlane.normal()) > abs(q2 & bestFitPlane.normal()))
            q1 = q2;

        vector u(bestFitPlane.normal() ^ q1);
        vector v(bestFitPlane.normal() ^ u);

        DynamicList<plane> clockwisePlanes;
        List<scalar> helpList(6);
        helpList[0] = 0.0;
        helpList[1] = 1.0;
        helpList[2] = 0.0;
        helpList[3] = -1.0;
        helpList[4] = 0.0;
        helpList[5] = 1.0;

        // loop over parts of plane to find and sort points
        DynamicVectorList commCellsInSections;
        for (label i = 0; i < 4; i++)
        //~ for (label i = 0; i < 8; i++)
        {
            scalar uStep(helpList[i + 1] -  helpList[i]);
            scalar vStep(helpList[i + 2] -  helpList[i + 1]);
            DynamicPointList pointsInSection;
            //~ for (scalar j = 0.0; j < 5.0; j += 1.0)
            for (scalar j = 0.0; j < 15.0; j += 1.0)
            {
                plane uPlane(contactCenter_, u*(helpList[i] + uStep*j/15.0) + v*(helpList[i + 1] + vStep*j/15.0));
                plane vPlane(contactCenter_, u*(helpList[i] + uStep*(j+1)/15.0) + v*(helpList[i + 1] + vStep*(j+1)/15.0));

                forAll (commCellsPosInPlane, celli)
                {
                    if (uPlane.sideOfPlane(commCellsPosInPlane[celli]) == 0 && vPlane.sideOfPlane(commCellsPosInPlane[celli]) == 1)
                    {
                        pointsInSection.append(commCellsPosInPlane[celli]);
                    }
                }

                if (pointsInSection.size() > SMALL)
                {
                    vector max(vector::zero);
                    scalar distance(0);
                    forAll (pointsInSection, pointI)
                    {
                         scalar magnitude = mag(pointsInSection[pointI]-contactCenter_);
                         if (magnitude > distance)
                         {
							 magnitude = distance;
							 max = pointsInSection[pointI];
						 }
					}
                    commCellsInSections.append(max);
                }
            }
        }
        // calculate contact area
        for (label i = 0; i + 1 < commCellsInSections.size(); i++)
        {
            //~ Info << commCellsInSections[i] << endl;
            vector AC(commCellsInSections[i] - contactCenter_);
            vector BC(commCellsInSections[i + 1] - contactCenter_);
            vector crossPr(AC ^ BC);
            area += mag(crossPr)/2;
        }

        if (commCellsInSections.size() > 2)
        {
            vector AC(commCellsInSections[commCellsInSections.size() - 1] - contactCenter_);
            vector BC(commCellsInSections[0] - contactCenter_);
            vector crossPr(AC ^ BC);
            area += mag(crossPr)/2;
        }
    }

    if (normalVec == vector::zero)
    {
        normalVec = normalVector;
	}

	if ((normalVector & normalVec) < 0)
	{
		normalVec = normalVec * (-1);
	}

    Tuple2<scalar,vector> returnValue(area,normalVec);

    return returnValue;
}
//---------------------------------------------------------------------------//
// ************************************************************************* //
