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

Description
    Algorithmic details of the particle-wall virtual mesh.

    Flood-fill structure (all traversals)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Breadth-first frontier expansion over the bbMatrix lattice:

      nextToCheck : sub-volumes of the current frontier
      auxToCheck  : unvisited (toCheck == true) neighbours gathered
                    for the next frontier

      while (nextToCheck not empty)
          for each sv in frontier (skipping already-visited ones)
              checkSubVolume(sv)          // centroid pointInside test,
                                          // toCheck := false
              if (sv.isCBody) -> contact found / volume counted
              append unvisited neighbours of sv to auxToCheck
          swap frontier <-> auxToCheck

    The two autoPtr-frontier lists are swapped (not copied) between
    iterations; the toCheck flag on each sub-volume guarantees single
    visits and thus termination within one lattice scan.

    Neighbourhoods: corner (26-neighbourhood) expansion is used by
    detectFirstContactPoint to locate any contact point quickly;
    face (6-neighbourhood) expansion is used by
    detectFirstFaceContactPoint and evaluateContact to stay inside the
    connected contact patch.

    First contact point (detectFirstContactPoint)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The seed is the sub-volume containing the starting point plus its
    corner neighbours. The first classified-inside sub-volume
    terminates the search; its centre becomes the new starting point
    (persisted by the caller for the next time step) and the volume is
    reset so a subsequent evaluation re-floods from a clean state.

    Contact evaluation (evaluateContact)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Floods from the (updated) starting point with face neighbours,
    counting classified-inside sub-volumes and accumulating their
    centres into the (unweighted) contact centre. The returned
    volume is volumeCount*subVolumeV. Because only inside
    sub-volumes spread the frontier, the flood is confined to the
    connected patch containing the seed.

    Exact second pass (virtualMesh dict entry exactSubVolume, set
    wall planes only): after the flood every face neighbour of an
    inside sub-volume is classified, so the contact band (inside
    sub-volumes with an outside or missing face neighbour) is read
    off directly. Band leaves are clipped by the dominant body
    surface plane (geomModel::getBoxSurfacePlane) and the fluid
    side of every wall half-space (planePolyClip, exact); interior
    leaves keep their full sub-volume; the contact centre becomes
    volume-weighted. Fallbacks: a band leaf without a reliable body
    plane keeps its full sub-volume; a clipped-to-zero total falls
    back to the legacy count. Plane virtual meshes (wall-contact
    area evaluation) carry no wall planes and stay legacy.

    Wetted-area evaluation (evaluateContactAreaExact, plane
    virtual meshes): per band leaf of the preceding
    evaluateContact flood, the polygon where the kept material
    (body side of the surface plane, behind the wall) meets the
    wall plane is the wetted footprint of that leaf; its polygon
    area is the exact contribution. Interior leaves and leaves
    without a reliable body plane keep their full face area; a
    band leaf whose kept region does not reach the wall plane
    wets nothing. Semantics: the area is the planar cross-section
    of the contact patch on the wall (the sphere analytic path
    measures pi*contactRad^2 likewise), not the curved cap.

    Overflow cap (maxVSIter, from virtualMeshTools)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    iterMax = min(nLatticeSubVolumes, maxSubVolumes), computed in
    double to avoid 32-bit label overflow. It is checked after each
    frontier swap: if the cap is reached while the frontier is still
    non-empty, the traversal was genuinely truncated and a warning is
    issued (detection: no contact reported; evaluation: truncated
    volume, force possibly underestimated). An empty frontier at the
    cap means the lattice was scanned completely, which is legitimate
    and silent.

    Sketch of the frontier expansion
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        frontier:      +----------------------------+
                      | sv | sv | sv | ...          |  checkSubVolume
                      +----------------------------+  (pointInside test)
                            | face/corner neighbours
                            v
        next frontier: unvisited (toCheck) neighbours only

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/
#include "virtualMeshWall.H"

#include "virtualMeshTools.H"

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
bbMatrix_
(
    vMeshWallInfo.subVolumeNVector,
    vMeshWallInfo.bBox,
    vMeshWallInfo.charCellSize,
    vMeshWallInfo.subVolumeV
)
{}

virtualMeshWall::~virtualMeshWall()
{
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

    insideLeaves_.clear();

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
                insideLeaves_.append(nextToCheck()[sV]);
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
            // the inside-leaf set is incomplete: keep the exact area
            // pass from walking a partial band (it refuses on the
            // empty set and the caller keeps its legacy area)
            insideLeaves_.clear();
            break;
        }
    }
    if (volumeCount > 0)
    {
        contactCenter_ /= volumeCount;
    }

    // exact second pass: the flood has classified every face neighbour
    // of every inside leaf, so the contact band (inside leaves with an
    // outside or missing face neighbour) can be read off directly.
    // band leaves are clipped by the body surface plane and every wall
    // half-space; interior leaves keep their full sub-volume. plane
    // virtual meshes carry no wall planes and stay on the legacy count
    if (    volumeCount > 0
        &&  virtualMeshLevel::getExactSubVolume()
        &&  wallPlanes_.size() > 0)
    {
        scalar exactVolume(0);
        vector exactCenter(vector::zero);

        forAll(insideLeaves_,lI)
        {
            subVolumeProperties& cSubVolume = bbMatrix_[insideLeaves_[lI]];
            vector svI(insideLeaves_[lI]);

            if (isOnContactBand(svI))
            {
                evaluateContactLeafExact(cSubVolume);
                exactVolume += cSubVolume.clippedVolume_;
                exactCenter +=
                    cSubVolume.clippedVolume_*cSubVolume.clippedCentroid_;
            }
            else
            {
                exactVolume += bbMatrix_.getSubVolumeV();
                exactCenter += bbMatrix_.getSubVolumeV()*cSubVolume.center;
            }
        }

        if (exactVolume > VSMALL)
        {
            contactCenter_ = exactCenter/exactVolume;

            return exactVolume;
        }

        // detected-but-clipped-to-zero: the flood counted inside
        // leaves but every band clip removed its region, so a plane
        // estimate failed on the band rather than the contact being
        // empty. fall back to the legacy count (contactCenter_
        // already carries the unweighted mean)
        // Note (MI): see how often this fires, maybe, we shall
        //            debug levels and write this only for debug levels > 1
        WarningInFunction
            << "virtualMeshWall::evaluateContact: the flood found "
            << volumeCount << " inside sub-volumes but the exact band "
            << "clips summed to zero -- falling back to the legacy "
            << "sub-volume count. "
            << "Virtual mesh bBox: " << bbMatrix_.getBBox()
            << ", matrixSize: " << bbMatrix_.getMatrixSize()
            << ", startingPoint: " << vMeshWallInfo_.getStartingPoint()
            << endl;
    }

    return volumeCount*bbMatrix_.getSubVolumeV();
}
//---------------------------------------------------------------------------//
void virtualMeshWall::evaluateContactLeafExact
(
    subVolumeProperties& subVolume
)
{
    if (subVolume.clippedVolume_ > -0.5)
    {
        // already evaluated (a previous evaluateContact this check)
        return;
    }

    subVolume.clippedVolume_ = 0;
    subVolume.clippedCentroid_ = vector::zero;

    const boundBox leaf(leafBox(subVolume));

    point bP(vector::zero);
    vector bN(vector::zero);

    if (cGeomModel_.getBoxSurfacePlane(leaf, bP, bN))
    {
        // kept material = body side of the surface plane and the
        // fluid side of every wall half-space
        List<planePolyClip::halfSpace> planes(wallPlanes_);
        planes.append(planePolyClip::halfSpace(bP, bN));

        planePolyClip::volumeAndCentroid
        (
            leaf,
            planes,
            subVolume.clippedVolume_,
            subVolume.clippedCentroid_
        );
    }
    else
    {
        // no reliable body plane in this leaf: full sub-volume
        subVolume.clippedVolume_ = bbMatrix_.getSubVolumeV();
        subVolume.clippedCentroid_ = subVolume.center;
    }
}
//---------------------------------------------------------------------------//
bool virtualMeshWall::isOnContactBand
(
    vector& svI
)
{
    List<vector> nbrSVI(bbMatrix_.faceNeighbourSubVolumes(svI));

    if (nbrSVI.size() < 6)
    {
        // missing face neighbour (out of matrix): the seam is
        // interior to the contact, treat it as outside
        return true;
    }

    forAll(nbrSVI,nSV)
    {
        if (!bbMatrix_[nbrSVI[nSV]].isCBody)
        {
            return true;
        }
    }

    return false;
}
//---------------------------------------------------------------------------//
boundBox virtualMeshWall::leafBox
(
    const subVolumeProperties& subVolume
)
{
    // uniform lattice: edge charCellSize/levelOfDivision in every
    // direction, centred on the sub-volume centre
    const scalar svEdge
    (
        bbMatrix_.getCharCellSize()
       /virtualMeshLevel::getLevelOfDivision()
    );

    return boundBox
    (
        subVolume.center - 0.5*svEdge*vector::one,
        subVolume.center + 0.5*svEdge*vector::one
    );
}
//---------------------------------------------------------------------------//
scalar virtualMeshWall::evaluateContactAreaExact
(
    const planePolyClip::halfSpace& wallPlane
)
{
    // exact wetted area on the wall plane, from the inside leaves
    // of the preceding evaluateContact flood (persisted in
    // insideLeaves_): per band leaf, the polygon where the kept
    // material (body side of the surface plane, behind the wall)
    // meets the wall plane is the wetted footprint of that leaf
    if (    !virtualMeshLevel::getExactSubVolume()
        ||  insideLeaves_.size() == 0)
    {
        return -1;
    }

    scalar exactArea(0);
    label nBandNoPolygon(0);

    const scalar svEdge
    (
        bbMatrix_.getCharCellSize()
       /virtualMeshLevel::getLevelOfDivision()
    );
    const scalar leafFaceArea(svEdge*svEdge);

    forAll(insideLeaves_,lI)
    {
        subVolumeProperties& cSubVolume = bbMatrix_[insideLeaves_[lI]];
        vector svI(insideLeaves_[lI]);

        if (!isOnContactBand(svI))
        {
            // interior leaf: full footprint in the wall plane
            exactArea += leafFaceArea;
            continue;
        }

        const boundBox leaf(leafBox(cSubVolume));

        point bP(vector::zero);
        vector bN(vector::zero);

        if (cGeomModel_.getBoxSurfacePlane(leaf, bP, bN))
        {
            // kept material = body side of the surface plane and
            // behind the wall; the hs1 closing face is the wetted
            // footprint of this leaf in the wall plane
            scalar V(0);
            vector C(vector::zero);
            DynamicList<point> face;

            planePolyClip::volumeCentroidAndFace
            (
                leaf,
                wallPlane,
                true,
                planePolyClip::halfSpace(bP, bN),
                V, C, face
            );

            if (face.size() >= 3)
            {
                exactArea += planePolyClip::polygonArea(face);
            }
            else if (V > VSMALL)
            {
                // kept region fully behind the wall plane (the
                // plane does not cut the leaf): full face area
                exactArea += leafFaceArea;
            }
            else
            {
                // no kept material behind the wall in this leaf
                nBandNoPolygon++;
            }
        }
        else
        {
            // no reliable body plane in this leaf: full face area
            exactArea += leafFaceArea;
        }
    }

    if (exactArea < VSMALL)
    {
        // detected-but-clipped-to-zero: same rationale as the
        // volume pass -- the flood found inside leaves but no
        // leaf wets the wall plane, which a plane estimate failed
        // rather than a genuinely dry patch. keep the legacy
        // count-based area
        WarningInFunction
            << "virtualMeshWall::evaluateContactAreaExact: the flood "
            << "found " << insideLeaves_.size() << " inside sub-volumes "
            << "but the band footprints summed to zero ("
            << nBandNoPolygon
            << " band leaves without a wall-plane polygon) -- "
            << "keeping the legacy count-based area. "
            << "Virtual mesh bBox: " << bbMatrix_.getBBox()
            << ", matrixSize: " << bbMatrix_.getMatrixSize()
            << ", startingPoint: " << vMeshWallInfo_.getStartingPoint()
            << endl;

        return -1;
    }

    return exactArea;
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
    subVolume.resetClip();
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
