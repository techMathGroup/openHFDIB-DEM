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

Description
    Geometric resolution of contact between an immersed body of
    arbitrary shape and a planar wall, evaluated on a structured
    "virtual mesh".

    Fundamentals
    ~~~~~~~~~~~~
    - The virtual mesh is a regular Cartesian lattice of sub-volumes
      (bbMatrix) superimposed on a wall-aligned bounding box of the
      body. It is independent of the CFD mesh: its resolution is
      controlled solely by the virtualMesh dict entries (level,
      charCellSize) via the sub-volume edge
      charCellSize/levelOfDivision.

    - A sub-volume is classified with a single centroid test,
      geomModel::pointInside(); a sub-volume whose centroid lies
      inside the body marks the contact. There is no adaptive
      refinement: accuracy is set by the uniform sub-volume edge,
      which is typically finer than the CFD mesh near the surface.

    - The contact volume is obtained by a flood fill over the
      sub-volume lattice, seeded at the starting point and expanding
      through face neighbours of classified sub-volumes. The reported
      volume is (number of inside sub-volumes) x subVolumeV.

    - Two flood variants exist: corner-neighbour expansion for
      locating a first contact point (detectFirstContactPoint) and
      face-neighbour expansion for evaluating the connected contact
      volume/area (detectFirstFaceContactPoint, evaluateContact).
      Each sub-volume carries a toCheck flag; a sub-volume is queued
      at most once, so every flood terminates after at most one full
      lattice scan.

    - Wall-alignment of the box: the bounding boxes handed to this
      class have already been projected onto the wall plane (see
      wallContactInfo::findContactAreas / correctSMBBforWall), so the
      lattice itself resolves the contact patch in the plane tangent
      to the wall.

    - Pseudo-2D support: in quasi-2D setups the box is clipped to a
      single sub-volume layer in the empty direction and the flooded
      volume/area are rescaled by emptyScale =
      originalSpan/subVolumeEdge. This is exact for bodies extruded
      uniformly along the empty direction and removes the
      bodyThickness/subVolumeEdge factor from every flood.

    - Overflow safety: every flood is bounded by a visit cap
      (lattice sub-volume count, further limited by the virtualMesh
      dict entry maxSubVolumes). The cap is truncation-only: with the
      single-queue guard each flood visits each sub-volume at most
      once, so a complete lattice scan can never exceed the cap and
      the cap can never truncate a legitimate search.

SourceFiles
    virtualMeshWall.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef virtualMeshWall_H
#define virtualMeshWall_H

#include "dictionary.H"
#include "fvCFD.H"
#include "geomModel.H"
#include "bbMatrix.H"
#include "clockTime.H"
#include "subVolumeProperties.H"
// #include ""


// #include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct virtualMeshWallInfo
{
    virtualMeshWallInfo():
    bBox(boundBox()),
    startingPoint(vector::zero),
    subVolumeNVector(vector::zero),
    charCellSize(0),
    subVolumeV(0),
    emptyScale(1)

    {}
    virtualMeshWallInfo(
        boundBox bBox,
        point startingPoint,
        vector subVolumeNVector,
        scalar charCellSize,
        scalar subVolumeV,
        scalar emptyScale = 1):
    bBox(bBox),
    startingPoint(startingPoint),
    subVolumeNVector(subVolumeNVector),
    charCellSize(charCellSize),
    subVolumeV(subVolumeV),
    emptyScale(emptyScale)
    {}

    const point& getStartingPoint()
    {
        return startingPoint;
    }

    const scalar& getSVVolume()
    {
        return subVolumeV;
    }

    // Pseudo-2D support: when the bBox was clipped in the empty direction
    // to a single sub-volume layer, the flood finds only the 2D footprint
    // of the contact patch. For bodies extruded in the empty direction the
    // patch is uniform there, so multiplying the flooded volume by
    // emptyScale = originalExtent/subVolumeEdge restores the physical
    // volume (and, via volume/subVolumeV, the physical area) exactly.
    const scalar& getEmptyScale()
    {
        return emptyScale;
    }

    const scalar& getCharCellSize()
    {
        return charCellSize;
    }

    void shiftBBox
    (
        vector shiftVec
    )
    {
        bBox.min() += shiftVec;
        bBox.max() += shiftVec;
    }

    boundBox& getbBox()
    {
        return bBox;
    }

    void setStartingPoint
    (
        point sPoint
    )
    {
        startingPoint = sPoint;
    }

    boundBox bBox;
    point startingPoint;
    vector subVolumeNVector;
    scalar charCellSize;
    scalar subVolumeV;
    scalar emptyScale;
};

/*---------------------------------------------------------------------------*\
                        Class virtualMeshWall Declaration
\*---------------------------------------------------------------------------*/
class virtualMeshWall
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;
    typedef HashSet<vector,Hash<vector>> vectorHashSet;

    geomModel& cGeomModel_;

    virtualMeshWallInfo& vMeshWallInfo_;

    bbMatrix bbMatrix_;

    DynamicPointList edgeSubVolumesPoints_;

    // Append a neighbour sub-volume to the next flood-fill layer only if it
    // has not been visited yet (toCheck == true). Without this guard the
    // floods re-queue already-processed sub-volumes and their frontiers grow
    // without bound on large virtual meshes.
    void checkAndAppend
    (
        const vector& svI,
        DynamicVectorList& auxToCheck
    );

    // Same as checkAndAppend but using face neighbours
    void checkAndAppendFace
    (
        const vector& svI,
        DynamicVectorList& auxToCheck
    );

    public:
    
    vector contactCenter_ = vector::zero;

    // Constructors ------------------------------------------------

    virtualMeshWall
    (
        virtualMeshWallInfo& vMeshWallInfo,
        geomModel& cGeomModel
    );


    // Destructors -------------------------------------------------
    ~virtualMeshWall();

    bool detectFirstContactPoint();

    bool detectFirstFaceContactPoint();

    scalar evaluateContact();

    void checkSubVolume(subVolumeProperties& subVolume);

    void resetSubVolume(subVolumeProperties& subVolume);

    vector shiftEdgeSubVolume
    (
        vector& subVolumeIndex
    );

    label getInternalSV();

    DynamicPointList& getEdgeSVPoints()
    {
        return edgeSubVolumesPoints_;
    }

    vector& getContactCenter()
    {
        return contactCenter_;
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
