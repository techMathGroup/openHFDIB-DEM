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

Description
    Geometric resolution of contact between two immersed bodies of
    arbitrary shape, evaluated on an adaptive octree "virtual mesh".

    Fundamentals
    ~~~~~~~~~~~~
    - The virtual mesh is an axis-aligned octree superimposed on the
      intersection of the two bodies' bounding boxes. It is independent
      of the CFD mesh: its refinement is controlled solely by the
      virtualMesh dict entries (level, charCellSize) and by the geometry
      of the contacting pair.

    - Each octree node (subVolume) carries, for each of the two bodies
      separately, a volumeType classification -- INSIDE, OUTSIDE or
      MIXED -- obtained through geomModel::getVolumeType(). A node
      classified as MIXED (i.e. cut by the body surface) is refined
      into its eight octants; refinement stops once the node volume
      drops below the sub-volume resolution subVolumeV.

    - The contact volume is assembled as the sum of the volumes of all
      leaf sub-volumes classified INSIDE with respect to both bodies
      simultaneously. At leaf level, the axis-aligned bounding boxes
      obtained from geomModel::limitFinalSubVolume() replace the raw
      octant boxes, so the leaf estimate follows the local surface
      orientation rather than the octree frame.

    - The contact normal and area are derived from the set of "edge"
      sub-volumes (leaves that remain MIXED for both bodies): a
      covariance-based best-fit plane through their midpoints provides
      the normal, and an angular sector sweep of the points projected
      into that plane constructs a convex hull whose area estimates the
      contact area. For non-convex bodies the edge sub-volumes are first
      clustered into connected subContacts and the estimate is carried
      out per cluster.

    - Temporal coherence: the virtualMeshInfo (root sub-volume,
      resolution, starting point, pseudo-2D emptyScale) is owned by the
      calling prtSubContactInfo and persists across time steps. The
      cached starting point (i) warm-starts the contact search and
      (ii) underpins the matching of a new sub-contact to its
      predecessor, which preserves the tangential-force history.

    - Pseudo-2D support: in quasi-2D setups the pair bounding box is
      clipped to a single sub-volume layer in the empty direction before
      the octree is built, and the resulting volume/area are rescaled by
      emptyScale = originalSpan/subVolumeEdge. This is exact for bodies
      extruded uniformly along the empty direction and reduces the
      traversal cost by the clipped factor.

    - Overflow safety: the traversal is bounded by a visit cap
      (8*rootVolume/subVolumeV, further limited by the virtualMesh
      dict entry maxSubVolumes). The cap is truncation-only: a
      complete scan of any admissible octree stays below it, so the
      cap can never truncate a legitimate search.

SourceFiles
    virtualMesh.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef virtualMesh_H
#define virtualMesh_H

#include "dictionary.H"
#include "fvCFD.H"
#include "geomModel.H"
#include "clockTime.H"
#include "subVolume.H"
#include "subContact.H"
#include "virtualMeshLevel.H"


// #include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct virtualMeshInfo
{
    virtualMeshInfo():
    sV(boundBox()),
    subVolumeV(0),
    emptyScale(1)
    {}

    virtualMeshInfo(const virtualMeshInfo& vMeshInfo):
    sV(vMeshInfo.sV),
    subVolumeV(vMeshInfo.subVolumeV),
    emptyScale(vMeshInfo.emptyScale)
    {
        if (vMeshInfo.startingPoint.valid())
        {
            startingPoint.reset(new point(vMeshInfo.startingPoint()));  //OF.com: set -> reset
        }
    }

    virtualMeshInfo(boundBox& bBox, scalar subVolumeV, scalar emptyScale = 1):
    sV(bBox),
    subVolumeV(subVolumeV),
    emptyScale(emptyScale)
    {}

    void operator= (const virtualMeshInfo& vMeshInfo)
    {
        sV = vMeshInfo.sV;
        subVolumeV = vMeshInfo.subVolumeV;
        emptyScale = vMeshInfo.emptyScale;

        if (vMeshInfo.startingPoint.valid())
        {
            startingPoint.reset(new point(vMeshInfo.startingPoint()));
        }
    }

    const point& getStartingPoint()
    {
        if (!startingPoint.valid())
        {
            startingPoint.reset(new point(sV.midpoint()));
        }

        return startingPoint();
    }

    subVolume& getSubVolume()
    {
        return sV;
    }

    // Pseudo-2D support: when the bBox was clipped in the empty direction
    // to a single sub-volume layer, the traversal finds only the 2D
    // footprint of the contact patch. For bodies extruded in the empty
    // direction the patch is uniform there, so multiplying the traversed
    // volume/area by emptyScale = originalExtent/subVolumeEdge restores
    // the physical values exactly. 1 otherwise.
    const scalar& getEmptyScale() const
    {
        return emptyScale;
    }

    subVolume sV;
    scalar subVolumeV;
    scalar emptyScale;
    autoPtr<point> startingPoint;
};

/*---------------------------------------------------------------------------*\
                        Class virtualMesh Declaration
\*---------------------------------------------------------------------------*/
class virtualMesh
{
    private:

    typedef DynamicList<label>  DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point>  DynamicPointList;
    typedef HashSet<vector,Hash<vector>> vectorHashSet;

    geomModel& cGeomModel_;
    geomModel& tGeomModel_;

    virtualMeshInfo& vMeshInfo_;

    point contactCenter_;
    DynamicPointList edgeSubVolumesPoints_;

    // Truncation-only visit cap (mirrors maxVSIter in virtualMeshWall):
    // a complete scan of this octree visits at most
    // 8*rootVolume/subVolumeV sub-volumes (leaves can be as small as
    // subVolumeV/8 and internal nodes add ~14% at most), so the cap can
    // only fire on a genuinely truncated traversal, never on a legitimate
    // one. iterMax_ is additionally limited by maxSubVolumes.
    label iterCount_;
    label iterMax_;
    bool truncated_;

    //- Warn (once per traversal) that the visit cap was reached
    void warnTruncated(const word& where);

    public:


    // Constructors ------------------------------------------------

    virtualMesh
    (
        virtualMeshInfo& vMeshInfo,
        geomModel& cGeomModel_,
        geomModel& tGeomModel_
    );


    // Destructors -------------------------------------------------
    ~virtualMesh();

    bool detectFirstContactPoint();

    // Get contact subVolume
    bool detectFirstVolumeInContact
    (
        subVolume& sV,
        bool& startPointFound
    );

    scalar evaluateContact();

    void inspectSubVolume(
        subVolume& subVolume,
        scalar& contactVolume,
        vector& contactCenter,
        DynamicPointList& edgePoints
    );

    DynamicPointList& getEdgeSVPoints()
    {
        return edgeSubVolumesPoints_;
    }

    point& getContactCenter()
    {
        vMeshInfo_.startingPoint() = contactCenter_;
        return contactCenter_;
    }

    Tuple2<scalar,vector> get3DcontactNormalAndSurface(bool nonConvex);

    Tuple2<scalar,vector> get3DcontactNormalAndSurface(DynamicPointList edgeSubVolumesPoints);

    std::vector<subContact> findsubContacts(subVolume& sV);

    bool canCombineSubContacts(
        subContact& main,
        subContact& comp
    ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
