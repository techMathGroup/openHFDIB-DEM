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

Class
    Foam::facetBinning

Description
    uniform-grid broad phase over the facets of one triSurface: a
    candidate-gather replacement for the per-leaf octree walk of
    indexedOctree::findBox in the wall exact path.

    - one O(facets) build pass appends every facet to all bins its
      bounding box overlaps
    - bins are cubes of edge cbrt(surface volume/nFacets): mean
      occupancy is O(1) and a slab-sized facet only duplicates
      into the bins its bbox truly spans (target-count sizing,
      NOT max-facet-extent sizing - one slab facet must not
      collapse the grid
    - a query gathers from the bins the query box overlaps using
      the SAME floor-based bin mapping as the build, so closed-
      touching geometries (facet bbox corner exactly on a bin
      edge, query box exactly on a bin edge) share a bin and the
      gather is a superset of the findBox candidates; the
      platform treeDataPrimitivePatch::overlaps predicate (the
      one findBox itself applies, on the raw unexpanded query
      box) then prunes the gather exactly

    the returned candidate SET is therefore equal to
    indexedOctree::findBox(treeBoundBox(queryBox)) by
    construction; only the ORDER differs (findBox returns hash
    order, this returns sorted - consumers are order-independent
    since the smallest-index tie-break of getShapesSurfacePlane)

    degenerate fallback: empty surfaces, zero-extent bboxes and
    grids that collapse (total bins below nFacets/8) or overload
    (mean bucket occupancy above occupancyMax) report invalid();
    the caller falls back to tree.findBox.

    invalidation is the caller's duty: rebuild after any mutation
    of the surface. each rank owns its body surfaces and mutates
    them single-threaded between DEM substeps, so one binning
    per surface instance is safe; concurrent queries are not
    supported.

SourceFiles
    facetBinning.C

Contributors
    Martin Isoz (2019-*)
\*---------------------------------------------------------------------------*/

#ifndef facetBinning_H
#define facetBinning_H

#include "fvCFD.H"
#include "triSurface.H"
#include "treeDataTriSurface.H"
#include <algorithm>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class facetBinning Declaration
\*---------------------------------------------------------------------------*/

class facetBinning
{
    // Private Data

        //- the binned surface
        const triSurface& surf_;

        //- surface bounding box (hull of the facet bboxes)
        boundBox bBox_;

        //- cube bin edge: cbrt(volume/nFacets)
        scalar binSize_;

        //- bins per direction
        label nBins_[3];

        //- flat storage, index i + j*nBins_[0] + k*nBins_[0]*nBins_[1]
        List<DynamicList<label>> bins_;

        //- sum of the bin list sizes (the build memory footprint)
        label totalEntries_;

        //- bins covered by the worst single facet
        label maxFacetBins_;

        //- degenerate for binning: the caller must fall back to
        //- tree.findBox
        bool valid_;

public:

    //- mean bucket occupancy above which the grid is declared
    //- degenerate (measured occupancy: 3.9-6.1 across sphere,
    //- boulder, slab and particle STLs with cube-edge sizing;
    //- 16 is a loose bound against gross overload only)
    static constexpr scalar occupancyMax = 16;

    //- bin count below nFacets/8 declares the grid collapsed
    static constexpr label occupancyBinsPerFacet = 8;

    // Constructors

        //- bin surf; check valid() before querying
        facetBinning(const triSurface& surf);

    // Member Functions

        //- degenerate (empty, zero-extent, collapsed or overloaded
        //- grid): fall back to tree.findBox
        bool valid() const
        {
            return valid_;
        }

        //- sorted candidate set of facets overlapping queryBox.
        //  shapes is the octree's own tree().shapes(), so the
        //  filter is the platform predicate with its tolerance;
        //  set-equal to indexedOctree::findBox(treeBoundBox(queryBox))
        labelList facetsNear
        (
            const boundBox& queryBox,
            const treeDataTriSurface& shapes
        ) const;

        //- number of binned facets
        label nFacets() const
        {
            return surf_.size();
        }

        //- total number of bins
        label nBins() const
        {
            return nBins_[0]*nBins_[1]*nBins_[2];
        }

        //- the cube bin edge
        scalar binSize() const
        {
            return binSize_;
        }

        //- the binned extent
        const boundBox& bounds() const
        {
            return bBox_;
        }

        //- mean bucket occupancy (totalEntries/nBins)
        scalar occupancy() const;

        //- bins covered by the worst single facet
        label maxFacetBins() const
        {
            return maxFacetBins_;
        }

        //- upper bound on the per-facet bin coverage (surface
        //  extent over bin size, per component): the memory
        //  bound of the build
        label facetBinsBound() const;

    private:

        //- flat bin index of a position, clamped to the grid.
        //  shared by the build (facet bbox corners) and the
        //  query (query box corners) so closed-touching
        //  geometries share a bin
        label binIndex(const point& p, const label dir) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
