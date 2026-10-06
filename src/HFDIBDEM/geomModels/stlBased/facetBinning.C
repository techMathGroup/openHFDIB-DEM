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

Contributors
    Martin Isoz (2019-*)
\*---------------------------------------------------------------------------*/

#include "facetBinning.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::facetBinning::facetBinning(const triSurface& surf)
:
surf_(surf),
bBox_(boundBox::invertedBox),
binSize_(0),
nBins_{1, 1, 1},
bins_(),
totalEntries_(0),
maxFacetBins_(0),
valid_(false)
{
    const label nFacets(surf_.size());

    if (nFacets == 0)
    {
        return;
    }

    // hull of the facet bboxes; every facet lands in >= 1 bin
    forAll(surf_, fI)
    {
        bBox_.add(surf_.points(), surf_[fI]);
    }

    const vector span(bBox_.span());

    // target-count sizing: cube edge so the mean occupancy is
    // O(1) given nFacets and the surface volume. cbrt(volume)
    // also shrinks with a thin direction, so a plate or slab
    // surface does not explode the bin count;
    binSize_ = cbrt(span[0]*span[1]*span[2])/cbrt(scalar(nFacets));

    if (binSize_ < VSMALL)
    {
        // zero-volume surface (fully flat): no binning
        return;
    }

    for (label dir = 0; dir < 3; dir++)
    {
        // ceil: a thin direction still gets one bin layer
        nBins_[dir] = max(label(1), label(ceil(span[dir]/binSize_)));
    }

    // the bin count must scale with the facet count; below
    // nFacets/8 the grid has collapsed and every query gathers
    // a large fraction of the surface - not worth the memory
    if (nBins() < nFacets/occupancyBinsPerFacet)
    {
        return;
    }

    bins_.setSize(nBins());

    // one O(facets) pass: append every facet to all bins its
    // bbox overlaps
    forAll(surf_, fI)
    {
        const triFace& f(surf_[fI]);
        const pointField& sPts(surf_.points());

        const point fMin
        (
            min(min(sPts[f[0]], sPts[f[1]]), sPts[f[2]])
        );
        const point fMax
        (
            max(max(sPts[f[0]], sPts[f[1]]), sPts[f[2]])
        );

        const label iMin(binIndex(fMin, 0));
        const label iMax(binIndex(fMax, 0));
        const label jMin(binIndex(fMin, 1));
        const label jMax(binIndex(fMax, 1));
        const label kMin(binIndex(fMin, 2));
        const label kMax(binIndex(fMax, 2));

        maxFacetBins_ = max
        (
            maxFacetBins_,
            (iMax - iMin + 1)*(jMax - jMin + 1)*(kMax - kMin + 1)
        );

        for (label k = kMin; k <= kMax; k++)
        {
            for (label j = jMin; j <= jMax; j++)
            {
                for (label i = iMin; i <= iMax; i++)
                {
                    bins_[i + j*nBins_[0] + k*nBins_[0]*nBins_[1]].append(fI);
                    totalEntries_++;
                }
            }
        }
    }

    // mean occupancy above the bound: the grid overloads (the
    // extreme facets dominate every bucket) - not worth the
    // memory, the caller falls back to tree.findBox
    if (occupancy() > occupancyMax)
    {
        bins_.clear();
        return;
    }

    valid_ = true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::facetBinning::binIndex(const point& p, const label dir) const
{
    // floor of the position in bin units, clamped to the grid.
    // the build (facet bbox corners) and the query (query box
    // corners) share this mapping, so closed-touching
    // geometries share a bin and the gather is a superset of
    // the true bbox overlaps (a shared coordinate maps to a
    // shared bin by monotonicity of floor: floor(max(a,b))
    // equals max(floor(a),floor(b)) exactly in IEEE 754)
    label i(floor((p[dir] - bBox_.min()[dir])/binSize_));

    return min(max(i, label(0)), nBins_[dir] - 1);
}

Foam::scalar Foam::facetBinning::occupancy() const
{
    if (nBins() == 0)
    {
        return GREAT;
    }
    return scalar(totalEntries_)/scalar(nBins());
}

Foam::label Foam::facetBinning::facetBinsBound() const
{
    // ceil of the surface extent in bins, per direction: the
    // memory bound of the build (the worst facet cannot cover
    // more bins than the whole grid)
    const vector span(bBox_.span());

    label bound(1);
    for (label dir = 0; dir < 3; dir++)
    {
        bound *= max(label(1), label(ceil(span[dir]/binSize_)));
    }
    return bound;
}

// * * * * * * * * binning query * * * * * * * * * * * * * * * * * * * * * * //

Foam::labelList Foam::facetBinning::facetsNear
(
    const boundBox& queryBox,
    const treeDataTriSurface& shapes
) const
{
    if (!valid_)
    {
        return labelList(0);
    }

    DynamicList<label> candidates;

    const label iMin(binIndex(queryBox.min(), 0));
    const label iMax(binIndex(queryBox.max(), 0));
    const label jMin(binIndex(queryBox.min(), 1));
    const label jMax(binIndex(queryBox.max(), 1));
    const label kMin(binIndex(queryBox.min(), 2));
    const label kMax(binIndex(queryBox.max(), 2));

    for (label k = kMin; k <= kMax; k++)
    {
        for (label j = jMin; j <= jMax; j++)
        {
            for (label i = iMin; i <= iMax; i++)
            {
                const DynamicList<label>& bin
                (
                    bins_[i + j*nBins_[0] + k*nBins_[0]*nBins_[1]]
                );
                forAll(bin, e)
                {
                    candidates.append(bin[e]);
                }
            }
        }
    }

    // deduplicate (a facet spanning bins appears once per
    // bin), then filter exactly with the platform predicate
    // findBox itself applies - on the raw unexpanded query box.
    // the gather is a superset, so the surviving set is
    // exactly the findBox candidate set
    sort(candidates);
    candidates.setSize
    (
        std::unique(candidates.begin(), candidates.end())
      - candidates.begin()
    );

    const treeBoundBox searchBox(queryBox);
    label nKeep(0);
    forAll(candidates, i)
    {
        if (shapes.overlaps(candidates[i], searchBox))
        {
            candidates[nKeep++] = candidates[i];
        }
    }
    candidates.setSize(nKeep);

    return labelList(candidates);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
