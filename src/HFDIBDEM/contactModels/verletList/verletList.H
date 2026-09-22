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
    class for Verlet lists to optimize contact detection

SourceFiles
    verletList.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef verletList_H
#define verletList_H

#include "dictionary.H"
#include "fvCFD.H"
#include "verletPoint.H"
#include "verletContact.H"
#include "immersedBody.H"

#include <list>
#include <map>
#include <unordered_set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct hashFunction
{
  size_t operator()(const std::pair<label, label> &x) const
  {
    return x.first ^ x.second;
  }
};

/*---------------------------------------------------------------------------*\
                        Class verletList Declaration
\*---------------------------------------------------------------------------*/

class verletList
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

    typedef std::pair<label, label> cPair;
    typedef std::map<cPair, verletContact> cPairVcontactMap;

    typedef std::unordered_set<cPair, hashFunction> cPairHasSet;

    // List of all verlet boxes
    std::list<std::shared_ptr<verletBox>> verletBoxes_;

    List<std::list<std::shared_ptr<verletPoint>>> verletLists_;
    // List for pairs of intersected Bodies in specific coordinate
    List<cPairVcontactMap> cntNeighList_;
    // List of possible contact du to intersection of bounding boxes.
    // Pair is added to this list when it is present in all 3 cntNeighList_.
    cPairHasSet posCntList_;

    void swapVerletPoints(
        std::shared_ptr<verletPoint> a,
        std::shared_ptr<verletPoint> b,
        label coord
    );

    void addCPairToCntNList(
        cPair cPair,
        label coord,
        std::shared_ptr<verletPoint> a,
        std::shared_ptr<verletPoint> b
    );

    public:

    // Constructors ------------------------------------------------
    verletList();

    // Destructors -------------------------------------------------
    ~verletList();

    // Add new body to be tracked by Verlet list
    void addBodyToVList(immersedBody& ib);

    // Remove body from Verlet list
    void removeBodyFromVList(immersedBody& ib);

    // More efficient initial soting
    void initialSorting();

    // Update positions of bounding points
    // sort the list
    // solve potential intersections
    void update(PtrList<immersedBody>& ibs);

    // Conservative contact-potential test (adaptiveStepDEM): fills
    // potentialBodies with the bodies (both pair members, static ones
    // included) of every pair whose bounding boxes, each inflated by the
    // per-body sweep distance sweepDist[i], overlap in every non-empty
    // coordinate. Callers use it to sort bodies into the sub-cycled /
    // single-stepped sets.
    // Note (MI): bodies listed in hardFlags are always sub-cycled
    // Note (MI): the test never mutates the live Verlet list
    //            (displaced intervals are built as copies)
    void computePotentialBodies
    (
        PtrList<immersedBody>& ibs,
        const HashTable<scalar,label,Hash<label>>& sweepDist,
        const HashSet<label,Hash<label>>& hardFlags,
        HashSet<label,Hash<label>>& potentialBodies
    );

    // Pair-level companion of computePotentialBodies (adaptive DEM
    // time-step sizing): fills inflatedPairs with every pair whose
    // bounding boxes, each inflated by the per-body sweep distance,
    // overlap in every non-empty coordinate. Pairs with a hard-flagged
    // member are skipped (they are inflated by GREAT and would pair
    // with every body in the case); pairs of two static bodies are
    // skipped. The callers' set classification is not affected.
    void computeInflatedPairs
    (
        PtrList<immersedBody>& ibs,
        const HashTable<scalar,label,Hash<label>>& sweepDist,
        const HashSet<label,Hash<label>>& hardFlags,
        cPairHasSet& inflatedPairs
    );

    cPairHasSet::iterator begin()
    {
        return posCntList_.begin();
    }

    cPairHasSet::iterator end()
    {
        return posCntList_.end();
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
