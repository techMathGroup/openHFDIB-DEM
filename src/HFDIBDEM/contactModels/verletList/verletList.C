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

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/
#include "verletList.H"

using namespace Foam;

//---------------------------------------------------------------------------//
verletList::verletList()
:
verletLists_(3),
cntNeighList_(3)
{}
verletList::~verletList()
{}
//---------------------------------------------------------------------------//
void verletList::swapVerletPoints(
    std::shared_ptr<verletPoint> a,
    std::shared_ptr<verletPoint> b,
    label coord
)
{
    if(a->getBodyId() != b->getBodyId())
    {
        if(!a->isMin() && b->isMin())
        {
            cPair newPair(
                min(a->getBodyId(), b->getBodyId()),
                max(a->getBodyId(), b->getBodyId())
            );
            addCPairToCntNList(newPair, coord, a, b);
        }
        else if(a->isMin() && !b->isMin())
        {
            cPair cPair(
                min(a->getBodyId(), b->getBodyId()),
                max(a->getBodyId(), b->getBodyId())
            );

            if(cntNeighList_[coord].find(cPair) != cntNeighList_[coord].end())
            {
                if (a->getBodyId() < b->getBodyId())
                {
                    cntNeighList_[coord][cPair].removeContactBox(a->getParentBox(), b->getParentBox());
                }
                else
                {
                    cntNeighList_[coord][cPair].removeContactBox(b->getParentBox(), a->getParentBox());
                }

                if (cntNeighList_[coord][cPair].getContactBoxes().size() == 0)
                {
                    cntNeighList_[coord].erase(cPair);
                }
            }

            if(posCntList_.find(cPair) != posCntList_.end() &&
                ((cntNeighList_[0].find(cPair) == cntNeighList_[0].end() && emptyDim != 0) ||
                (cntNeighList_[1].find(cPair) == cntNeighList_[1].end() && emptyDim != 1) ||
                (cntNeighList_[2].find(cPair) == cntNeighList_[2].end() && emptyDim != 2)))
            {
                posCntList_.erase(cPair);
            }
        }
    }
}
//---------------------------------------------------------------------------//
void verletList::addCPairToCntNList(
    cPair cPair,
    label coord,
    std::shared_ptr<verletPoint> a,
    std::shared_ptr<verletPoint> b
)
{
    if(cntNeighList_[coord].find(cPair) == cntNeighList_[coord].end())
    {
        verletContact newContact(cPair);
        if (a->getBodyId() < b->getBodyId())
        {
            newContact.addContactBox(a->getParentBox(), b->getParentBox());
        }
        else
        {
            newContact.addContactBox(b->getParentBox(), a->getParentBox());
        }

        cntNeighList_[coord].insert({cPair, newContact});
    }
    else
    {
        if (a->getBodyId() < b->getBodyId())
        {
            cntNeighList_[coord][cPair].addContactBox(a->getParentBox(), b->getParentBox());
        }
        else
        {
            cntNeighList_[coord][cPair].addContactBox(b->getParentBox(), a->getParentBox());
        }
    }

    if(a->getIsStatic() && b->getIsStatic())
    {
        return;
    }
    if(posCntList_.find(cPair) == posCntList_.end() &&
        (cntNeighList_[0].find(cPair) != cntNeighList_[0].end() || emptyDim == 0) &&
        (cntNeighList_[1].find(cPair) != cntNeighList_[1].end() || emptyDim == 1) &&
        (cntNeighList_[2].find(cPair) != cntNeighList_[2].end() || emptyDim == 2))
    {
        posCntList_.insert(cPair);
    }
}
//---------------------------------------------------------------------------//
void verletList::addBodyToVList(immersedBody& ib)
{
    List<std::shared_ptr<boundBox>> bBoxes = ib.getGeomModel().getBBoxes();

    forAll(bBoxes, bBox)
    {
        verletBoxes_.push_back(verletBox::create(
            ib.getBodyId(),
            bBoxes[bBox],
            ib.getbodyOperation()==0
        ));

        verletBoxes_.back()->setVerletPoints();

        forAll(verletLists_, vListI)
        {
            if(vListI != emptyDim)
            {
                verletLists_[vListI].push_back(
                    verletBoxes_.back()->getMinPoint()
                );

                verletLists_[vListI].push_back(
                    verletBoxes_.back()->getMaxPoint()
                );
            }
        }
    }
}
//---------------------------------------------------------------------------//
void verletList::removeBodyFromVList(immersedBody& ib)
{
    forAll (verletLists_, coordI)
    {
        verletLists_[coordI].remove_if(
            [&ib](std::shared_ptr<verletPoint>& vPoint)
            {
                return vPoint->getBodyId() == ib.getBodyId();
            }
        );

        for (auto it = cntNeighList_[coordI].begin();
            it != cntNeighList_[coordI].end();)
        {
            if (it->first.first == ib.getBodyId()
                ||
                it->first.second == ib.getBodyId())
            {
                it = cntNeighList_[coordI].erase(it);
            }
            else
            {
                ++it;
            }
        }
    }

    verletBoxes_.remove_if(
        [&ib](std::shared_ptr<verletBox>& vBox)
        {
            return vBox->getBodyId() == ib.getBodyId();
        }
    );

    for (auto iter = posCntList_.begin(); iter != posCntList_.end();)
    {
        if (iter->first == ib.getBodyId()
            ||
            iter->second == ib.getBodyId())
        {
            iter = posCntList_.erase(iter);
        }
        else
        {
            ++iter;
        }
    }
}
//---------------------------------------------------------------------------//
void verletList::initialSorting()
{
    for (label coord = 0; coord < 3; ++coord)
    {
        verletLists_[coord].sort([&coord](
            std::shared_ptr<verletPoint>& vPoint1,
            std::shared_ptr<verletPoint>& vPoint2)
        {
            return vPoint1->getPoint()[coord] < vPoint2->getPoint()[coord];
        });

        std::list<std::shared_ptr<verletPoint>> openedIb;
        for (auto iter = verletLists_[coord].begin();
            iter != verletLists_[coord].end(); ++iter)
        {
            if((*iter)->isMin())
            {
                label curIb = (*iter)->getBodyId();
                for (auto oIbIter = openedIb.begin();
                    oIbIter != openedIb.end(); ++oIbIter)
                {
                    cPair newPair(
                        min(curIb, (*oIbIter)->getBodyId()), max(curIb, (*oIbIter)->getBodyId())
                    );
                    addCPairToCntNList(newPair, coord, *iter, *oIbIter);
                }

                openedIb.push_back(*iter);
            }
            else
            {
                label curIb = (*iter)->getBodyId();
                openedIb.remove_if(
                    [&curIb](std::shared_ptr<verletPoint>& vPoint)
                    {
                        return vPoint->getBodyId() == curIb;
                    }
                );
            }
        }
    }
}
//---------------------------------------------------------------------------//
void verletList::update(PtrList<immersedBody>& ibs)
{
    forAll(ibs, ibi)
    {
        ibs[ibi].getGeomModel().getBBoxes();
    }

    // use insertion sort to sort the neighbour list
    // besides the first sorting, this is O(N) efficient
    for (label coord = 0; coord < 3; ++coord)
    {
        if (!verletLists_[coord].empty())
        {
            auto it2 = verletLists_[coord].begin();
            auto it1 = it2++;

            while(true)
            {
                if((*it1)->getPoint()[coord]
                    > (*it2)->getPoint()[coord])
                {
                    swapVerletPoints(*it1, *it2, coord);
                    std::swap(*it1, *it2);

                    if (it1 != verletLists_[coord].begin())
                    {
                        it2 = it1--;
                        continue;
                    }
                }
                it1 = it2++;
                if (it2 == verletLists_[coord].end())
                {
                    break;
                }
            }
        }
    }
}
//---------------------------------------------------------------------------//
void verletList::computePotentialBodies
(
    PtrList<immersedBody>& ibs,
    const HashTable<scalar,label,Hash<label>>& sweepDist,
    const HashSet<label,Hash<label>>& hardFlags,
    HashSet<label,Hash<label>>& potentialBodies
)
{
    potentialBodies.clear();

    // per-coordinate overlap of sweep-inflated intervals; a pair is a
    // potential contact if it overlaps in every non-empty coordinate
    cPairHasSet overlapCnt[3];

    for (label coord = 0; coord < 3; ++coord)
    {
        if (coord == emptyDim) continue;

        // build displaced interval endpoints (copies - the live Verlet
        // points alias the shared bboxes and must not be touched)
        struct sweptPoint
        {
            label bodyId;
            scalar pos;
            bool isMin;
        };
        std::vector<sweptPoint> swept;
        swept.reserve(2*verletBoxes_.size());

        for (const auto& vBox : verletBoxes_)
        {
            const label bId(vBox->getBodyId());

            if (!ibs[bId].getIsActive()) continue;

            scalar s(0);
            if (hardFlags.found(bId))
            {
                s = GREAT;
            }
            else if (sweepDist.found(bId))
            {
                s = sweepDist[bId];
            }

            swept.push_back({bId, vBox->getBBox()->min()[coord] - s, true});
            swept.push_back({bId, vBox->getBBox()->max()[coord] + s, false});
        }

        std::sort(swept.begin(), swept.end(),
            [](const sweptPoint& a, const sweptPoint& b)
            {
                return a.pos < b.pos;
            });

        // sweep: an opened (min not yet closed) body overlaps every later
        // min before its own max is reached
        std::unordered_set<label> opened;
        for (const sweptPoint& p : swept)
        {
            if (p.isMin)
            {
                for (const label oId : opened)
                {
                    if (oId == p.bodyId) continue;

                    overlapCnt[coord].insert(
                        cPair(min(oId, p.bodyId), max(oId, p.bodyId)));
                }
                opened.insert(p.bodyId);
            }
            else
            {
                opened.erase(p.bodyId);
            }
        }
    }

    // intersect the coordinate sets. NOTE: every non-empty dimension
    // participates - an empty overlap set in an active coordinate means
    // no pair overlaps in that coordinate and empties the result. Only
    // the empty (2D) dimension is excluded. The first active coordinate
    // seeds the candidates (unconditionally, even if it is empty).
    cPairHasSet candidates;
    bool seeded(false);
    for (label coord = 0; coord < 3; ++coord)
    {
        if (coord == emptyDim) continue;

        if (!seeded)
        {
            candidates = overlapCnt[coord];
            seeded = true;
        }
        else
        {
            cPairHasSet intersected;
            for (const cPair& pair : overlapCnt[coord])
            {
                if (candidates.find(pair) != candidates.end())
                {
                    intersected.insert(pair);
                }
            }
            candidates = std::move(intersected);
        }
    }

    for (const cPair& pair : candidates)
    {
        if (!ibs[pair.first].getIsActive() || !ibs[pair.second].getIsActive())
        {
            continue;
        }

        bool cStatic(ibs[pair.first].getbodyOperation() == 0);
        bool tStatic(ibs[pair.second].getbodyOperation() == 0);

        if (cStatic && tStatic) continue;

        // report both pair members (static ones included - the
        // caller filters static bodies out of the movement sets)
        potentialBodies.insert(pair.first);
        potentialBodies.insert(pair.second);
    }
}
//---------------------------------------------------------------------------//
void verletList::computeInflatedPairs
(
    PtrList<immersedBody>& ibs,
    const HashTable<scalar,label,Hash<label>>& sweepDist,
    const HashSet<label,Hash<label>>& hardFlags,
    cPairHasSet& inflatedPairs
)
{
    inflatedPairs.clear();

    // same sweep-and-prune as computePotentialBodies, but the output is
    // the pair list itself and hard-flagged bodies are excluded from
    // the sweep: their GREAT-inflated intervals would pair them with
    // every body and pin the DEM time step
    cPairHasSet overlapCnt[3];

    for (label coord = 0; coord < 3; ++coord)
    {
        if (coord == emptyDim) continue;

        struct sweptPoint
        {
            label bodyId;
            scalar pos;
            bool isMin;
        };
        std::vector<sweptPoint> swept;
        swept.reserve(2*verletBoxes_.size());

        for (const auto& vBox : verletBoxes_)
        {
            const label bId(vBox->getBodyId());

            if (!ibs[bId].getIsActive()) continue;
            if (hardFlags.found(bId)) continue;

            scalar s(0);
            if (sweepDist.found(bId))
            {
                s = sweepDist[bId];
            }

            swept.push_back({bId, vBox->getBBox()->min()[coord] - s, true});
            swept.push_back({bId, vBox->getBBox()->max()[coord] + s, false});
        }

        std::sort(swept.begin(), swept.end(),
            [](const sweptPoint& a, const sweptPoint& b)
            {
                return a.pos < b.pos;
            });

        std::unordered_set<label> opened;
        for (const sweptPoint& p : swept)
        {
            if (p.isMin)
            {
                for (const label oId : opened)
                {
                    if (oId == p.bodyId) continue;

                    overlapCnt[coord].insert(
                        cPair(min(oId, p.bodyId), max(oId, p.bodyId)));
                }
                opened.insert(p.bodyId);
            }
            else
            {
                opened.erase(p.bodyId);
            }
        }
    }

    // intersect the coordinate sets (same logic as
    // computePotentialBodies)
    cPairHasSet candidates;
    bool seeded(false);
    for (label coord = 0; coord < 3; ++coord)
    {
        if (coord == emptyDim) continue;

        if (!seeded)
        {
            candidates = overlapCnt[coord];
            seeded = true;
        }
        else
        {
            cPairHasSet intersected;
            for (const cPair& pair : overlapCnt[coord])
            {
                if (candidates.find(pair) != candidates.end())
                {
                    intersected.insert(pair);
                }
            }
            candidates = std::move(intersected);
        }
    }

    for (const cPair& pair : candidates)
    {
        if (!ibs[pair.first].getIsActive() || !ibs[pair.second].getIsActive())
        {
            continue;
        }

        bool cStatic(ibs[pair.first].getbodyOperation() == 0);
        bool tStatic(ibs[pair.second].getbodyOperation() == 0);

        if (cStatic && tStatic) continue;

        inflatedPairs.insert(pair);
    }
}
//---------------------------------------------------------------------------//
