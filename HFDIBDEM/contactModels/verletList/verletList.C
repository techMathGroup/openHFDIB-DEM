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

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
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
void verletList::mergeSortList(IDLList<verletPoint>& listToSort, label coord)
{
    if(listToSort.size() > 1)
    {
        label half = floor(listToSort.size() / 2);
        IDLList<verletPoint> secondList;

        for (label i = 0; i < half; ++i)
        {
            secondList.append(listToSort.removeHead());
        }

        mergeSortList(listToSort, coord);
        mergeSortList(secondList, coord);

        IDLList<verletPoint> resultList;

        while(listToSort.size() > 0 && secondList.size() > 0)
        {
            if(listToSort.first()->getPoint()[coord] 
                < secondList.first()->getPoint()[coord])
            {
                resultList.append(listToSort.removeHead());
            }
            else
            {
                resultList.append(secondList.removeHead());
            }
        }

        while(listToSort.size() > 0)
        {
            resultList.append(listToSort.removeHead());
        }
        while(secondList.size() > 0)
        {
            resultList.append(secondList.removeHead());
        }
        listToSort = resultList;
    }
}
//---------------------------------------------------------------------------//
void verletList::swapVerletPoints(verletPoint* a, verletPoint* b, label coord)
{
    if(a->getBodyID() != b->getBodyID())
    {
        if(!a->isMin() && b->isMin())
        {
            Tuple2<label, label> newPair(
                min(a->getBodyID(), b->getBodyID()), 
                max(a->getBodyID(), b->getBodyID())
            );
            addCPairToCntNList(newPair, coord);
        }
        else if(a->isMin() && !b->isMin())
        {
            Tuple2<label, label> cPair(
                min(a->getBodyID(), b->getBodyID()), 
                max(a->getBodyID(), b->getBodyID())
            );
            if(cntNeighList_[coord].found(cPair))
            {
                cntNeighList_[coord].erase(cPair);
            }
            if(posCntList_.found(cPair))
            {
                posCntList_.erase(cPair);
            }
        }
    }

    verletLists_[coord].swapDown(a);
}
//---------------------------------------------------------------------------//
void verletList::addCPairToCntNList(Tuple2<label, label> cPair, label coord)
{
    if(!cntNeighList_[coord].found(cPair))
    {
        cntNeighList_[coord].insert(cPair);
    }

    if(!posCntList_.found(cPair) && 
        (cntNeighList_[0].found(cPair) || emptyDim == 0) &&
        (cntNeighList_[1].found(cPair) || emptyDim == 1) &&
        (cntNeighList_[2].found(cPair) || emptyDim == 2))
    {
        posCntList_.insert(cPair);
    }
}
//---------------------------------------------------------------------------//
void verletList::addBodyToVList(immersedBody& ib)
{
    List<boundBox*> bBoxes = ib.getGeomModel().getBBoxes();
    PtrList<verletPoint> newVerletPoints;
    forAll(bBoxes, bBox)
    {
        forAll(verletLists_, vListI)
        {
            if(vListI != emptyDim)
            {
                verletLists_[vListI].append(new verletPoint(
                    ib.getBodyId(),
                    bBoxes[bBox]->min(),
                    true
                ));

                verletLists_[vListI].append(new verletPoint(
                    ib.getBodyId(),
                    bBoxes[bBox]->max(),
                    false
                ));
            }
        }
    }
}
//---------------------------------------------------------------------------//
void verletList::initialSorting()
{
    for (label coord = 0; coord < 3; ++coord)
    {
        mergeSortList(verletLists_[coord], coord);

        labelHashSet openedIb;
        for (auto i = verletLists_[coord].begin(); 
            i != verletLists_[coord].end(); ++i)
        {
            if((*i).isMin())
            {
                label curIb = (*i).getBodyID();
                for (auto oIbIter = openedIb.begin();
                    oIbIter != openedIb.end(); ++oIbIter)
                {
                    // Create a Tuple2 for bodies that are newly 
                    // intersected always keep first the body 
                    // with lower bodyID
                    Tuple2<label, label> newPair(
                        min(curIb, oIbIter.key()), max(curIb, oIbIter.key())
                    );
                    addCPairToCntNList(newPair, coord);
                }

                openedIb.insert(curIb);
            }
            else
            {
                label curIb = (*i).getBodyID();
                if(openedIb.found(curIb))
                {
                    openedIb.erase(curIb);
                }
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
            IDLList<verletPoint>::link* it = verletLists_[coord].first();
            while(it != verletLists_[coord].last())
            {
                if(it != verletLists_[coord].first())
                {
                    verletPoint* currPoint = 
                        static_cast<verletPoint*>(it);

                    verletPoint* prevPoint = 
                        static_cast<verletPoint*>(it->prev_);

                    if(prevPoint->getPoint()[coord]
                        > currPoint->getPoint()[coord])
                    {
                        swapVerletPoints(prevPoint, currPoint, coord);
                        it = it->prev_;
                    }
                }

                it = it->next_;
            }
        }
    }
}
//---------------------------------------------------------------------------//
