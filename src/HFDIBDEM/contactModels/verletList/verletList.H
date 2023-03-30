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
    class for Verlet lists to optimize contact detection

SourceFiles
    verletList.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef verletList_H
#define verletList_H

#include "dictionary.H"
#include "fvCFD.H"
#include "verletPoint.H"
#include "verletContact.H"
#include "immersedBody.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

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

    typedef HashTable<verletContact, Tuple2<label, label>,Hash<Tuple2<label, label>>>
        cPairHasTable;

    typedef HashSet<Tuple2<label, label>,Hash<Tuple2<label, label>>>
        cPairHasSet;

    List<IDLList<verletPoint>> verletLists_;
    // List for pairs of intersected Bodies in specific coordinate
    List<cPairHasTable> cntNeighList_;
    // List of possible contact du to intersection of bounding boxes.
    // Pair is added to this list when it is present in all 3 cntNeighList_.
    cPairHasSet posCntList_;

    void mergeSortList(IDLList<verletPoint>& listToSort, label coord);

    void swapVerletPoints(verletPoint* a, verletPoint* b, label coord);

    void addCPairToCntNList(
        Tuple2<label, label> cPair,
        label coord,
        verletPoint* a,
        verletPoint* b
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
