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
    class for Verlet contact

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef verletContact_H
#define verletContact_H

#include "verletPoint.H"

#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class verletContact Declaration
\*---------------------------------------------------------------------------*/

class verletContact
{
    private:

    label cBodyID_;
    label tBodyID_;

    std::vector<Tuple2<verletPoint*, verletPoint*>> contactPoints_;

    public:

    // Constructors ------------------------------------------------
    verletContact
    (
        Tuple2<label, label> cPair
    )
    :
    cBodyID_(cPair.first()),
    tBodyID_(cPair.second())
    {};

    // Destructors -------------------------------------------------
    ~verletContact(){};

    void addContactPoint
    (
        verletPoint* cPoint,
        verletPoint* tPoint
    )
    {
        // Check if the contact point is already in the list
        for (Tuple2<verletPoint*, verletPoint*>& cP : contactPoints_)
        {
            if
            (
                cP.first() == cPoint
             && cP.second() == tPoint
            )
            {
                return;
            }
        }

        contactPoints_.push_back(Tuple2<verletPoint*, verletPoint*>(cPoint, tPoint));
    }

    void removeContactPoint
    (
        verletPoint* cPoint,
        verletPoint* tPoint
    )
    {
        for (auto it = contactPoints_.begin(); it != contactPoints_.end(); ++it)
        {
            if
            (
                it->first() == cPoint
             && it->second() == tPoint
            )
            {
                contactPoints_.erase(it);
                return;
            }
        }
    }

    std::vector<Tuple2<verletPoint*, verletPoint*>>& getContactPoints()
    {
        return contactPoints_;
    }

}; // End class verletContact

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
