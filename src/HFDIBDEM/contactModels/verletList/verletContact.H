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

#include <list>
#include <memory>

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

    std::list<std::pair<std::shared_ptr<verletBox>, std::shared_ptr<verletBox>>> contactBoxes_;

    public:

    // Constructors ------------------------------------------------
    verletContact():
    cBodyID_(-1),
    tBodyID_(-1)
    {};

    verletContact
    (
        std::pair<label, label> cPair
    )
    :
    cBodyID_(cPair.first),
    tBodyID_(cPair.second)
    {};

    // Destructors -------------------------------------------------
    ~verletContact(){};

    void addContactBox
    (
        std::shared_ptr<verletBox> cBox,
        std::shared_ptr<verletBox> tBox
    )
    {
        // Check if the contact point is already in the list
        for (std::pair<std::shared_ptr<verletBox>, std::shared_ptr<verletBox>>& cB : contactBoxes_)
        {
            if (cB.first == cBox && cB.second == tBox)
            {
                return;
            }
        }

        contactBoxes_.push_back(std::make_pair(cBox, tBox));
    }

    void removeContactBox
    (
        std::shared_ptr<verletBox> cBox,
        std::shared_ptr<verletBox> tBox
    )
    {
        contactBoxes_.remove_if(
            [&cBox, &tBox] (const std::pair<std::shared_ptr<verletBox>, std::shared_ptr<verletBox>>& cB)
            { return cB.first == cBox && cB.second == tBox; });
    }

    std::list<std::pair<std::shared_ptr<verletBox>, std::shared_ptr<verletBox>>>& getContactBoxes()
    {
        return contactBoxes_;
    }

}; // End class verletContact

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
