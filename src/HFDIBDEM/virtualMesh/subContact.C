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
#include "subContact.H"

using namespace Foam;

//---------------------------------------------------------------------------//
subContact::subContact()
{}

subContact::~subContact()
{}
//---------------------------------------------------------------------------//
void subContact::addSubVolume(std::shared_ptr<subVolume> sV)
{
    if (subVolumes_.size() == 0)
    {
        boundBox_ = *sV;
    }

    subVolumes_.push_back(sV);
    tmp<pointField> points = boundBox_.points();
    points->append(sV->points());
    boundBox_ = boundBox(points,false);
    volume_ += sV->volume();
}
//---------------------------------------------------------------------------//
bool subContact::canCombine(subVolume& sV)
{
    if (!boundBox_.overlaps(sV))
    {
        return false;
    }

    for (auto& subVolume : subVolumes_)
    {
        if (subVolume->overlaps(sV))
        {
            return true;
        }
    }

    return false;
}
//---------------------------------------------------------------------------//
DynamicList<point> subContact::getEdgePoints() const
{
    DynamicList<point> edgePoints;
    for (auto& subVolume : subVolumes_)
    {
        if (subVolume->isEdge())
        {
            edgePoints.append(subVolume->midpoint());
        }
    }

    return edgePoints;
}
//---------------------------------------------------------------------------//
