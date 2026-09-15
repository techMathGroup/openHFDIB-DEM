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
#include "subContact.H"

using namespace Foam;

//---------------------------------------------------------------------------//
subContact::subContact()
:
volume_(0)
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
