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
#include "subVolume.H"

using namespace Foam;

//---------------------------------------------------------------------------//
subVolume::subVolume()
{
}

subVolume::subVolume
(
    const boundBox bb
)
:
treeBoundBox(bb),
parentSV_(nullptr),
cVolumeInfo_(volumeType::unknown),
tVolumeInfo_(volumeType::unknown)
{
}

subVolume::subVolume
(
    const boundBox bb,
    const std::shared_ptr<subVolume> parentSV,
    const volumeType cVolumeType,
    const volumeType tVolumeType
)
:
treeBoundBox(bb),
parentSV_(parentSV),
cVolumeInfo_(cVolumeType),
tVolumeInfo_(tVolumeType)
{
}

subVolume::~subVolume()
{
}
//---------------------------------------------------------------------------//
ibSubVolumeInfo& subVolume::cVolumeInfo()
{
    return cVolumeInfo_;
}
//---------------------------------------------------------------------------//
ibSubVolumeInfo& subVolume::tVolumeInfo()
{
    return tVolumeInfo_;
}
//---------------------------------------------------------------------------//
List<subVolume>& subVolume::childSubVolumes()
{
    std::shared_ptr<subVolume> parentSV = std::make_shared<subVolume>(*this);
    if (childSubVolumes_.size() == 0)
    {
        for (direction octant = 0; octant < 8; octant++)
        {
            childSubVolumes_.append
            (
                subVolume
                (
                    subBbox(octant),
                    parentSV,
                    cVolumeInfo_.volumeType_,
                    tVolumeInfo_.volumeType_
                )
            );
        }
    }

    return childSubVolumes_;
}
//---------------------------------------------------------------------------//
std::shared_ptr<subVolume>& subVolume::parentSV()
{
    return parentSV_;
}
//---------------------------------------------------------------------------//
ibSubVolumeInfo& subVolume::getVolumeInfo(bool cIb)
{
    if (cIb)
    {
        return cVolumeInfo_;
    }
    else
    {
        return tVolumeInfo_;
    }
}
//---------------------------------------------------------------------------//
const ibSubVolumeInfo& subVolume::getVolumeInfo(bool cIb) const
{
    if (cIb)
    {
        return cVolumeInfo_;
    }
    else
    {
        return tVolumeInfo_;
    }
}
//---------------------------------------------------------------------------//
