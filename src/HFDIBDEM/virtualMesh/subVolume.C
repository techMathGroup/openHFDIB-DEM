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
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "subVolume.H"

using namespace Foam;

//---------------------------------------------------------------------------//
subVolume::subVolume()
:
treeBoundBox(boundBox()),
parentSV_(),
cVolumeInfo_(volumeType::UNKNOWN),
tVolumeInfo_(volumeType::UNKNOWN),
isEdge_(false)
{
}

subVolume::subVolume
(
    const boundBox bb
)
:
treeBoundBox(bb),
parentSV_(),
cVolumeInfo_(volumeType::UNKNOWN),
tVolumeInfo_(volumeType::UNKNOWN),
isEdge_(false)
{
}

subVolume::subVolume
(
    const boundBox bb,
    const std::shared_ptr<subVolume>& parentSV,
    const volumeType cVolumeType,
    const volumeType tVolumeType
)
:
treeBoundBox(bb),
parentSV_(parentSV),
cVolumeInfo_(cVolumeType),
tVolumeInfo_(tVolumeType),
isEdge_(false)
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
//~ List<subVolume>& subVolume::childSubVolumes()
List<autoPtr<subVolume>>& subVolume::childSubVolumes()
{
    if (childSubVolumes_.size() == 0)
    {
        std::shared_ptr<subVolume> parentSV = std::make_shared<subVolume>(*this);
        for (direction octant = 0; octant < 8; octant++)
        {
            childSubVolumes_.append
            (
                autoPtr<subVolume>
                (
                    new subVolume
                    (
                        subBbox(octant),
                        parentSV,
                        cVolumeInfo_.volumeType_,
                        tVolumeInfo_.volumeType_
                    )
                )
            );
        }
    }

    return childSubVolumes_;
}
//---------------------------------------------------------------------------//
bool subVolume::hasChildSubVolumes() const
{
    return childSubVolumes_.size() > 0;
}
//---------------------------------------------------------------------------//
std::shared_ptr<subVolume> subVolume::parentSV() const
{
    return parentSV_.lock();
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
void subVolume::setAsEdge()
{
    isEdge_ = true;
}
//---------------------------------------------------------------------------//
bool subVolume::isEdge() const
{
    return isEdge_;
}
//---------------------------------------------------------------------------//
// --- Copy constructor (deep-copies List<autoPtr<subVolume>>)
subVolume::subVolume(const subVolume& other)
:
    // copy the treeBoundBox part
    treeBoundBox(static_cast<const treeBoundBox&>(other)),
    // shared_ptr can be copied directly (shared ownership of parent)
    parentSV_(other.parentSV_),
    // create an empty list; we'll size & fill below
    childSubVolumes_(),
    // cheap copies (ibSubVolumeInfo is now copyable)
    cVolumeInfo_(other.cVolumeInfo_),
    tVolumeInfo_(other.tVolumeInfo_),
    isEdge_(other.isEdge_)
{
    // Deep copy children: clone each pointed subVolume if present
    const label n = other.childSubVolumes_.size();
    childSubVolumes_.setSize(n);

    for (label i = 0; i < n; ++i)
    {
        if (other.childSubVolumes_[i].valid())
        {
            childSubVolumes_[i].reset
            (
                new subVolume(*other.childSubVolumes_[i]) // recursive copy
            );
        }
        else
        {
            childSubVolumes_[i].clear(); // keep invalid (null) autoPtr
        }
    }
}

// --- Copy assignment (deep-copies List<autoPtr<subVolume>>)
subVolume& subVolume::operator=(const subVolume& other)
{
    if (this != &other)
    {
        // copy base
        static_cast<treeBoundBox&>(*this) = static_cast<const treeBoundBox&>(other);

        // copy trivials
        parentSV_   = other.parentSV_;
        cVolumeInfo_ = other.cVolumeInfo_;
        tVolumeInfo_ = other.tVolumeInfo_;
        isEdge_      = other.isEdge_;

        // rebuild children (strong exception-safety not strictly required here,
        // but this is reasonably safe with OpenFOAM containers)
        const label n = other.childSubVolumes_.size();
        childSubVolumes_.setSize(n);

        for (label i = 0; i < n; ++i)
        {
            if (other.childSubVolumes_[i].valid())
            {
                childSubVolumes_[i].reset
                (
                    new subVolume(*other.childSubVolumes_[i])
                );
            }
            else
            {
                childSubVolumes_[i].clear();
            }
        }
    }
    return *this;
}
