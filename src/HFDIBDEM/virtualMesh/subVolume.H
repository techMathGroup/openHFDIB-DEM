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

Class
    subVolume

Description
    Octree node of the virtual mesh: one cubic sub-volume
    with lazy-split children and per-body volume classification.

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef subVolume_H
#define subVolume_H

#include "treeBoundBox.H"
#include "volumeType.H"

#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct ibSubVolumeInfo
{
    ibSubVolumeInfo():
        volumeType_(volumeType::UNKNOWN)
    {}

    explicit ibSubVolumeInfo(volumeType vType)
    :
        volumeType_(volumeType::UNKNOWN)
    {
        switch(vType)
        {
            case volumeType::UNKNOWN:
            case volumeType::MIXED:
                volumeType_ = volumeType::UNKNOWN;
                break;
            case volumeType::INSIDE:
                volumeType_ = volumeType::INSIDE;
                break;
            case volumeType::OUTSIDE:
                volumeType_ = volumeType::OUTSIDE;
                break;
            default:
                volumeType_ = volumeType::UNKNOWN;
                break;
        }
    }

    volumeType volumeType_;
    autoPtr<labelList> shapesIn_;
    
    // --- copy ctor (deep copy of shapesIn_) ---
    ibSubVolumeInfo(const ibSubVolumeInfo& other)
    :
        volumeType_(other.volumeType_),
        shapesIn_
        (
            other.shapesIn_.valid()
            ? autoPtr<labelList>(new labelList(*other.shapesIn_))
            : autoPtr<labelList>()
        )
    {}

    // --- copy assignment (deep copy of shapesIn_) ---
    ibSubVolumeInfo& operator=(const ibSubVolumeInfo& other)
    {
        if (this != &other)
        {
            volumeType_ = other.volumeType_;

            // reset current storage
            shapesIn_.clear();

            if (other.shapesIn_.valid())
            {
                shapesIn_.reset(new labelList(*other.shapesIn_));
            }
            // else shapesIn_ remains empty (invalid)
        }
        return *this;
    }

    // --- move operations (transfer ownership of autoPtr) ---
    ibSubVolumeInfo(ibSubVolumeInfo&&) noexcept = default;
    ibSubVolumeInfo& operator=(ibSubVolumeInfo&&) noexcept = default;
};

/*---------------------------------------------------------------------------*\
                        Class subVolume Declaration
\*---------------------------------------------------------------------------*/
class subVolume : public treeBoundBox
{
    private:

    std::weak_ptr<subVolume> parentSV_;
    List<autoPtr<subVolume>> childSubVolumes_;                          //This is a List of autoPtr, which requires manual deep-copy

    ibSubVolumeInfo cVolumeInfo_;
    ibSubVolumeInfo tVolumeInfo_;

    bool isEdge_;

    public:

    // Constructors ------------------------------------------------
    subVolume();

    subVolume
    (
        const boundBox bb
    );

    subVolume
    (
        const boundBox bb,
        const std::shared_ptr<subVolume>& parentSV,
        const volumeType cVolumeType,
        const volumeType tVolumeType
    );

    // Destructors -------------------------------------------------
    ~subVolume();

    // Member Functions --------------------------------------------

    ibSubVolumeInfo& cVolumeInfo();
    ibSubVolumeInfo& tVolumeInfo();

    //~ List<subVolume>& childSubVolumes();
    List<autoPtr<subVolume>>& childSubVolumes();
    
    bool hasChildSubVolumes() const;

    std::shared_ptr<subVolume> parentSV() const;

    ibSubVolumeInfo& getVolumeInfo(bool cIb);
    const ibSubVolumeInfo& getVolumeInfo(bool cIb) const;

    void setAsEdge();
    bool isEdge() const;
        
    // Make copy/move operations available -------------------------
    subVolume(const subVolume& other);                                  //defined in .C
    subVolume& operator=(const subVolume& other);                       //defined in .C
    subVolume(subVolume&&) noexcept = default;
    subVolume& operator=(subVolume&&) noexcept = default;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
