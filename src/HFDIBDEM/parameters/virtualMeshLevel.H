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
    contactModel

Description
    Static class for storing virtual mesh level

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef virtualMeshLevel_H
#define virtualMeshLevel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class virtualMeshLevel Declaration
\*---------------------------------------------------------------------------*/


class virtualMeshLevel
{
    public:

    static void setVirtualMeshLevel(scalar virtualMeshLevel,scalar charCellSize)
    {
        virtualMeshLevel_ = virtualMeshLevel;
        levelOfDivision_ = int(pow(2,virtualMeshLevel_));
        charCellSize_ = charCellSize;
    }

    static scalar getVirtualMeshLevel()
    {
        return virtualMeshLevel_;
    }

    static scalar getLevelOfDivision()
    {
        return levelOfDivision_;
    }

    static scalar getCharCellSize()
    {
        return charCellSize_;
    }

    static void setMaxSubVolumes(const scalar maxSubVolumes)
    {
        maxSubVolumes_ = maxSubVolumes;
    }

    static scalar getMaxSubVolumes()
    {
        return maxSubVolumes_;
    }

    static scalar virtualMeshLevel_;
    static scalar levelOfDivision_;
    static scalar charCellSize_;

    // Upper bound on the number of sub-volumes a single wall virtual mesh
    // may span. Exceeding it triggers a FatalError at virtual-mesh creation
    // (see wallSubContactInfo), and bounds the flood-fill iteration caps
    // (see virtualMeshWall). Configurable via the virtualMesh dict entry
    // maxSubVolumes.
    static scalar maxSubVolumes_;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
