/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | |___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  | |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___/ |_____|_|  |_|
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

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#include "virtualMeshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar clipEmptyDirection
(
    boundBox& bB,
    point& startingPoint
)
{
    if (case3D || emptyDim < 0 || emptyDim > 2)
    {
        return 1.0;
    }

    scalar sv =
        virtualMeshLevel::getCharCellSize()
       /virtualMeshLevel::getLevelOfDivision();

    scalar origSpan = bB.span()[emptyDim];

    if (origSpan <= 1.5*sv)
    {
        return 1.0;
    }

    scalar mid = 0.5*(bB.min()[emptyDim] + bB.max()[emptyDim]);
    bB.min()[emptyDim] = mid - 0.5*sv;
    bB.max()[emptyDim] = mid + 0.5*sv;

    if
    (
        startingPoint[emptyDim] < bB.min()[emptyDim]
     || startingPoint[emptyDim] > bB.max()[emptyDim]
    )
    {
        startingPoint[emptyDim] = mid;
    }

    return origSpan/sv;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void checkVMLeafCount
(
    const scalar& nLeaves,
    const boundBox& bB,
    const word& what
)
{
    if (nLeaves <= max(virtualMeshLevel::getMaxSubVolumes(), scalar(1)))
    {
        return;
    }

    FatalErrorIn("void Foam::checkVMLeafCount()")
        << "The " << what << " virtual mesh could span " << nLeaves
        << " sub-volumes, exceeding the limit of "
        << virtualMeshLevel::getMaxSubVolumes() << endl
        << "    virtual mesh bBox: " << bB << endl
        << "    virtualMesh level: " << virtualMeshLevel::getVirtualMeshLevel()
        << " (levelOfDivision: "
        << virtualMeshLevel::getLevelOfDivision() << ")" << endl
        << "    virtualMesh charCellSize: "
        << virtualMeshLevel::getCharCellSize() << endl
        << "    sub-volume volume: "
        << pow(virtualMeshLevel::getCharCellSize()
        /virtualMeshLevel::getLevelOfDivision(),3) << endl
        << "Reduce the extent of the overlapping bounding box (e.g. body "
        << "thickness in the empty direction for pseudo-2D cases), lower "
        << "virtualMesh level, increase virtualMesh charCellSize, or raise "
        << "maxSubVolumes if this is intentional."
        << exit(FatalError);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void checkVMSize
(
    const vector& subVolumeNVector,
    const boundBox& bB,
    const word& what
)
{
    scalar nSV =
        max(subVolumeNVector.x(), scalar(1))
       *max(subVolumeNVector.y(), scalar(1))
       *max(subVolumeNVector.z(), scalar(1));

    checkVMLeafCount(nSV, bB, what);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label maxVSIter(const scalar& nSV)
{
    scalar nSVc = min(nSV, max(virtualMeshLevel::getMaxSubVolumes(), scalar(1)));

    return label(min(nSVc, scalar(labelMax)));
}

label maxVSIter(const vector& matrixSize)
{
    return maxVSIter
    (
        matrixSize.x()*matrixSize.y()*matrixSize.z()
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
