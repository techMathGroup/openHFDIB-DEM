/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
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

InNamespace
    Foam

Description
    shared helpers for the virtual meshes (pseudo-2D clipping and
    overflow guards), used by both the wall and the particle-particle
    contact paths

    Note (MI): current idea is to potentially use this class for 
               future additional helpers

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef virtualMeshTools_H
#define virtualMeshTools_H

#include "fvCFD.H"
#include "boundBox.H"
#include "virtualMeshLevel.H"
#include "outputHFDIBDEM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Pseudo-2D: clip a virtual mesh bBox to a single sub-volume layer in the
// empty direction and return the volume/area rescaling factor. The contact
// patch of an extruded body is uniform along the empty direction, so
// flooding one layer and multiplying the result by
// emptyScale = originalSpan/subVolumeEdge reproduces the 3D volume (and,
// after division by subVolumeV, the 3D patch area) exactly, while removing
// the bodyThickness/subVolumeEdge factor from every traversal. Bodies
// already thinner than ~1.5 sub-volume edges (incl. the zero-span case)
// are left untouched with emptyScale = 1. The starting point is clamped
// into the clipped slab. No-op for 3D cases.
scalar clipEmptyDirection(boundBox& bB, point& startingPoint);

// Stop before a virtual mesh that could never be processed in a reasonable
// time (or memory) gets allocated. The bound is
// virtualMeshLevel::maxSubVolumes_ (virtualMesh dict entry maxSubVolumes).
// The sub-volume count is computed in double to avoid 32-bit label
// overflow.
void checkVMSize
(
    const vector& subVolumeNVector,
    const boundBox& bB,
    const word& what
);

// Same guard for mesh-less (octree) virtual meshes, where the worst-case
// leaf count replaces the matrix sub-volume count.
void checkVMLeafCount
(
    const scalar& nLeaves,
    const boundBox& bB,
    const word& what
);

// Upper bound on the number of sub-volumes a single traversal over a
// virtual mesh may visit. With the push-guard each sub-volume is visited
// at most once, so a full scan is a guaranteed terminating bound that never
// truncates a legitimate search. The bound is additionally limited by
// virtualMeshLevel::maxSubVolumes_ so that an accidentally huge virtual
// mesh cannot monopolise a time step. Computed in double to avoid 32-bit
// label overflow on large meshes.
label maxVSIter(const scalar& nSV);

// Matrix-flavoured overload (kept for the wall flood-fills).
label maxVSIter(const vector& matrixSize);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
