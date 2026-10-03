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
    subVolumeProperties

Description
    Classification flags of a single virtual-mesh sub-volume
    (contacting/target body containment, edge status), plus the
    exact-clip cache of the wall virtual-mesh path.

SourceFiles
    subVolume.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef subVolumeProperties_H
#define subVolumeProperties_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class ibContactClass Declaration
\*---------------------------------------------------------------------------*/
class subVolumeProperties
{
    public:

    subVolumeProperties(point center):
    center(center)
    {
        resetClip();
    }

    // const point center;
    const point center;
    bool isCBody=false;
    bool isTBody=false;
    bool toCheck=true;
    bool isOnEdge=false;

    // exact-clip cache (wall path): -1 = unevaluated, 0 = evaluated
    // empty; guard with > -0.5 as on the prt path (subVolume.H);
    // reset by resetSubVolume so a re-checked leaf re-evaluates
    scalar clippedVolume_ = -1;
    vector clippedCentroid_ = vector::zero;

    void resetClip()
    {
        clippedVolume_ = -1;
        clippedCentroid_ = vector::zero;
    }

    // point pointInMesh;
    // autoPtr<point> pointInMesh;

    // List<vector> faceNeighbours;
    // autoPtr<List<vector>> faceNeighbours;

    bool isContact()
    {
        return isCBody&&isTBody;
    };

    bool isInBody()
    {
        return  (isCBody||isTBody);
    }
    // Destructors -------------------------------------------------
    ~subVolumeProperties(){}
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
// ************************************************************************* //
