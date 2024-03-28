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
    contactModel

Description
    methods to solve particle-particle contact

SourceFiles
    virtualMesh.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef virtualMeshWall_H
#define virtualMeshWall_H

#include "dictionary.H"
#include "fvCFD.H"
#include "geomModel.H"
#include "bbMatrix.H"
#include "clockTime.H"
#include "subVolumeProperties.H"
// #include ""


// #include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct virtualMeshWallInfo
{
    virtualMeshWallInfo():
    bBox(boundBox()),
    startingPoint(vector::zero),
    subVolumeNVector(vector::zero),
    charCellSize(0),
    subVolumeV(0)

    {}
    virtualMeshWallInfo(
        boundBox bBox,
        point startingPoint,
        vector subVolumeNVector,
        scalar charCellSize,
        scalar subVolumeV):
    bBox(bBox),
    startingPoint(startingPoint),
    subVolumeNVector(subVolumeNVector),
    charCellSize(charCellSize),
    subVolumeV(subVolumeV)
    {}

    const point& getStartingPoint()
    {
        return startingPoint;
    }

    const scalar& getSVVolume()
    {
        return subVolumeV;
    }

    const scalar& getCharCellSize()
    {
        return charCellSize;
    }

    void shiftBBox
    (
        vector shiftVec
    )
    {
        bBox.min() += shiftVec;
        bBox.max() += shiftVec;
    }

    boundBox& getbBox()
    {
        return bBox;
    }

    void setStartingPoint
    (
        point sPoint
    )
    {
        startingPoint = sPoint;
    }

    boundBox bBox;
    point startingPoint;
    vector subVolumeNVector;
    scalar charCellSize;
    scalar subVolumeV;    
};

/*---------------------------------------------------------------------------*\
                        Class virtualMeshWall Declaration
\*---------------------------------------------------------------------------*/
class virtualMeshWall
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;
    typedef HashSet<vector,Hash<vector>> vectorHashSet;

    geomModel& cGeomModel_;

    virtualMeshWallInfo& vMeshWallInfo_;

    bbMatrix bbMatrix_;

    DynamicPointList edgeSubVolumesPoints_;

    public:
    
    vector contactCenter_ = vector::zero;

    // Constructors ------------------------------------------------

    virtualMeshWall
    (
        virtualMeshWallInfo& vMeshWallInfo,
        geomModel& cGeomModel
    );


    // Destructors -------------------------------------------------
    ~virtualMeshWall();

    bool detectFirstContactPoint();

    bool detectFirstFaceContactPoint();

    scalar evaluateContact();

    void checkSubVolume(subVolumeProperties& subVolume);

    void resetSubVolume(subVolumeProperties& subVolume);

    vector shiftEdgeSubVolume
    (
        vector& subVolumeIndex
    );

    label getInternalSV();

    DynamicPointList& getEdgeSVPoints()
    {
        return edgeSubVolumesPoints_;
    }

    vector& getContactCenter()
    {
        return contactCenter_;
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
