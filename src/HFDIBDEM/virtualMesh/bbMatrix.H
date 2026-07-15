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

#ifndef bbMatrix_H
#define bbMatrix_H

#include "dictionary.H"
#include "fvCFD.H"
#include "ibContactClass.H"
#include "subVolumeProperties.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ibContactClass Declaration
\*---------------------------------------------------------------------------*/
class bbMatrix
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

    const vector matrixSize_;
    const boundBox bBox_;
    const scalar charCellSize_;
    const scalar subVolumeV_;
    List<List<List<autoPtr<subVolumeProperties>>>> bbMatrix_;

    public:

    // Constructors ------------------------------------------------

    bbMatrix
    (
    const vector matrixSize,
    const boundBox bBox,
    const scalar& charCellSize,
    const scalar& subVolumeV
    );

    // Destructors -------------------------------------------------
    ~bbMatrix();

//-----------------------------operator-------------------------------//
    subVolumeProperties& operator[](const vector& sVI)
    {
        if (!bbMatrix_[sVI[0]][sVI[1]][sVI[2]].valid())
        {
            bbMatrix_[sVI[0]][sVI[1]][sVI[2]].reset(                    //OF.com: set -> reset
                new subVolumeProperties(getPointInMesh(sVI))
            );
        }

        return bbMatrix_[sVI[0]][sVI[1]][sVI[2]]();
    }
//-----------------------------operator-------------------------------//
//-----------------------------Func-------------------------------//
    vector getPointInMesh
    (
        const vector& subVolumeIndex
    );

    vector getSVIndexForPoint
    (
        const point& pointInDomain
    );

    vector getSVIndexForPoint_Wall
    (
        point pointInDomain
    );

    vector getFirstSubVolumeIndex
    (
        point& subVolumePoint,
        bool& isInMatrix
    );

    List<vector> faceNeighbourSubVolumes
    (
        vector& subVolumeIndex
    );

    List<vector> edgeNeighbourSubVolumes
    (
        vector& subVolumeIndex
    );

    List<vector> cornerNeighbourSubVolumes
    (
        vector& subVolumeIndex
    );

    scalar getSubVolumeV()
    {
        return subVolumeV_;
    }

    scalar getCharCellSize()
    {
        return charCellSize_;
    }

    boundBox getBBox()
    {
        return bBox_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
