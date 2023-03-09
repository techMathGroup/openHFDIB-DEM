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
    wallContactInfo.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef spectatorMesh_H
#define spectatorMesh_H

#include "dictionary.H"
#include "fvCFD.H"
#include "ibContactClass.H"
#include "sMProperties.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ibContactClass Declaration
\*---------------------------------------------------------------------------*/
class spectatorMesh
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

    public:
    
    const vector matrixSize_;
    const boundBox bBox_;
    const scalar charCellSize_;

    List<List<List<autoPtr<sMProperties>>>> centroidMatrix_;
    List<List<List<autoPtr<sMProperties>>>> vertexMatrix_;
    // Constructors ------------------------------------------------

    spectatorMesh
    (
    const vector matrixSize,
    const boundBox bBox,
    const scalar charCellSize
    );

    // Destructors -------------------------------------------------
    ~spectatorMesh();

//-----------------------------operator-------------------------------//
    sMProperties& operator[](const vector& sVI)
    {
        if (!centroidMatrix_[sVI[0]][sVI[1]][sVI[2]].valid())
        {
            centroidMatrix_[sVI[0]][sVI[1]][sVI[2]].set(
                new sMProperties(getCentroidPoint(sVI))
            );
        }
        return centroidMatrix_[sVI[0]][sVI[1]][sVI[2]]();
    }
//-----------------------------operator-------------------------------// 
   sMProperties& operator()(const vector& sVI)
    {
        if (!vertexMatrix_[sVI[0]][sVI[1]][sVI[2]].valid())
        {
            vertexMatrix_[sVI[0]][sVI[1]][sVI[2]].set(
                new sMProperties(getVertexPoint(sVI))
            );
        }
        return vertexMatrix_[sVI[0]][sVI[1]][sVI[2]]();
    }
//-----------------------------operator-------------------------------//
//-----------------------------Func-------------------------------//
    point getCentroidPoint
    (
        const vector& elementIndex
    );

    point getVertexPoint
    (
        const point& elementIndex
    );

    vector getSMCentroidIndex
    (
        const point& pointInDomain
    );

    vector getSMVertexIndex
    (
        const point& pointInDomain
    );

    List<vector> faceNeighbourElements
    (
        vector& elementIndex
    );    

    List<vector> elementVertexIndexies
    (
        vector& elementIndex
    );

    bool isPointInElementBB
    (
        point& pointToCheck,
        vector& elementIndex
    );

    boundBox getElementBB
    (
        vector& elementIndex
    );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
