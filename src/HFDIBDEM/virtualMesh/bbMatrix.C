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
#include "bbMatrix.H"

#include "virtualMeshLevel.H"

using namespace Foam;

//---------------------------------------------------------------------------//
bbMatrix::bbMatrix
(
    const vector matrixSize,
    const boundBox bBox,
    const scalar& charCellSize,
    const scalar& subVolumeV
)
:
matrixSize_(matrixSize),
bBox_(bBox),
charCellSize_(charCellSize),
subVolumeV_(subVolumeV)
{
    // Note (MI): here, implementation was changed by PK, probrably
    //            for clarity? -> ASK
    
    
    // bbMatrix_ = List<List<List<autoPtr<subVolumeProperties>>>>(matrixSize_[0],
        // List<List<autoPtr<subVolumeProperties>>>(matrixSize_[1],
        // List<autoPtr<subVolumeProperties>>(matrixSize_[2])));

    bbMatrix_.setSize(matrixSize_[0]);
    forAll(bbMatrix_, i)
    {
        bbMatrix_[i].setSize(matrixSize_[1]);
        forAll(bbMatrix_[i], j)
        {
            bbMatrix_[i][j].setSize(matrixSize_[2]);
        }
    }        
}
bbMatrix::~bbMatrix()
{
}
//---------------------------------------------------------------------------//
vector bbMatrix::getPointInMesh
(
	const vector& subVolumeIndex
)
{
    return vector(
        bBox_.min()[0] + (2*subVolumeIndex[0]+1)*(charCellSize_/virtualMeshLevel::getLevelOfDivision())*0.5,
        bBox_.min()[1] + (2*subVolumeIndex[1]+1)*(charCellSize_/virtualMeshLevel::getLevelOfDivision())*0.5,
        bBox_.min()[2] + (2*subVolumeIndex[2]+1)*(charCellSize_/virtualMeshLevel::getLevelOfDivision())*0.5
    );
}
//---------------------------------------------------------------------------//
vector bbMatrix::getSVIndexForPoint
(
    const point& pointInDomain
)
{
    vector subVolumeIndex(vector::zero);
    for(label i = 0;i<3;i++)
    {
        subVolumeIndex[i] = floor((pointInDomain[i]-bBox_.min()[i])/charCellSize_*virtualMeshLevel::getLevelOfDivision());
        if(subVolumeIndex[i] >= matrixSize_[i] || subVolumeIndex[i] < 0 )
        {
            subVolumeIndex = (vector::one)*(-1);
            break;
        }
    }
    return subVolumeIndex;
}
//---------------------------------------------------------------------------//
vector bbMatrix::getSVIndexForPoint_Wall
(
    point pointInDomain
)
{
    for(label i = 0;i<3;i++)
    {
        if(mag(pointInDomain[i])<SMALL)
        {
            pointInDomain[i] = 0.0;
        }
    } 
    vector subVolumeIndex(vector::zero);
    for(label i = 0;i<3;i++)
    {
        subVolumeIndex[i] = floor((pointInDomain[i]-bBox_.min()[i])/charCellSize_*virtualMeshLevel::getLevelOfDivision());
        if(subVolumeIndex[i] >= matrixSize_[i] || subVolumeIndex[i] < 0 )
        {
            return(getSVIndexForPoint(bBox_.midpoint()));
            break;
        }
    }
    return subVolumeIndex;
}
//---------------------------------------------------------------------------//
vector bbMatrix::getFirstSubVolumeIndex
(
    point& subVolumePoint,
    bool& isInMatrix
)
{
    vector subVolumeIndex(vector::zero);
    isInMatrix = true;
    for(label i = 0;i<3;i++)
    {
        subVolumeIndex[i] = floor((subVolumePoint[i]-bBox_.min()[i])
            /charCellSize_*virtualMeshLevel::getLevelOfDivision());

        if(subVolumeIndex[i] >= matrixSize_[i] || subVolumeIndex[i] < 0 )
        {
            subVolumeIndex = (vector::one)*(-1);
            isInMatrix = false;
            break;
        }
    }
    return subVolumeIndex;
}
//---------------------------------------------------------------------------//
List<vector> bbMatrix::faceNeighbourSubVolumes
(
    vector& subVolumeIndex
)
{
    const label i(subVolumeIndex[0]);
    const label j(subVolumeIndex[1]);
    const label k(subVolumeIndex[2]);

    vector neighbour1(i-1,j,k);
    vector neighbour2(i+1,j,k);
    vector neighbour3(i,j-1,k);
    vector neighbour4(i,j+1,k);
    vector neighbour5(i,j,k-1);
    vector neighbour6(i,j,k+1);

    const List<vector> NLL{{neighbour1,neighbour2,neighbour3,neighbour4,neighbour5,neighbour6}};
    List<vector> cellNeighbours;

    forAll(NLL, neig)
    {
        if((NLL[neig][0] >= 0 && NLL[neig][0] < matrixSize_[0])&&
            (NLL[neig][1] >= 0 && NLL[neig][1] < matrixSize_[1])&&
            (NLL[neig][2] >= 0 && NLL[neig][2] < matrixSize_[2]))
        {
            cellNeighbours.append(NLL[neig]);
        }
    }

	return cellNeighbours;
}
//---------------------------------------------------------------------------//
List<vector> bbMatrix::edgeNeighbourSubVolumes //returns only edgeNeighbours without faceNeighbours
(
    vector& subVolumeIndex
)
{
    const label i(subVolumeIndex[0]);
    const label j(subVolumeIndex[1]);
    const label k(subVolumeIndex[2]);

    vector edgeNeighbour1(i-1,j,k-1);
    vector edgeNeighbour2(i+1,j,k-1);
    vector edgeNeighbour3(i+1,j,k+1);
    vector edgeNeighbour4(i-1,j,k+1);

    vector edgeNeighbour5(i-1,j+1,k);
    vector edgeNeighbour6(i-1,j-1,k);
    vector edgeNeighbour7(i+1,j-1,k);
    vector edgeNeighbour8(i+1,j+1,k);

    vector edgeNeighbour9(i,j+1,k+1);
    vector edgeNeighbour10(i,j+1,k-1);
    vector edgeNeighbour11(i,j-1,k-1);
    vector edgeNeighbour12(i,j-1,k+1);

    const List<vector> NLL{{edgeNeighbour1,edgeNeighbour2,edgeNeighbour3,edgeNeighbour4,
                            edgeNeighbour5,edgeNeighbour6,edgeNeighbour7,edgeNeighbour8,
                            edgeNeighbour9,edgeNeighbour10,edgeNeighbour11,edgeNeighbour12}};
    List<vector> cellNeighbours;

    forAll(NLL, neig)
    {
        if((NLL[neig][0] >= 0 && NLL[neig][0] < matrixSize_[0])&&
            (NLL[neig][1] >= 0 && NLL[neig][1] < matrixSize_[1])&&
            (NLL[neig][2] >= 0 && NLL[neig][2] < matrixSize_[2]))
        {
            cellNeighbours.append(NLL[neig]);
        }
    }

	return cellNeighbours;
}
//---------------------------------------------------------------------------//
List<vector> bbMatrix::cornerNeighbourSubVolumes
(
    vector& subVolumeIndex
)
{
    label i(subVolumeIndex[0]);
    label j(subVolumeIndex[1]);
    label k(subVolumeIndex[2]);

    vector cornerNeighbour1(i+1,j+1,k+1);
    vector cornerNeighbour2(i+1,j-1,k+1);
    vector cornerNeighbour3(i-1,j-1,k+1);
    vector cornerNeighbour4(i-1,j+1,k+1);

    vector cornerNeighbour5(i+1,j+1,k-1);
    vector cornerNeighbour6(i+1,j-1,k-1);
    vector cornerNeighbour7(i-1,j-1,k-1);
    vector cornerNeighbour8(i-1,j+1,k-1);

    List<vector> NLL{{cornerNeighbour1,cornerNeighbour2,cornerNeighbour3,cornerNeighbour4,
                    cornerNeighbour5,cornerNeighbour6,cornerNeighbour7,cornerNeighbour8}};

    List<vector> cellNeighbours;
    forAll(NLL, neig)
    {
        if((NLL[neig][0] >= 0 && NLL[neig][0] < matrixSize_[0])&&
            (NLL[neig][1] >= 0 && NLL[neig][1] < matrixSize_[1])&&
            (NLL[neig][2] >= 0 && NLL[neig][2] < matrixSize_[2]))
        {
            cellNeighbours.append(NLL[neig]);
        }
    }

	return cellNeighbours;
}
//---------------------------------------------------------------------------//
