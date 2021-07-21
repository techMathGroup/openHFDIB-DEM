/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________
                       | | | ||  ___|  _  \_   _| ___ \     H ybrid
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /     F ictitious
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \     D omain
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /     I mmersed
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/      B oundary
      | |
      |_|
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
    Martin Isoz (2019-*), Martin Šourek (2019-*), 
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "shapeBased.H"

using namespace Foam;

//---------------------------------------------------------------------------//
shapeBased::shapeBased
(
    const  dynamicFvMesh&   mesh,
    contactType cType,
    scalar  thrSurf,
    Vector<label> geometricD
)
:
geomModel(mesh,cType,thrSurf,geometricD)
{}
//---------------------------------------------------------------------------//
vector shapeBased::addModelReturnRandomPosition
(
    const bool allActiveCellsInMesh,
    const boundBox  cellZoneBounds,       
    Random&          randGen
)
{
    // Note (MI): this function will always return acceptable random
    //            position IF
    //            -> body boundBox is completely inside active boundBox
    //            AND
    //            -> active boundBox is completely contained in the mesh
    //
    // Note (MI): the check if body boundBox is inside active boundBox
    //            is simple and probably unecessary
    // Note (MI): an efficient check if all the active boundBox is inside
    //            mesh is an open issue at the moment
    vector ranVec(vector::zero);
    
    label nGeometricD(0);
    forAll (geometricD_, direction)
    {
        if (geometricD_[direction] == 1)
        {
            nGeometricD++;
        }
    }
    
    meshSearch searchEng(mesh_);
    
    // get its center of mass
    vector CoM(getCoM());
    
    const vector validDirs = (geometricD_ + Vector<label>::one)/2;
    vector dirCorr(cmptMultiply((vector::one - validDirs),CoM));
    dirCorr += cmptMultiply((vector::one - validDirs),0.5*(mesh_.bounds().max() + mesh_.bounds().min()));
    
    // get the body boundBox
    boundBox bodyBounds(getBounds());
    // compute the max scales to stay in active bounding box
    vector maxScales(cellZoneBounds.max() - bodyBounds.max());
    maxScales -= cellZoneBounds.min() - bodyBounds.min();
    maxScales *= 0.5*0.9;//0.Y is there just to be sure 
    
    Info << "-- addModelMessage-- " << "acceptable movements: " << maxScales << endl;
    
    scalar ranNum = 0;
    for (int i=0;i<3;i++)
    {
        ranNum = 2.0*maxScales[i]*randGen.scalar01() - 1.0*maxScales[i];
        ranVec[i] = ranNum;
    }
    
    ranVec = cmptMultiply(validDirs,ranVec);//translate only with respect to valid directions
    ranVec += dirCorr;
    
    return ranVec;
}
//---------------------------------------------------------------------------//
