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
#include "shapeBased.H"

using namespace Foam;

//---------------------------------------------------------------------------//
shapeBased::shapeBased
(
    const  fvMesh&   mesh,
    const contactType cType,
    scalar  thrSurf
)
:
geomModel(mesh,cType,thrSurf)
{}
//---------------------------------------------------------------------------//
vector shapeBased::addModelReturnRandomPosition
(
    const bool allActiveCellsInMesh,
    const boundBox  cellZoneBounds,
    Random&          randGen
)
{
    vector ranVec(vector::zero);

    meshSearch searchEng(mesh_);

    // get its center of mass
    vector CoM(getCoM());

    const vector validDirs = (geometricD + vector::one)/2;
    vector dirCorr(cmptMultiply((vector::one - validDirs),CoM));
    dirCorr += cmptMultiply((vector::one - validDirs),0.5*(mesh_.bounds().max() + mesh_.bounds().min()));

    // get the body boundBox
    boundBox bodyBounds(getBounds());
    // compute the max scales to stay in active bounding box
    vector maxScales(cellZoneBounds.max() - bodyBounds.max());
    maxScales -= cellZoneBounds.min() - bodyBounds.min();
    maxScales *= 0.5*0.9;//0.Y is there just to be sure

    InfoH << addModel_Info << "-- addModelMessage-- "
        << "acceptable movements: " << maxScales << endl;

    scalar ranNum = 0;
    for (int i=0;i<3;i++)
    {
        ranNum = 2.0*maxScales[i]*randGen.sample01<scalar>() - 1.0*maxScales[i];
        ranVec[i] = ranNum;
    }

    ranVec = cmptMultiply(validDirs,ranVec);                            //translate only with respect to valid directions
    ranVec += dirCorr;

    return ranVec;
}
//---------------------------------------------------------------------------//
