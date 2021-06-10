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
#include "geomModel.H"

using namespace Foam;

//---------------------------------------------------------------------------//
geomModel::geomModel
(
    const  dynamicFvMesh&   mesh,
    scalar  thrSurf,
    Vector<label> geometricD
)
:
mesh_(mesh),
owner_(0),
cellToStartInCreateIB_(0),
thrSurf_(thrSurf),
geometricD_(geometricD),
intSpan_(2.0),
sdBasedLambda_(false),
curMeshBounds_(mesh_.points(),false),
lastIbPoints_(0)
{    
    ibPartialVolume_.setSize(Pstream::nProcs()); 
}
geomModel::~geomModel()
{
}
