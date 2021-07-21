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
#include "contactInfo.H"

using namespace Foam;

//---------------------------------------------------------------------------//
contactInfo::contactInfo
(
    geomModel& geomModel,
    List<DynamicLabelList>& surfCells,
    List<DynamicLabelList>& intCells,
    scalar kN,
    scalar kt,
    scalar gammaN,
    scalar gammat,
    scalar mu,
    scalar adhN,
    scalar adhEqui
)
:
geomModel_(geomModel),
isInWallContact_(false),
surfCells_(surfCells),
intCells_(intCells),
inContactWithStatic_(false),
timeStepsInContWStatic_(0),
kN_(kN),
kt_(kt),
gammaN_(gammaN),
gammat_(gammat),
mu_(mu),
adhN_(adhN),
adhEqui_(adhEqui)
{    
    wallContactFaces_.setSize(Pstream::nProcs());
}
contactInfo::~contactInfo()
{
}
//---------------------------------------------------------------------------//
