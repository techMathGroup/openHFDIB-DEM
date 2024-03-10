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
    DEM parameters

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef parameters_H
#define parameters_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "materialProperties.H"
#include "interAdhesion.H"
#include "wallMatInfo.H"
#include "wallPlaneInfo.H"
#include "virtualMeshLevel.H"
#include "cyclicPlaneInfo.H"
#include "contactModelInfo.H"
#include "solverInfo.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Parameters Definition
\*---------------------------------------------------------------------------*/

HashTable<materialInfo,string,Hash<string>> materialProperties::materialProps_ =
    HashTable<materialInfo,string,Hash<string>>();

HashTable<scalar,string,Hash<string>> interAdhesion::matInterAdh_ =
    HashTable<scalar,string,Hash<string>>();

HashTable<materialInfo,string,Hash<string>> wallMatInfo::wallMatInfo_ =
    HashTable<materialInfo,string,Hash<string>>();

HashTable<List<vector>,string,Hash<string>> wallPlaneInfo::wallPlaneInfo_ =
    HashTable<List<vector>,string,Hash<string>>();

HashTable<List<vector>,string,Hash<string>> cyclicPlaneInfo::cyclicPlaneInfo_ =
    HashTable<List<vector>,string,Hash<string>>();

std::map<string, string> cyclicPlaneInfo::cyclicPatchPairs_ =
    std::map<string, string>();

scalar virtualMeshLevel::virtualMeshLevel_ = 1;
scalar virtualMeshLevel::levelOfDivision_ = 1;
scalar virtualMeshLevel::charCellSize_ = 1;
scalar contactModelInfo::LcCoeff_ = 4.0;
label  solverInfo::nSolidsTreshnold_ = 1000;
bool solverInfo::onlyDEM_ = false;
bool contactModelInfo::useChenRotationalModel_ = false;
bool contactModelInfo::useMindlinRotationalModel_ = false;
bool contactModelInfo::increasedDamping_ = false;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //