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
scalar virtualMeshLevel::maxSubVolumes_ = 5e6;
scalar contactModelInfo::LcCoeff_ = 4.0;
label  solverInfo::nSolidsThreshold_ = 1000;
bool solverInfo::onlyDEM_ = false;
bool contactModelInfo::useChenRotationalModel_ = false;
bool contactModelInfo::useMindlinRotationalModel_ = false;
bool contactModelInfo::increasedDamping_ = false;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //