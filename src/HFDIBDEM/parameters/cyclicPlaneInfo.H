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
    Static class for storing wall material information

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef cyclicPlaneInfo_H
#define cyclicPlaneInfo_H

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class wallPlaneInfo Declaration
\*---------------------------------------------------------------------------*/


class cyclicPlaneInfo
{
    public:

    static const HashTable<List<vector>,string,Hash<string>>& getCyclicPlaneInfo()
    {
        return cyclicPlaneInfo_;
    }

    static const vector getCyclicTransVec(const string& cyclicPatch)
    {
        const string& cyclicPatchNeighbour = cyclicPatchPairs_[cyclicPatch];

        vector transVec = cyclicPlaneInfo_[cyclicPatchNeighbour][1] - cyclicPlaneInfo_[cyclicPatch][1];

        Info << "transVec1: " << transVec << endl;

        transVec = cyclicPlaneInfo_[cyclicPatch][0] * (transVec & cyclicPlaneInfo_[cyclicPatch][0]);

        Info << "transVec2: " << transVec << endl;

        return transVec;
    }

    static void insert
    (
        const string& cyclicPatch,
        const vector& patchNVec,
        const vector& planePoint,
        const string& cyclicPatchNeighbour
    )
    {
        List<vector> planeInfo{{patchNVec / mag(patchNVec),planePoint}};
        cyclicPlaneInfo_.insert(cyclicPatch,planeInfo);

        if (cyclicPatchPairs_.find(cyclicPatch) == cyclicPatchPairs_.end())
        {
            cyclicPatchPairs_.insert(std::pair<string,string>(cyclicPatch,cyclicPatchNeighbour));
        }
    };

    //- Hash table for storing wall material information
    static HashTable<List<vector>,string,Hash<string>> cyclicPlaneInfo_;
    static std::map<string, string> cyclicPatchPairs_;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
