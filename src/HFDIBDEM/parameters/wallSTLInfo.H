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
    Static class for storing stlBased collision patch information

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef wallSTLInfo_H
#define wallSTLInfo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "word.H"
#include "HashTable.H"
#include "string.H"
#include "dictionary.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class wallSTLInfo Declaration
\*---------------------------------------------------------------------------*/


class wallSTLInfo
{
    public:

        static const HashTable<List<word>,string,Hash<string>>& getWallSTLInfo()
        {
            return wallSTLInfo_;
        }

        static void wallSTLInfo_insert
        (
            const string& wallPatch,
            const word& stlName,
            const word& material,
            const word& bodyGeom
        )
        {
            List<word> stlInfo{{stlName,material,bodyGeom}};
            wallSTLInfo_.insert(wallPatch,stlInfo);
        };

        //- cover request per stlBased patch; the dict is stored
        //- verbatim and parsed later by the immersedBody machinery
        static const HashTable<dictionary,string,Hash<string>>&
        getWallSTLCoverInfo()
        {
            return wallSTLInfoCover_;
        }

        static void wallSTLInfo_insertCover
        (
            const string& wallPatch,
            const dictionary& coverDict
        )
        {
            wallSTLInfoCover_.insert(wallPatch,coverDict);
        };

        //- Hash table for storing stlBased collision patch info
        static HashTable<List<word>,string,Hash<string>> wallSTLInfo_;

        //- Hash table for storing stlBased cover requests
        static HashTable<dictionary,string,Hash<string>> wallSTLInfoCover_;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
