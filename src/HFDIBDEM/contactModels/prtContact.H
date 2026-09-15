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
    methods to solve particle-particle contact

SourceFiles
    prtContact.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef prtPrtContact_H
#define prtPrtContact_H

#include "dictionary.H"
#include "fvCFD.H"
// #include "fvMesh.H"

#include "ibContactClass.H"
#include "prtContactVars.H"
#include "prtContactInfo.H"
#include "periodicBody.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Namespace contactModel functions Declaration
\*---------------------------------------------------------------------------*/

namespace contactModel
{
    typedef DynamicList<label>  DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point>  DynamicPointList;

    void getContacts(
        const fvMesh&   mesh,
        prtContactInfo& prtcInfo
    );

    void getContacts_Cluster(
        const fvMesh&   mesh,
        prtContactInfo& prtcInfo
    );

    bool detectPrtPrtContact
    (
        const fvMesh& mesh,
        ibContactClass& cClass,
        ibContactClass& tClass,
        prtSubContactInfo& subCInfo
    );

    void getContacts_ArbShape(
        const fvMesh&   mesh,
        prtContactInfo& prtcInfo
    );

    bool detectPrtPrtContact_ArbShape
    (
        const fvMesh& mesh,
        ibContactClass& cClass,
        ibContactClass& tClass,
        prtSubContactInfo& subCInfo
    );

    bool detectPrtPrtContact_Sphere
    (
        const fvMesh& mesh,
        ibContactClass& cClass,
        ibContactClass& tClass
    );

    bool detectPrtPrtContact_Cluster
    (
        const fvMesh& mesh,
        ibContactClass& cClass,
        ibContactClass& tClass,
        prtSubContactInfo& subCInfo
    );

    void getPrtContactVars(
        const fvMesh&   mesh,
        ibContactClass& cClass,
        ibContactClass& tClass,
        prtSubContactInfo& subCInfo
    );

    void getPrtContactVars_ArbShape(
        const fvMesh&   mesh,
        ibContactClass& cClass,
        ibContactClass& tClass,
        prtSubContactInfo& subCInfo
    );

    void getPrtContactVars_Sphere
    (
        const fvMesh&   mesh,
        ibContactClass& cClass,
        ibContactClass& tClass,
        prtSubContactInfo& subCInfo
    );

    void getPrtContactVars_Cluster
    (
        const fvMesh&   mesh,
        ibContactClass& cClass,
        ibContactClass& tClass,
        prtSubContactInfo& subCInfo
    );

    bool solvePrtContact(
        const fvMesh&   mesh,
        prtContactInfo& cInfo,
        prtSubContactInfo& subCInfo,
        scalar deltaT
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
