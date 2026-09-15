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
    wallContact.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef wallContact_H
#define wallContact_H

#include "dictionary.H"
#include "fvCFD.H"
// #include "fvMesh.H"

#include "ibContactClass.H"
#include "materialInfo.H"
#include "wallContactInfo.H"
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

    bool detectWallContact
    (
        const fvMesh&   mesh,
        ibContactClass& ibClass,
        wallContactInfo& wallCntInfo
    );

    bool detectWallContact_ArbShape
    (
        const fvMesh&   mesh,
        ibContactClass& ibClass,
        wallContactInfo& wallCntInfo
    );

    bool detectWallContact_Sphere
    (
        const fvMesh&   mesh,
        ibContactClass& ibClass,
        wallContactInfo& wallCntInfo
    );

    bool detectWallContact_Cluster
    (
        const fvMesh&   mesh,
        ibContactClass& ibClass,
        wallContactInfo& wallCntInfo
    );

    void getWallContactVars
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo,
        const scalar deltaT,
        wallSubContactInfo& sWC
    );

    void getWallContactVars_ArbShape
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo,
        const scalar deltaT,
        wallSubContactInfo& sWC
    );

    void getWallContactVars_Sphere
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo,
        const scalar deltaT,
        wallSubContactInfo& sCW
    );

    void getWallContactVars_Cluster
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo,
        const scalar deltaT,
        wallSubContactInfo& sWC
    );

    DynamicList<Tuple2<label,string>> getContactFacesSphere
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo,
        label faceLabel,
        DynamicList<Tuple2<label,string>>& checkedFaces
    );

    // DynamicList<Tuple2<label,string>> getContactFacesArbShape
    // (
    //     const fvMesh&   mesh,
    //     wallContactInfo& wallCntInfo,
    //     label faceLabel,
    //     DynamicList<Tuple2<label,string>>& checkedFaces,
    //     HashTable<bool,label,Hash<label>>& pointTable,
    //     DynamicVectorList& center,
    //     DynamicVectorList& normal,
    //     DynamicLabelList& centerPoints,
    //     DynamicScalarList& area
    // );

    bool solveWallContact
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo,
        scalar deltaT,
        wallSubContactInfo& sWC
    );

    // scalar getInterVolume_ArbShape
    // (
    //     const fvMesh&   mesh,
    //     ibContactClass& cClass,
    //     ibContactVars& cVars,
    //     vector nVec,
    //     vector center
    // );

    scalar getInterVolume_Sphere
    (
        const fvMesh&   mesh,
        ibContactClass& cClass,
        ibContactVars& cVars,
        plane& cPlane
    );

    scalar sphereContactArea
    (
        const fvMesh&   mesh,
        ibContactClass& cClass,
        ibContactVars& cVars,
        plane& cPlane
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
