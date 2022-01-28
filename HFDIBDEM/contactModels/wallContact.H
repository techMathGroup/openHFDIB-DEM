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
    methods to solve particle-particle contact

SourceFiles
    wallContact.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef wallContact_H
#define wallContact_H

#include "dictionary.H"
#include "fvCFD.H"
// #include "fvMesh.H"

#include "ibContactClass.H"
#include "materialInfo.H"
#include "wallContactInfo.H"

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

    void detectWallContact
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo
    );

    template <contactType cT>
    void detectWallContact
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo
    );

    template <>
    void detectWallContact <sphere>
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo
    );

    void getWallContactVars
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo,
        const scalar deltaT
    );

    template <contactType cT>
    void getWallContactVars
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo,
        const scalar deltaT
    );

    template <>
    void getWallContactVars <sphere>
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo,
        const scalar deltaT
    );

    DynamicList<Tuple2<label,string>> getContactFacesSphere
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo,
        label faceLabel,
        DynamicList<Tuple2<label,string>>& checkedFaces
    );

    DynamicList<Tuple2<label,string>> getContactFacesArbShape
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo,
        label faceLabel,
        DynamicList<Tuple2<label,string>>& checkedFaces,
        HashTable<bool,label,Hash<label>>& pointTable,
        DynamicVectorList& center,
        DynamicVectorList& normal,
        DynamicLabelList& centerPoints,
        DynamicScalarList& area
    );

    void solveWallContact
    (
        const fvMesh&   mesh,
        wallContactInfo& wallCntInfo,
        scalar deltaT
    );

    template <contactType cT>
    scalar getInterVolume
    (
        const fvMesh&   mesh,
        ibContactClass& cClass,
        ibContactVars& cVars,
        vector nVec,
        vector center
    );

    template <>
    scalar getInterVolume <sphere>
    (
        const fvMesh&   mesh,
        ibContactClass& cClass,
        ibContactVars& cVars,
        vector nVec,
        vector center
    );

    scalar sphereContactArea
    (
        const fvMesh&   mesh,
        ibContactClass& cClass,
        ibContactVars& cVars,
        vector nVec,
        vector center
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "wallContact.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
