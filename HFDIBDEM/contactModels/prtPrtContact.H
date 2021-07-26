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
    prtPrtContact.C

Contributors
    Martin Isoz (2019-*), Martin Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef prtPrtContact_H
#define prtPrtContact_H

#include "dictionary.H"
#include "fvCFD.H"
#include "dynamicFvMesh.H"

#include "contactInfo.H"
//#include "wallContact.H"
#include "prtPrtContactVars.H"

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

    bool detectPrtPrtContact
    (
        const dynamicFvMesh& mesh,
        contactInfo& cInfo,
        contactInfo& tInfo
    );

    template <contactType cT, contactType tT>
    bool detectPrtPrtContact
    (
        const dynamicFvMesh& mesh,
        contactInfo& cInfo,
        contactInfo& tInfo
    );

    template <>
    bool detectPrtPrtContact <sphere,sphere>
    (
        const dynamicFvMesh& mesh,
        contactInfo& cInfo,
        contactInfo& tInfo
    );

    template <contactType cT, contactType tT>
    void getPrtPrtContactVars(
        const dynamicFvMesh&   mesh,
        contactInfo& cInfo,
        contactInfo& tInfo,
        vector geometricD,
        prtPrtContactVars& prtPrtCntVars
    );

    Tuple2<scalar,vector> get3DcontactInfo(
        const dynamicFvMesh&   mesh,
        DynamicLabelList commonCells,
        DynamicVectorList edgePoints,
        vector normalVector,
        vector contactCenter,
        label owner
    );

    Tuple2<scalar,vector>  get2DcontactInfo(
        const dynamicFvMesh&   mesh,
        DynamicLabelList commonCells,
        vector normalVector,
        vector contactCenter,
        vector geometricD
    );

    template <>
    void getPrtPrtContactVars  <sphere,sphere>
    (
        const dynamicFvMesh&   mesh,
        contactInfo& cInfo,
        contactInfo& tInfo,
        vector geometricD,
        prtPrtContactVars& prtPrtCntVars
    );

    void solvePrtContact(
        const dynamicFvMesh&   mesh,
        contactInfo& cInfo,
        contactVars* cVars,
        contactInfo& tInfo,
        contactVars* tVars,
        vector geometricD,
        scalar deltaT,
        Tuple2<Tuple2<vector,vector>,Tuple2<vector,vector>>& outVars
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "prtPrtContact.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
