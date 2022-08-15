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

    bool detectPrtPrtContact
    (
        const fvMesh& mesh,
        ibContactClass& cClass,
        ibContactClass& tClass
    );

    bool detectPrtPrtContact_ArbShape
    (
        const fvMesh& mesh,
        ibContactClass& cClass,
        ibContactClass& tClass
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
        ibContactClass& tClass
    );

    void getPrtContactVars(
        const fvMesh&   mesh,
        ibContactClass& cClass,
        ibContactClass& tClass,
        prtContactVars& prtCntVars
    );

    void getPrtContactVars_ArbShape(
        const fvMesh&   mesh,
        ibContactClass& cClass,
        ibContactClass& tClass,
        prtContactVars& prtCntVars
    );

    Tuple2<scalar,vector> get3DcontactVars(
        const fvMesh&   mesh,
        DynamicLabelList& commonCells,
        DynamicVectorList& edgePoints,
        vector normalVector,
        vector contactCenter,
        label owner
    );

    Tuple2<scalar,vector>  get2DcontactVars(
        const fvMesh&   mesh,
        DynamicLabelList& commonCells,
        vector normalVector,
        vector contactCenter
    );

    void getPrtContactVars_Sphere
    (
        const fvMesh&   mesh,
        ibContactClass& cClass,
        ibContactClass& tClass,
        prtContactVars& prtCntVars
    );

    void getPrtContactVars_Cluster
    (
        const fvMesh&   mesh,
        ibContactClass& cClass,
        ibContactClass& tClass,
        prtContactVars& prtCntVars
    );

    bool solvePrtContact(
        const fvMesh&   mesh,
        prtContactInfo& cInfo,
        scalar deltaT
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
