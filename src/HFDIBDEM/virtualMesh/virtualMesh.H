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
    virtualMesh.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef virtualMesh_H
#define virtualMesh_H

#include "dictionary.H"
#include "fvCFD.H"
#include "geomModel.H"
#include "clockTime.H"
#include "subVolume.H"
#include "subContact.H"


// #include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct virtualMeshInfo
{
    virtualMeshInfo():
    sV(boundBox()),
    subVolumeV(0)
    {}

    virtualMeshInfo(const virtualMeshInfo& vMeshInfo):
    sV(vMeshInfo.sV),
    subVolumeV(vMeshInfo.subVolumeV)
    {
        if (vMeshInfo.startingPoint.valid())
        {
            startingPoint.reset(new point(vMeshInfo.startingPoint()));  //OF.com: set -> reset
        }
    }

    virtualMeshInfo(boundBox& bBox, scalar subVolumeV):
    sV(bBox),
    subVolumeV(subVolumeV)
    {}

    void operator= (const virtualMeshInfo& vMeshInfo)
    {
        sV = vMeshInfo.sV;
        subVolumeV = vMeshInfo.subVolumeV;

        if (vMeshInfo.startingPoint.valid())
        {
            startingPoint.reset(new point(vMeshInfo.startingPoint()));
        }
    }

    const point& getStartingPoint()
    {
        if (!startingPoint.valid())
        {
            startingPoint.reset(new point(sV.midpoint()));
        }

        return startingPoint();
    }

    subVolume& getSubVolume()
    {
        return sV;
    }

    subVolume sV;
    scalar subVolumeV;
    autoPtr<point> startingPoint;
};

/*---------------------------------------------------------------------------*\
                        Class virtualMesh Declaration
\*---------------------------------------------------------------------------*/
class virtualMesh
{
    private:

    typedef DynamicList<label>  DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point>  DynamicPointList;
    typedef HashSet<vector,Hash<vector>> vectorHashSet;

    geomModel& cGeomModel_;
    geomModel& tGeomModel_;

    virtualMeshInfo& vMeshInfo_;

    point contactCenter_;
    DynamicPointList edgeSubVolumesPoints_;

    public:


    // Constructors ------------------------------------------------

    virtualMesh
    (
        virtualMeshInfo& vMeshInfo,
        geomModel& cGeomModel_,
        geomModel& tGeomModel_
    );


    // Destructors -------------------------------------------------
    ~virtualMesh();

    bool detectFirstContactPoint();

    // Get contact subVolume
    bool detectFirstVolumeInContact
    (
        subVolume& sV,
        bool& startPointFound
    );

    scalar evaluateContact();

    void inspectSubVolume(
        subVolume& subVolume,
        scalar& contactVolume,
        vector& contactCenter,
        DynamicPointList& edgePoints
    );

    DynamicPointList& getEdgeSVPoints()
    {
        return edgeSubVolumesPoints_;
    }

    point& getContactCenter()
    {
        vMeshInfo_.startingPoint() = contactCenter_;
        return contactCenter_;
    }

    Tuple2<scalar,vector> get3DcontactNormalAndSurface(bool nonConvex);

    Tuple2<scalar,vector> get3DcontactNormalAndSurface(DynamicPointList edgeSubVolumesPoints);

    std::vector<subContact> findsubContacts(subVolume& sV);

    bool canCombineSubContacts(
        subContact& main,
        subContact& comp
    ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
