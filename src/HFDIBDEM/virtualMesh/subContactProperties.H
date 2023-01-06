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

#ifndef subContactProperties_H
#define subContactProperties_H

#include "subVolumeProperties.H"
#include "bbMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ibContactClass Declaration
\*---------------------------------------------------------------------------*/
class subContactProperties
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;
    typedef HashSet<vector,Hash<vector>> vectorHashSet;

    public:

    label contactID_;

    labelHashSet oldContactFaces_;

    vectorHashSet contactEdgePoints_;

    bool isContact_;

    autoPtr<bbMatrix> bbMatrix_;

    vector contactCenter_;

    DynamicLabelList boundaryCells_;

    DynamicLabelList boundaryFaces_;

    DynamicLabelList firstContactFaces_;

    DynamicVectorList firstContactPoints_;

    DynamicVectorList contactSubVolumesCOMRecord_;

    DynamicVectorList edgeSubVolumesCOMrecord_;

    DynamicVectorList searchedPoints_;

    point firstContactSubVolume_;

    scalar contactVolume_;

    scalar contactArea_;

    vector contactNormal_;

    label questionsAsked_;

    labelHashSet checkedVertexLabels_;

    string patchName_;

    void clear()
    {
        boundaryCells_.clear();
        boundaryFaces_.clear();
        contactSubVolumesCOMRecord_.clear();
        edgeSubVolumesCOMrecord_.clear();
        contactCenter_ = vector::zero;
        contactVolume_ =0.0;
        contactArea_ = 0.0;
        contactNormal_ = vector::zero;
        isContact_ = false;
        searchedPoints_.clear();
        checkedVertexLabels_.clear();
        contactEdgePoints_.clear();
        oldContactFaces_.clear();

    }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
