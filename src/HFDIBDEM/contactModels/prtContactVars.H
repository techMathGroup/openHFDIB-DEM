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
    Foam

Description
    class for storage of particle-particle contact variables

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef prtContactVars_H
#define prtContactVars_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DynamicList.H"
#include "pointField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class prtContactVars Declaration
\*---------------------------------------------------------------------------*/

struct prtContactVars
{
    vector contactCenter_ = vector::zero;
    scalar contactVolume_ = 0.0f;
    vector contactNormal_ = vector::zero;
    scalar contactArea_ = 0.0f;

    // psp: true when this contact was synthesised from the frozen entry
    // normal of an embedded pass-through pair rather than evaluated from
    // the geometric contact hull
    bool synthesized_ = false;

    // contact plane footprint(s): every pointField is a triangle fan,
    // element 0 is the fan centre and the remaining points follow
    // the contact patch boundary. Empty when no footprint is available
    DynamicList<pointField> contactPolygons_;

    // surface-loading coefficients of the current sub-step (published by
    // prtSubContactInfo::getFNSl for the application's shape and record
    // writers; all zero for pairs on the standard law)
    bool surfaceLoading_ = false;
    scalar slPMean_ = 0;    // F/A [Pa]
    scalar slPStr_ = 0;     // strength branch cap [Pa]
    scalar slPEdge_ = 0;    // sample flow pressure for the edge cap [Pa]
    scalar slPHug_ = 0;     // impedance-matched pressure [Pa]
    scalar slBlendF_ = 0;   // shock fraction f [-]
    scalar slZEff_ = 0;     // effective impedance [kg/m2/s]
    scalar slW_ = 0;        // mean overlap depth [m]
    scalar slVn_ = 0;       // relative normal velocity [m/s]
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
