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
    Optional observer of the DEM sub-stepping in openHFDIBDEM::updateDEM.
    All hooks are empty by default. An application registers an observer
    through openHFDIBDEM::setDEMObserver() to read the contact state
    (e.g. for contact pressure logging) or to switch off pair contacts,
    without copying updateDEM. Without a registered observer the solver
    behaves exactly as before.

SourceFiles
    openHFDIBDEM.C
\*---------------------------------------------------------------------------*/

#ifndef demObserver_H
#define demObserver_H

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class prtSubContactInfo;

class demObserver
{
    public:

    virtual ~demObserver() {}

    // Called on the rank that solved the body-body (sub)contact, right
    // after a resolved contact; deltaTDEM is the DEM sub-step.
    virtual void prtContactSolved
    (
        prtSubContactInfo& subCInfo,
        const scalar deltaTDEM
    )
    {}

    // Called on all ranks once per DEM sub-step, after the contact forces
    // were applied and before the closing half-step velocity update
    // (collective calls are allowed). Ignoring a body pair here removes
    // its contact from the next sub-step on.
    virtual void demSubStepContactsSolved(const scalar deltaTDEM)
    {}
};

} // End namespace Foam

#endif

// ************************************************************************* //
