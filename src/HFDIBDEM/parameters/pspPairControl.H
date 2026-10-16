/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  | |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___/  |_____|_|  |_|
      | |                     H ybrid F ictitious D omain - I mmersed B oundary
      |_|                                        and D iscrete E leement M ethod
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
    Static registry of body pairs whose contact is externally controlled
    (psp pass-through support). An application (e.g. the pass-through
    monitor of pspHFDIBDEMFoam) registers the frozen entry normal of a
    monitored pair; the contact model then keeps the pair in a resolved
    DEM contact while one body is fully embedded in the other (the
    geometric contact hull degenerates there) by synthesising the contact
    variables from the frozen normal and the contained volume.

Contributors
    Martin Isoz (2019-*)
\*---------------------------------------------------------------------------*/

#ifndef pspPairControl_H
#define pspPairControl_H

#include "vector.H"

#include <map>
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class pspPairControl Declaration
\*---------------------------------------------------------------------------*/

class pspPairControl
{
    public:

        struct forcedNormal
        {
            label bodyA_ = -1;       // pair ids, bodyA < bodyB
            label bodyB_ = -1;
            label penetrator_ = -1;  // pair member that penetrated the other
            vector e_ = vector::zero; // unit entry normal; the normal force on
                                     // the penetrator acts along +e_
        };

        // register (or overwrite) the frozen entry normal of a pair;
        // penetrator is the pair member that entered the other one
        static void setForcedNormal
        (
            label bodyA,
            label bodyB,
            label penetrator,
            const vector& e
        )
        {
            forcedNormal fn;
            fn.bodyA_ = min(bodyA, bodyB);
            fn.bodyB_ = max(bodyA, bodyB);
            fn.penetrator_ = penetrator;
            fn.e_ = e/mag(e);
            forcedNormals()[std::make_pair(fn.bodyA_, fn.bodyB_)] = fn;
        }

        // false if no frozen normal is registered for the pair
        static bool getForcedNormal(label bodyA, label bodyB, forcedNormal& out)
        {
            const auto it
            (
                forcedNormals().find
                (
                    std::make_pair(min(bodyA, bodyB), max(bodyA, bodyB))
                )
            );
            if (it == forcedNormals().end())
            {
                return false;
            }
            out = it->second;
            return true;
        }

        static bool isForcedPair(label bodyA, label bodyB)
        {
            forcedNormal fn;
            return getForcedNormal(bodyA, bodyB, fn);
        }

    private:

        static std::map<std::pair<label, label>, forcedNormal>& forcedNormals()
        {
            static std::map<std::pair<label, label>, forcedNormal> fn;
            return fn;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
