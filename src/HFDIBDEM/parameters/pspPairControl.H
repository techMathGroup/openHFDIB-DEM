/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| |  \/ |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \| '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
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
    Static registry of body-pair states for the psp pass-through support.

    An application (e.g. the pass-through monitor of pspHFDIBDEMFoam) can
    register the frozen entry normal of a monitored pair with
    setForcedNormal; such pairs are flagged appControlled and the
    application keeps the responsibility for switching the pair's contact
    off after a pass-through.

    In addition, the registry drives an automatic per-pair state machine
    (updated from the contact model) that applies to ALL body pairs:
    while one body is fully embedded in another, the contact model keeps
    a resolved DEM contact alive by synthesising the contact variables
    from the frozen entry normal and the contained volume, and grows the
    synthesised overlap volume linearly with the (signed) penetration
    displacement since the embed onset times the maximum registered
    contact area of the episode - the embedded contact thus acts as a
    spring that stiffens with depth and relaxes on rebound. The
    synthesised contact carries a STATIC footprint: its area (the max
    geometric area of the episode), its shape (the polygon set of the
    last geometric contact patch) and its anchor (the patch centre of
    that contact) are frozen at the embed onset, and the reported
    contact area is the frozen one, so the derived pressure
    p = F/A_frozen carries the whole time dependence of the embedded
    load - the equivalent of a laser-shock-peening style surface
    pressure load. When the penetrator reaches the backside of the
    containing body with the movement direction unchanged, the pair is
    marked passedThrough and its contact force is switched off
    permanently (no exit-face repulsion).

Contributors
    Martin Isoz (2019-*)
\*---------------------------------------------------------------------------*/

#ifndef pspPairControl_H
#define pspPairControl_H

#include "vector.H"
#include "DynamicList.H"
#include "pointField.H"

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
            vector e_ = vector::zero; // unit entry normal (direction of the
                                     // penetrator's travel into the target);
                                     // the normal force on the penetrator
                                     // acts along -e_

            // --- automatic pass-through state (episode-scoped) -----------
            bool appControlled_ = false; // pair registered by an application
                                        // (setForcedNormal): the application
                                        // keeps the pass-through trigger
            vector lastNormal_ = vector::zero; // last standard contact normal
                                        // (points from the pair's second
                                        // member towards its first member)
            scalar maxArea_ = 0;      // max registered contact area of the
                                     // current contact episode
            vector lastCentre_ = vector::zero; // last geometric contact
                                     // patch centre (the static-footprint
                                     // anchor of a following embed)
            vector footprintCentre_ = vector::zero; // frozen footprint
                                     // centre of the embedded contact
            scalar footprintArea_ = 0; // frozen footprint area of the
                                     // embedded contact (static footprint:
                                     // the pressure carries the time
                                     // dependence, p = F/A_frozen)
            DynamicList<pointField> lastPolygons_; // polygon set of the
                                     // last geometric contact (the
                                     // footprint shape of a following
                                     // embed)
            DynamicList<pointField> footprintPolygons_; // frozen
                                     // footprint polygons of the embedded
                                     // contact (the last real patch,
                                     // translated to footprintCentre_)
            bool episodeActive_ = false; // a contact episode is ongoing
            bool embedded_ = false;   // penetrator fully inside target now
            bool passedThrough_ = false; // backside crossed with unchanged
                                     // direction: contact force off (sticky)
            vector relCoM0_ = vector::zero; // penetratorCoM - targetCoM at
                                     // the embed onset
            scalar sMax_ = 0;         // deepest signed penetration so far
        };

        // outcome of ending a contact episode (used for logging)
        enum episodeEndReason
        {
            noEpisode,        // pair had no active contact episode
            endPassedThrough, // episode already ended by backside penetration
            endWhileEmbedded, // contact lost while still embedded (body
                              // removed / switched off)
            endRepulsion      // standard end: bodies separated (overlap
                              // vanished without full immersion)
        };

        // register (or update) the frozen entry normal of a pair;
        // penetrator is the pair member that entered the other one. The
        // episode state of an existing record is preserved; the pair is
        // flagged appControlled (its pass-through trigger stays with the
        // application)
        static void setForcedNormal
        (
            label bodyA,
            label bodyB,
            label penetrator,
            const vector& e
        )
        {
            forcedNormal& fn(record(bodyA, bodyB));
            fn.penetrator_ = penetrator;
            fn.e_ = e/mag(e);
            fn.appControlled_ = true;
        }

        // false if no record is registered for the pair
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

        // --- automatic pass-through state machine -------------------------
        // create (or get) the record of a pair; ids are normalised to
        // (min, max)
        static forcedNormal& record(label bodyA, label bodyB)
        {
            const label a(min(bodyA, bodyB));
            const label b(max(bodyA, bodyB));
            forcedNormal& fn(forcedNormals()[std::make_pair(a, b)]);
            fn.bodyA_ = a;
            fn.bodyB_ = b;
            return fn;
        }

        // false if no record is registered for the pair (copy semantics as
        // getForcedNormal)
        static bool getRecord(label bodyA, label bodyB, forcedNormal& out)
        {
            return getForcedNormal(bodyA, bodyB, out);
        }

        // notify the registry that the geometric contact of the pair was
        // lost (the contact did not resolve in the current sub-step and the
        // contact record is being discarded). Returns the reason the
        // episode ended and resets the episode-scoped state; the
        // passedThrough_ flag is sticky (a pair that crossed the backside
        // stays force-free)
        static episodeEndReason endContactEpisode(label bodyA, label bodyB)
        {
            forcedNormal fn;
            if (!getForcedNormal(bodyA, bodyB, fn))
            {
                return noEpisode;
            }

            forcedNormal& r(record(bodyA, bodyB));

            if (!r.episodeActive_)
            {
                return noEpisode;
            }

            episodeEndReason reason(endRepulsion);
            if (r.passedThrough_)
            {
                reason = endPassedThrough;
            }
            else if (r.embedded_)
            {
                reason = endWhileEmbedded;
            }

            // reset the episode-scoped state; keep e_/penetrator_ (an
            // application-registered entry normal survives a repulsion
            // break-up so a monitored pair can re-embed) and the sticky
            // passedThrough_ flag
            r.lastNormal_ = vector::zero;
            r.maxArea_ = 0;
            r.lastCentre_ = vector::zero;
            r.footprintCentre_ = vector::zero;
            r.footprintArea_ = 0;
            r.lastPolygons_.clear();
            r.footprintPolygons_.clear();
            r.episodeActive_ = false;
            r.embedded_ = false;
            r.relCoM0_ = vector::zero;
            r.sMax_ = 0;

            return reason;
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
