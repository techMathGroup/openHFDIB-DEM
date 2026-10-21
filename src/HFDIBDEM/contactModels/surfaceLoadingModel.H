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
    Static registry of the surface-loading contact model: which body pairs use
    it and with which settings. Read once from DEM/surfaceLoading of
    constant/HFDIBDEMDict; the law itself lives in surfaceLoadingLaw.H.

        surfaceLoading
        {
            pairs           ((0 3));     // (flyer sample), ids as in bodiesInfo
            blend           pressureRatio; // pressureRatio | none
            ratioLo         0.3;
            ratioHi         3.0;
            impedance       secantHugoniot; // secantHugoniot | elasticImpedance
            stiffness       halfSpace;      // halfSpace | LcLayer
            stiffnessKappa  1.05;
            strengthCap     on;             // on | off
            edgeReg         0.05;           // d_reg/R_foot of the shape law
            writeInterval   5e-10;          // sub-ns loading record cadence [s]
            writeSnapshots  false;          // sparse per-cell field dumps
        }

    Absence of the block (or an empty pair list) leaves every contact on the
    standard law, so existing cases are unaffected.

Contributors
    Martin Isoz (2019-*)
\*---------------------------------------------------------------------------*/

#ifndef surfaceLoadingModel_H
#define surfaceLoadingModel_H

#include "dictionary.H"
#include "List.H"
#include "Tuple2.H"

#include <map>
#include <set>
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class surfaceLoadingModel Declaration
\*---------------------------------------------------------------------------*/

class surfaceLoadingModel
{
    public:

        struct settings
        {
            bool active = false;      // any pair registered
            bool blendOn = true;      // pressureRatio (true) or none (false)
            bool secant = true;       // secantHugoniot or elasticImpedance
            bool halfSpace = true;    // kappa*aY*sqrt(A) or aY*A/Lc
            bool capOn = true;        // strength ceiling applied
            scalar ratioLo = 0.3;
            scalar ratioHi = 3.0;
            scalar kappa = 1.05;
            scalar edgeReg = 0.05;
            scalar writeInterval = 5e-10;
            bool writeSnapshots = false;

            // pair key (min,max) -> flyer body id
            std::map<std::pair<label,label>, label> flyerOf;

            // pairs whose first contact has already been reported
            std::set<std::pair<label,label>> reported;
        };

        static settings& get()
        {
            static settings s;
            return s;
        }

        static bool isActive(const label a, const label b)
        {
            return get().flyerOf.find(key(a, b)) != get().flyerOf.end();
        }

        // the registered flyer of a pair; -1 when the pair is not active
        static label flyerOf(const label a, const label b)
        {
            const auto it(get().flyerOf.find(key(a, b)));

            return it == get().flyerOf.end() ? label(-1) : it->second;
        }

        // true only on the first call for a pair (used for the one-shot
        // first-contact report)
        static bool markReported(const label a, const label b)
        {
            return get().reported.insert(key(a, b)).second;
        }

        // read DEM/surfaceLoading; call once, after the materials are loaded
        static void read(const dictionary& demDic)
        {
            settings& s(get());

            if (!demDic.found("surfaceLoading"))
            {
                s.active = false;
                return;
            }

            const dictionary& d(demDic.subDict("surfaceLoading"));

            const word blend(d.lookupOrDefault<word>("blend", "pressureRatio"));
            if (blend == "pressureRatio")
            {
                s.blendOn = true;
            }
            else if (blend == "none")
            {
                s.blendOn = false;
            }
            else
            {
                FatalErrorInFunction
                    << "surfaceLoading: unknown blend '" << blend << "'"
                    << " (pressureRatio | none)" << exit(FatalError);
            }

            const word imp
            (
                d.lookupOrDefault<word>("impedance", "secantHugoniot")
            );
            if (imp == "secantHugoniot")
            {
                s.secant = true;
            }
            else if (imp == "elasticImpedance")
            {
                s.secant = false;
            }
            else
            {
                FatalErrorInFunction
                    << "surfaceLoading: unknown impedance '" << imp << "'"
                    << " (secantHugoniot | elasticImpedance)"
                    << exit(FatalError);
            }

            const word stiff(d.lookupOrDefault<word>("stiffness", "halfSpace"));
            if (stiff == "halfSpace")
            {
                s.halfSpace = true;
            }
            else if (stiff == "LcLayer")
            {
                s.halfSpace = false;
            }
            else
            {
                FatalErrorInFunction
                    << "surfaceLoading: unknown stiffness '" << stiff << "'"
                    << " (halfSpace | LcLayer)" << exit(FatalError);
            }

            const word cap(d.lookupOrDefault<word>("strengthCap", "on"));
            if (cap == "on" || cap == "true" || cap == "yes")
            {
                s.capOn = true;
            }
            else if (cap == "off" || cap == "false" || cap == "no")
            {
                s.capOn = false;
            }
            else
            {
                FatalErrorInFunction
                    << "surfaceLoading: unknown strengthCap '" << cap << "'"
                    << " (on | off)" << exit(FatalError);
            }

            s.ratioLo = d.lookupOrDefault<scalar>("ratioLo", 0.3);
            s.ratioHi = d.lookupOrDefault<scalar>("ratioHi", 3.0);
            s.kappa = d.lookupOrDefault<scalar>("stiffnessKappa", 1.05);
            s.edgeReg = d.lookupOrDefault<scalar>("edgeReg", 0.05);
            s.writeInterval =
                d.lookupOrDefault<scalar>("writeInterval", 5e-10);
            s.writeSnapshots =
                d.lookupOrDefault<bool>("writeSnapshots", false);

            if (s.ratioLo <= 0 || s.ratioHi <= s.ratioLo)
            {
                FatalErrorInFunction
                    << "surfaceLoading: need 0 < ratioLo < ratioHi, got "
                    << s.ratioLo << " and " << s.ratioHi
                    << exit(FatalError);
            }

            s.flyerOf.clear();

            if (d.found("pairs"))
            {
                const List<Tuple2<label,label>> pairs
                (
                    d.lookup("pairs")
                );

                forAll(pairs, i)
                {
                    const label f(pairs[i].first());
                    const label t(pairs[i].second());

                    if (f == t)
                    {
                        FatalErrorInFunction
                            << "surfaceLoading: pair " << f << "-" << t
                            << " has identical members" << exit(FatalError);
                    }

                    s.flyerOf[key(f, t)] = f;
                }
            }

            s.active = !s.flyerOf.empty();

            Info<< "surface-loading model: " << s.flyerOf.size()
                << " pair(s), blend "
                << (s.blendOn ? "pressureRatio" : "none")
                << ", impedance "
                << (s.secant ? "secantHugoniot" : "elasticImpedance")
                << ", stiffness "
                << (s.halfSpace ? "halfSpace" : "LcLayer")
                << ", strengthCap " << (s.capOn ? "on" : "off")
                << endl;
        }

    private:

        static std::pair<label,label> key(const label a, const label b)
        {
            return std::make_pair(min(a, b), max(a, b));
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
