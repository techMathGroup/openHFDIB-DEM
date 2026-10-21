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
    Regime-aware normal force law for flyer-induced surface loading.

        p_str  = min( C_f*Yc_f , C_s*Yc_s )      strength branch (flyer cap)
        p_H    = impedance match                 shock branch (no strength cap)
        R      = p_H/p_str
        f      = smoothstep( log10(R) over [log10(lo), log10(hi)] )
        p_ceil = (1-f)*p_str + f*p_H

        F = min( k_eff*w + f*Z_eff*A*v , p_ceil*A )

    with w = V_overlap/A, k_eff = kappa*aY*sqrt(A) (half-space flat punch) or
    aY*A/Lc (the legacy layer spring), Z_eff = p_H/v (secant impedance of the
    Hugoniot match) and the elastic limit Z = Z1*Z2/(Z1+Z2).

    Every function here is a free function over plain doubles: no OpenFOAM
    types, no registry, no mesh. The registry that reads the dictionary
    (surfaceLoadingModel.H) only builds the parameter struct; the law itself
    is exercised directly by tests/surfaceLoading/.

Contributors
    Martin Isoz (2019-*)
\*---------------------------------------------------------------------------*/

#ifndef surfaceLoadingLaw_H
#define surfaceLoadingLaw_H

#include "materialDynamics.H"

#include <algorithm>
#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceLoadingLaw
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// the two materials of an active pair, oriented: flyer loads sample
struct pairParams
{
    materialDynamics flyer;
    materialDynamics sample;

    double aY = 0.0;          // series contact modulus [Pa]
    double reduceM = 0.0;     // reduced mass of the pair [kg]
    double reduceBeta = 0.0;  // standard DEM restitution damping [-]

    bool blendOn = true;      // false = strength branch only (f = 0)
    double ratioLo = 0.3;     // lower edge of the smoothstep in p_H/p_str
    double ratioHi = 3.0;     // upper edge
    bool secant = true;       // secant Hugoniot impedance (else elastic)
    bool halfSpace = true;    // kappa*aY*sqrt(A) (else aY*A/Lc)
    double kappa = 1.05;      // flat-punch stiffness factor
    bool capOn = true;        // apply the p_ceil ceiling to F
};

// per-sub-step outcome; the p* fields are what the application reports
struct result
{
    double F = 0.0;        // normal force [N]
    double pMean = 0.0;    // F/A [Pa]
    double pHug = 0.0;     // impedance-matched pressure [Pa]
    double pStr = 0.0;     // strength-branch cap [Pa]
    double pEdge = 0.0;    // sample flow pressure for the edge redistribution
    double blendF = 0.0;   // shock fraction f [-]
    double zEff = 0.0;     // effective impedance [kg/m2/s]
    double w = 0.0;        // mean overlap depth V/A [m]
    double kEff = 0.0;     // spring stiffness used [N/m]
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Hugoniot pressure of one material at particle velocity u
inline double materialPressure(const materialDynamics& m, const double u)
{
    if (u <= 0.0)
    {
        return 0.0;
    }

    const bool high(m.twoSegment() && u >= m.uMelt);

    return m.rho*((high ? m.c0High : m.c0) + (high ? m.sHigh : m.s)*u)*u;
}

// Interface pressure of two Hugoniots, P1 = P2 with u1 + u2 = v.
//
// A plain bisection over [0, v] is NOT valid across a melt segment: the
// Hugoniot jumps down in pressure at u = uMelt, so P(u) is discontinuous
// and f(u) = P_a(u) - P_b(v-u) is not monotone. Across the kapton melt
// kink the equations then have two roots (one on each segment) and the
// physical state is the FIRST one reached from u = 0, because the loading
// path enters at zero particle velocity and moves up the low segment until
// it either meets the other Hugoniot or runs off the end of the segment.
//
// The bracket is therefore found by a scan for the first sign change of
// f(u) and only then refined by bisection.
inline double impedanceMatch
(
    const materialDynamics& a,
    const materialDynamics& b,
    const double v
)
{
    if (v <= 0.0)
    {
        return 0.0;
    }

    const int nScan(256);

    double uPrev(0.0);
    double fPrev(materialPressure(a, 0.0) - materialPressure(b, v));

    if (fPrev >= 0.0)
    {
        return materialPressure(a, 0.0);
    }

    for (int i = 1; i <= nScan; ++i)
    {
        const double u(v*double(i)/double(nScan));

        const double f(materialPressure(a, u) - materialPressure(b, v - u));

        if (f >= 0.0)
        {
            double lo(uPrev);
            double hi(u);

            for (int k = 0; k < 60; ++k)
            {
                const double m(0.5*(lo + hi));

                if
                (
                    materialPressure(a, m) - materialPressure(b, v - m)
                  < 0.0
                )
                {
                    lo = m;
                }
                else
                {
                    hi = m;
                }
            }

            return materialPressure(a, 0.5*(lo + hi));
        }

        uPrev = u;
    }

    // unreachable: at u = v the second material carries no pressure and
    // the first one is positive
    return materialPressure(a, v);
}

// Elastic series impedance p = [Z1 Z2/(Z1 + Z2)] v
inline double elasticImpedance
(
    const materialDynamics& a,
    const materialDynamics& b
)
{
    const double za(a.rho*a.c0);
    const double zb(b.rho*b.c0);

    return za*zb/(za + zb);
}

inline double smoothstep(const double t)
{
    const double x(std::min(std::max(t, 0.0), 1.0));

    return x*x*(3.0 - 2.0*x);
}

// shock fraction f = smoothstep of log10(R) over [lo, hi]
inline double blendFactor
(
    const double pHug,
    const double pStr,
    const double lo,
    const double hi
)
{
    if (pStr <= 0.0)
    {
        return 1.0;
    }

    return smoothstep
    (
        (std::log10(pHug/pStr) - std::log10(lo))
       /(std::log10(hi) - std::log10(lo))
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Normal force of one sub-contact.
//   v  relative normal velocity (positive = approaching)
//   A  contact area
//   V  overlap volume
//   Lc characteristic length (legacy layer spring only)
inline result evaluate
(
    const pairParams& p,
    const double v,
    const double A,
    const double V,
    const double Lc
)
{
    result r;

    if (A <= 0.0)
    {
        return r;
    }

    r.w = V/A;

    r.pStr = std::min
    (
        p.flyer.capFactor*p.flyer.Yc,
        p.sample.capFactor*p.sample.Yc
    );
    r.pEdge = p.sample.capFactor*p.sample.Yc;

    r.zEff = p.secant
        ? (v > 0.0 ? impedanceMatch(p.flyer, p.sample, v)/v
                   : elasticImpedance(p.flyer, p.sample))
        : elasticImpedance(p.flyer, p.sample);
    r.pHug = r.zEff*v;

    r.blendF = p.blendOn
        ? blendFactor(r.pHug, r.pStr, p.ratioLo, p.ratioHi)
        : 0.0;

    const double pCeil((1.0 - r.blendF)*r.pStr + r.blendF*r.pHug);

    r.kEff = p.halfSpace
        ? p.kappa*p.aY*std::sqrt(A)
        : p.aY*A/(Lc + 1.0e-15);

    // legacy layer spring written as the original expression so that
    // switching the model off is bit-identical to the previous law
    const double Fspring
    (
        p.halfSpace ? r.kEff*r.w : p.aY*V/(Lc + 1.0e-15)
    );

    double F(Fspring + r.blendF*r.pHug*A);

    // standard DEM restitution damping, same expression and clamp as
    // the classic law, so 'blend none + cap off + LcLayer' reproduces it
    // bit for bit; the ceiling is applied after it
    const double damp
    (
        p.reduceBeta
       *std::sqrt(p.aY*p.reduceM*A/(Lc + 1.0e-15))
    );

    double Fd(damp*v);

    if (Fd < 0.0 && -Fd > F)
    {
        Fd = -F;
    }

    F += Fd;

    if (p.capOn && F > pCeil*A)
    {
        F = pCeil*A;
    }

    if (F < 0.0)
    {
        F = 0.0;
    }

    r.F = F;
    r.pMean = F/A;

    return r;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceLoadingLaw
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
