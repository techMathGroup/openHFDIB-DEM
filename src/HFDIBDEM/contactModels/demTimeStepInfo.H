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
    pure scalar helpers for the adaptive DEM time-step estimation.
    all formulas are pair-local; wall contacts use reduceM equal to the
    particle mass and the wall analogue of the reduced moduli.

    the estimated stability step of a contact pair is
        dtPair = stepDEM * 2 * sqrt(reduceM / kEff)
    with the governing stiffness kEff = max(kN, kT) taken from the
    contact laws as implemented in prtSubContactInfo / wallSubContactInfo.

SourceFiles
    demTimeStepInfo.H

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef demTimeStepInfo_H
#define demTimeStepInfo_H

#include "scalar.H"
#include "mathematicalConstants.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class demTimeStepInfo Declaration
\*---------------------------------------------------------------------------*/

class demTimeStepInfo
{
    public:

    // Note (MI): tangTune_ = 200 was an empirical tuning factor (cc92e39; Ondra),
    // not from mindlin theory. if we want to plate with this, this
    // is the source of truth for all the models
    static constexpr scalar tangTune_ = 1.0;

    // stiffness of the implemented normal law (getFNe form):
    // FNe = aY*contactVolume/Lc, hence kN = aY*contactVolume/Lc
    static scalar contactKN
    (
        const scalar aY,
        const scalar contactVolume,
        const scalar Lc
    )
    {
        return aY*contactVolume/(Lc + SMALL);
    }

    // stiffness of the implemented tangential law (mindlin branch of
    // getFt): kT = empiricalConstant*8*aG*contactArea/Lc               <- empiricalConstant = tangTune_
    static scalar contactKT
    (
        const scalar aG,
        const scalar contactArea,
        const scalar Lc
    )
    {
        return tangTune_*8.0*aG*contactArea/(Lc + SMALL);
    }

    // governing stiffness used for step sizing; tangentialFactor scales
    // kT only (1 = full conservatism)
    static scalar effK
    (
        const scalar aY,
        const scalar contactVolume,
        const scalar aG,
        const scalar contactArea,
        const scalar Lc,
        const scalar tangentialFactor
    )
    {
        return max
        (
            contactKN(aY, contactVolume, Lc),
            tangentialFactor*contactKT(aG, contactArea, Lc)
        );
    }

    // stability step of a spring-mass pair (undamped, explicit leapfrog):
    // dt < 2*sqrt(reduceM/k). damping sits inside the stepDEM safety
    // factor (f(xi) = 1 by agreement).
    static scalar pairDtCrit(const scalar kEff, const scalar reduceM)
    {
        return 2.0*Foam::sqrt(reduceM/(kEff + VSMALL));
    }

    // hertz max normal overlap from the elastic energy balance:
    // deltaMax = (15*reduceM*v^2/(16*aY*sqrt(R)))^(2/5)
    // used only by the velocityAware area estimate and the hertz
    // contact-duration diagnostic
    static scalar hertzDeltaMax
    (
        const scalar reduceM,
        const scalar aY,
        const scalar R,
        const scalar v
    )
    {
        return Foam::pow
        (
            15.0*reduceM*sqr(v)/(16.0*aY*Foam::sqrt(R + SMALL)),
            0.4
        );
    }

    // first-touch contact area estimate for pre-contact pairs:
    // - conservative: aEst = areaCoeff*h^2 (h = min pair charCellSize)
    // - velocityAware: min of the above and the hertz peak area
    //   pi*R*deltaMax. a larger R gives a smaller deltaMax and a
    //   larger dt, so R must be the bodies' equivalent radius
    //   (mesh-based R would be the un-conservative direction here)
    static scalar aEst
    (
        const scalar h,
        const scalar R,
        const scalar reduceM,
        const scalar aY,
        const scalar v,
        const scalar areaCoeff,
        const bool velocityAware
    )
    {
        scalar aGeom(areaCoeff*sqr(h));

        if (!velocityAware)
        {
            return aGeom;
        }

        return min
        (
            aGeom,
            constant::mathematical::pi*R
                *hertzDeltaMax(reduceM, aY, R, v)
        );
    }

    // hertz contact duration (diagnostic only, never used for sizing):
    // tC = 2.87*(reduceM^2/(R*aY^2*v))^(1/5)
    static scalar hertzContactDt
    (
        const scalar reduceM,
        const scalar aY,
        const scalar R,
        const scalar v
    )
    {
        return 2.87*Foam::pow
        (
            sqr(reduceM)/(R*Foam::sqr(aY)*(v + SMALL)),
            0.2
        );
    }

    // rayleigh surface-wave step (diagnostic only, never used for
    // sizing): dtR = pi*R*sqrt(rho/aG)/(0.1631*nu + 0.8766)
    static scalar rayleighDt
    (
        const scalar rho,
        const scalar aG,
        const scalar nu,
        const scalar R
    )
    {
        return constant::mathematical::pi*R*Foam::sqrt(rho/(aG + SMALL))
            /(0.1631*nu + 0.8766);
    }

    // equivalent sphere radius of a body from its volume
    static scalar equivRadius(const scalar volume)
    {
        return Foam::cbrt
        (
            3.0*volume/(4.0*constant::mathematical::pi) + SMALL
        );
    }

    // reduced (harmonic) radius of a pair from the body radii
    static scalar harmonicRadius(const scalar rC, const scalar rT)
    {
        return rC*rT/(rC + rT + SMALL);
    }

    // smallest principal moment of inertia: the rotation-invariant
    // conservative Ieff.
    // Foam::eigenValues clips near-degenerate roots of the
    // characteristic cubic to the double root (verified:
    // diag(a,a,a) + c on (x,z) returns a instead of a-c for c
    // below ~1e-3*a) - exactly the two-equal-principal-moment
    // shapes (spheres, cylinders, discs) - and the smallest
    // eigenvalue then comes back too large, the unconservative
    // direction for the stability step. the Gershgorin row bound
    // is a guaranteed lower bound on every eigenvalue, so the
    // smaller of the two is returned: it can never exceed the
    // true smallest eigenvalue. the bound is used only when
    // positive - a non-positive bound carries no information
    // (frames with large off-diagonals) and would poison the
    // sqrt in rotationalDtCrit
    static scalar minEigenvalue(const symmTensor& I)
    {
        const vector ev(eigenValues(I));
        const scalar evMin(min(min(ev.x(), ev.y()), ev.z()));

        // gershgorin: every eigenvalue lies within
        // [T_ii - sum|T_ij|, T_ii + sum|T_ij|] of some diagonal
        const scalar gershMin
        (
            min
            (
                min
                (
                    I.xx() - mag(I.xy()) - mag(I.xz()),
                    I.yy() - mag(I.xy()) - mag(I.yz())
                ),
                I.zz() - mag(I.xz()) - mag(I.yz())
            )
        );

        if (gershMin > SMALL*mag(evMin))
        {
            return min(evMin, gershMin);
        }

        return evMin;
    }

    // rotational stability step of a tangential contact: the
    // tangential force kT*x on a lever arm r acts on the body
    // inertia about the resulting torque axis
    static scalar rotationalDtCrit
    (
        const scalar kT,
        const scalar leverArm,
        const scalar Ieff
    )
    {
        const scalar kEff(kT*sqr(leverArm));

        if (kEff < SMALL) return GREAT;                                 // unconstrained rotation

        return 2.0*Foam::sqrt(Ieff/kEff);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
