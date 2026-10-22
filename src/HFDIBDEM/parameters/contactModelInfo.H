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
    contactModel

Description
    Static class for storing contact model parameters and the
    shared tangential force law

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef contactModelInfo_H
#define contactModelInfo_H

#include "vector.H"
#include "demTimeStepInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class contactModelInfo Declaration
\*---------------------------------------------------------------------------*/


class contactModelInfo
{
    public:

    // dashpot prefactors of the dissipative branches: each
    // dashpot reads pref*zeta(e)*sqrt(m*k_branch); the branch
    // stiffnesses differ so the prefactors do too
    // (reducebeta_damping_consistency_20261021.md, Sec. 10-11)
    static scalar normalDashpotPrefactor()
    {
        return sqrt(3.0);
    }

    static scalar mindlinDashpotPrefactor()
    {
        return sqrt(3.0);
    }

    static scalar chenDashpotPrefactor()
    {
        return sqrt(5.0);
    }

    // shared tangential force law of the wall and particle
    // contact paths: the rotational model (mindlin/chen; an
    // unknown dict name defaults to mindlin through the
    // dict-reader) builds the spring increment, the coulomb
    // ceiling caps the applied force and the stored state
    // alike. FtPrev is the spring state, updated in place and
    // returned as the applied force; the call sites differ
    // only in where the inputs live
    static vector tangentialSpringForce
    (
        vector& FtPrev,
        const vector& relVeli,
        const vector& contactNormal,
        const scalar contactArea,
        const scalar Lc,
        const scalar aG,
        const scalar zetaE,
        const scalar reduceM,
        const scalar deltaT,
        const scalar FtCeil
    )
    {
        // project last Ft into a new direction
        vector FtLastP(FtPrev - (FtPrev & contactNormal)*contactNormal);
        // scale projected Ft to have same magnitude as FtLast
        vector FtLastS(mag(FtPrev)*(FtLastP/(mag(FtLastP)+SMALL)));

        // relative tangential velocity
        vector Vn((relVeli & contactNormal)*contactNormal);
        vector Vt(relVeli - Vn);

        // compute tangential force
        vector deltaFt(vector::zero);
        if (getUseMindlinRotationalModel())
        {
            // empirical tuning factor, see demTimeStepInfo
            scalar tangTune(demTimeStepInfo::tangTune_);
            scalar kT = tangTune*8.0*aG*(contactArea/(Lc+SMALL));
            deltaFt = kT*Vt*deltaT
                + mindlinDashpotPrefactor()*zetaE
                *sqrt(kT*reduceM)*Vt;
        }
        else if (getUseChenRotationalModel())
        {
            deltaFt = chenDashpotPrefactor()*zetaE
                *sqrt(aG*reduceM*Lc)*Vt;
            deltaFt += aG*Lc*Vt*deltaT;
        }

        // the spring can stretch by at most one ceiling per
        // sub-step: a fresh contact (or a slip reversal) builds
        // the tangential force up to the coulomb limit, never
        // across it within a single step
        if (mag(deltaFt) > FtCeil)
        {
            deltaFt *= FtCeil/(mag(deltaFt) + SMALL);
        }
        FtPrev = FtLastS - deltaFt;

        // coulomb cap feeds back into the stored state: the
        // spring stops stretching at the sliding ceiling, so
        // the applied force and the state can never disagree
        if (mag(FtPrev) > FtCeil)
        {
            if (mag(Vt) > vSlipMin_)
            {
                // sliding: kinetic friction opposes the slip;
                // the state is redirected with the applied
                // force, on re-stick the spring re-forms from
                // the sliding direction
                FtPrev = -FtCeil*Vt/(mag(Vt) + SMALL);
            }
            else
            {
                // below vSlipMin_ the slip direction is assumed
                // noise-dominated: keep the spring direction
                // (the cap firing from a shrinking ceiling must
                // not scramble the state)
                FtPrev *= FtCeil/(mag(FtPrev) + SMALL);
            }
        }

        return FtPrev;
    }

    static void setLcCoeff(scalar LcCoeff)
    {
        LcCoeff_ = LcCoeff;
    }

    static void setRotationModel(label modelID)
    {
        if (modelID == 0)
        {
            useChenRotationalModel_ = true;
            useMindlinRotationalModel_ = false;
        }
        else if (modelID == 1)
        {
            useChenRotationalModel_ = false;
            useMindlinRotationalModel_ = true;
        }
    }

    static scalar getLcCoeff()
    {
        return LcCoeff_;
    }

    static bool getUseChenRotationalModel()
    {
        return useChenRotationalModel_;
    }

    static bool getUseMindlinRotationalModel()
    {
        return useMindlinRotationalModel_;
    }

    static void setVSlipMin(scalar vSlipMin)
    {
        vSlipMin_ = vSlipMin;
    }

    static scalar getVSlipMin()
    {
        return vSlipMin_;
    }

    static scalar LcCoeff_;
    static bool useChenRotationalModel_;
    static bool useMindlinRotationalModel_;

    // slip-speed gate of the coulomb-cap direction redirect:
    // below it the slip direction is assumed noise-dominated
    static scalar vSlipMin_;



};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
