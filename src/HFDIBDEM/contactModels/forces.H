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
    class for storage of contact variables

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef forces_H
#define forces_H

#include "dictionary.H"
#include "fvCFD.H"

#include "geomModel.H"
#include "ibContactVars.H"
#include "wallContactVars.H"
#include "materialInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class forces Declaration
\*---------------------------------------------------------------------------*/
struct forces{
    vector F = vector::zero;
    vector T = vector::zero;

    forces()
    {}

    forces(const forces& rhs)
    {
        F = rhs.F;
        T = rhs.T;
    }

    forces(const vector& F, const vector& T):
        F(F),
        T(T)
    {}

    void operator=(const forces& rhs)
    {
        F = rhs.F;
        T = rhs.T;
    }

    void operator+=(const forces& rhs)
    {
        F += rhs.F;
        T += rhs.T;
    }

    void operator-=(const forces& rhs)
    {
        F -= rhs.F;
        T -= rhs.T;
    }
    
    void operator*=(const scalar& coef)
    {
        F *= coef;
        T *= coef;
    }

    void clear()
    {
        F = vector::zero;
        T = vector::zero;
    }
};

// symmetric multiplication by scalar
inline forces operator*(const forces& f, const scalar& coef) 
{
    return forces(f.F * coef, f.T * coef);
}

inline forces operator*(const scalar& coef, const forces& f) 
{
    return f * coef;  // use the previously defined operator
}

// addition operator
inline forces operator+(const forces& f, const forces& g) 
{
    return forces(f.F + g.F, f.T + g.T);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
