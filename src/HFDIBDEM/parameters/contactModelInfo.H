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
    Static class for storing virtual mesh level

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef contactModelInfo_H
#define contactModelInfo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class contactModelInfo Declaration
\*---------------------------------------------------------------------------*/


class contactModelInfo
{
    public:

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

    static void setIncreasedDamping(bool increasedDamping)
    {
        increasedDamping_ = increasedDamping;
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

    static bool getIncreasedDamping()
    {
        return increasedDamping_;
    }

    static scalar LcCoeff_;
    static bool useChenRotationalModel_;
    static bool useMindlinRotationalModel_;
    static bool increasedDamping_;



};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
