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

    openHFDIB-DEM is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

Description
    struct for interpolation point

SourceFiles

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef intPoint_H
#define intPoint_H

namespace Foam
{
    struct intPoint
    {
        // interpolation point
        point iPoint_ = vector::zero;

        // info about cell with interpolation point
        label iCell_ = 0;
        label iProc_ = -1;

        // info about origin surface cell
        label oProc_ = -1;
        label oLabel_ = -1;

        // is this point the last from origin cell
        bool last_ = true;

        // velocity in interpolation point
        vector iVel_ = vector::zero;

        intPoint():
        iPoint_(vector::zero),
        iCell_(0),
        iProc_(-1),
        oProc_(-1),
        oLabel_(-1),
        last_(true)
        {};

        intPoint(point iPoint, label iCell, label iProc, label oProc, label oLabel, bool last):
        iPoint_(iPoint),
        iCell_(iCell),
        iProc_(iProc),
        oProc_(oProc),
        oLabel_(oLabel),
        last_(last)
        {};

        ~intPoint(){};
    };
}
#endif
