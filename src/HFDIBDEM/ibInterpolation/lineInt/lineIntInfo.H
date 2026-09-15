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
    class for interpolation

SourceFiles

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef lineIntInfo_H
#define lineIntInfo_H

#include "interpolationInfo.H"
#include "intPoint.H"

#include "processorPolyPatch.H"



namespace Foam
{
    class lineIntInfo
    : public interpolationInfo
    {
        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        private:

        autoPtr<List<List<intPoint>>> intPointsPtr_;                    // Points used for interp.
        autoPtr<List<List<vector>>> intVelocPtr_;                       // Velocity at int points

        intPoint  findIntPoint
        (
            intPoint& fromP,
            point& endP
        );

        void correctIntPoint
        (
            point ibPoint,
            intPoint& cPoint
        );

        vector getClosestPoint
        (
            vector ibPoint,
            intPoint& cPoint
        );

        label getFaceInDir
        (
            const intPoint& retPoint,
            const label prevFace
        );

        bool pointInCell
        (
            point pToCheck,
            label cToCheck
        );

        // Syncing functions
        void sendAndRecvIntPoints
        (
            List<DynamicList<point>>& ibPointsToSend,
            List<DynamicList<vector>>& ibNormalsToSend,
            List<DynamicList<intPoint>>& intPointsToSend,
            List<DynamicList<point>>& ibPointsRecv,
            List<DynamicList<vector>>& ibNormalsRecv,
            List<DynamicList<intPoint>>& intPointsRecv
        );

        void returnSolvedIntPoints
        (
            List<DynamicList<intPoint>>& intPointsToCont,
            List<DynamicList<intPoint>>& intPointsRecv,
            List<DynamicList<intPoint>>& intPointsSolved
        );

        public:

        // Constructors ------------------------------------------------
        lineIntInfo
        (
            const  fvMesh&   mesh,
            std::shared_ptr<geomModel>& gModel
        );

        // Destructors -------------------------------------------------
        ~lineIntInfo();

        void setIntpInfo() final override;

        void resetIntpInfo(label newSize)
        {
            ibPointsPtr_.reset(new List<point>(newSize, vector::zero));
            ibNormalsPtr_.reset(new List<point>(newSize, vector::zero));
            intPointsPtr_.reset(new List<List<intPoint>>(newSize));
            intVelocPtr_.reset(new List<List<vector>>(newSize));
        }

        List<List<intPoint>>& getIntPoints()
        {
            return intPointsPtr_();
        }

        const List<List<intPoint>>& getIntPoints() const
        {
            return intPointsPtr_();
        }

        List<List<vector>>& getIntVeloc()
        {
            return intVelocPtr_();
        }

        const List<List<vector>>& getIntVeloc() const
        {
            return intVelocPtr_();
        }
    };
}
#endif
