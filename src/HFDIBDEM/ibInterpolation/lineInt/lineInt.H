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
    class for addition and removal of immersedBodies

SourceFiles
    lineInt.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef lineInt_H
#define lineInt_H

#include "dictionary.H"
#include "fvCFD.H"
#include "triSurfaceMesh.H"
#include "interpolationCellPoint.H"
#include "interpolationCell.H"

#include "ibInterpolation.H"
#include "lineIntInfo.H"

namespace Foam
{
    class lineInt
    : public ibInterpolation
    {
        private:
        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        autoPtr<interpolation<vector>> interpV_;
        dictionary& interpDict_;

        void correctVelocity
        (
            lineIntInfo& intpInfo,
            volVectorField& Ui,
            vectorField& ibPointsVal,
            const Foam::fvMesh& mesh
        );

        void getCurVelocity
        (
            List<List<intPoint>>& intPoints
        );

        List<label> getIntOrder
        (
            List<List<intPoint>>& intPoints
        );

        public:

        // Constructors ------------------------------------------------
        lineInt(dictionary& interpDict);

        // Destructors -------------------------------------------------
        ~lineInt();

        void resetInterpolator(volVectorField& V) override
        {
            interpV_ = interpolation<vector>::New(interpDict_, V);
        }

        void ibInterpolate
        (
            interpolationInfo& intpInfo,
            volVectorField& Ui,
            vectorField ibPointsVal,
            const Foam::fvMesh& mesh
        ) override;

        interpolationInfo* createIntpInfo
        (
            const Foam::fvMesh& mesh,
            std::shared_ptr<geomModel>& gModel
        ) override
        {
            return new lineIntInfo
            (
                mesh,
                gModel
            );
        };
    };
}
#endif
