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
    leastSquaresInt.C

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef leastSquaresInt_H
#define leastSquaresInt_H

#include "dictionary.H"
#include "fvCFD.H"
#include "triSurfaceMesh.H"
#include "ibInterpolation.H"
#include "leastSquaresIntInfo.H"

namespace Foam
{
    class leastSquaresInt
    : public ibInterpolation
    {
        private:
        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        const  fvMesh&   mesh_;

        // values corresponding to cell cell searching
        const  scalar   distFactor_;
        const  scalar   radiusFactor_;
        const  scalar   angleFactor_;
        const  scalar   maxCCRows_;

        void imposeDirichletCondition
        (
            leastSquaresIntInfo& intpInfo,
            volVectorField& Ui,
            vectorField& ibPointsVal,
            const Foam::fvMesh& mesh
        );

        void adjustPhi
        (
            leastSquaresIntInfo& intpInfo,
            surfaceScalarField& phi
        );

        public:

        // Constructors ------------------------------------------------
        leastSquaresInt(
            const fvMesh& mesh,
            scalar distFactor,
            scalar radiusFactor,
            scalar angleFactor,
            scalar maxCCRows
        );

        // Destructors -------------------------------------------------
        ~leastSquaresInt();

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
            return new leastSquaresIntInfo
            (
                mesh,
                gModel,
                distFactor_,
                radiusFactor_,
                angleFactor_,
                maxCCRows_
            );
        };
    };
}
#endif
