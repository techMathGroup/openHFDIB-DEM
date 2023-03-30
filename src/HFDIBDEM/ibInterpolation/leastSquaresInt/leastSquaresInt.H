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
    class for addition and removal of immersedBodies

SourceFiles
    leastSquaresInt.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
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
