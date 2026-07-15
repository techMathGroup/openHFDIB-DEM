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
    lineInt.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
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
