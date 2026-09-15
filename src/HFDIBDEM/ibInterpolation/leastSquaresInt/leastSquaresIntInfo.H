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
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/

#ifndef leastSquaresIntInfo_H
#define leastSquaresIntInfo_H

#include "interpolationInfo.H"

namespace Foam
{
    class leastSquaresIntInfo
    : public interpolationInfo
    {
        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        private:

        autoPtr<labelListList> cellCellsPtr_; // cellCell mapping
        PtrList<scalarRectangularMatrix> invDirichletMatrices_;
        autoPtr<DynamicLabelList> ibFacesPtr_;
        autoPtr<DynamicList<bool>> ibFaceFlipsPtr_;

        // values corresponding to cell cell searching
        const  scalar   distFactor_;
        const  scalar   radiusFactor_;
        const  scalar   angleFactor_;
        const  scalar   maxCCRows_;

        void findCellCells
        (
            const label cellId,
            labelList& cellCells,
            scalar& centerMeanDist
        );

        void getInvDirichletMatrix
        (
            const DynamicLabelList& ibCells
        );

        void getIbCellsFaces
        (
            const DynamicLabelList& ibCells
        );

        scalarSquareMatrix LUinvert
        (
            scalarSquareMatrix& matrix
        );

        // void setIntpInfo
        // (
        //     const List<DynamicLabelList>& surfCells
        // );

        public:

        // Constructors ------------------------------------------------
        leastSquaresIntInfo
        (
            const  fvMesh&   mesh,
            std::shared_ptr<geomModel>& gModel,
            scalar distFactor,
            scalar radiusFactor,
            scalar angleFactor,
            scalar maxCCRows
        );

        // Destructors -------------------------------------------------
        ~leastSquaresIntInfo();

        void resetIntpInfo(label newSize)
        {
            ibPointsPtr_.reset(new List<point>(newSize, vector::zero));
            ibNormalsPtr_.reset(new List<vector>(newSize, vector::zero));
            cellCellsPtr_.reset(new labelListList(newSize));
            ibFacesPtr_.reset(new DynamicLabelList());
            ibFaceFlipsPtr_.reset(new DynamicList<bool>());
        }

        void resetInvDirMat(label newSize)
        {
            invDirichletMatrices_.setSize(newSize);
        }

        void setIntpInfo() final override;

        List<vector>& getIbNormals()
        {
            return ibNormalsPtr_();
        }

        const List<vector>& getIbNormals() const
        {
            return ibNormalsPtr_();
        }

        labelListList& getCellCells()
        {
            return cellCellsPtr_();
        }

        const labelListList& getCellCells() const
        {
            return cellCellsPtr_();
        }

        PtrList<scalarRectangularMatrix>& getInvDirMats()
        {
            return invDirichletMatrices_;
        }

        const PtrList<scalarRectangularMatrix>& getInvDirMats() const
        {
            return invDirichletMatrices_;
        }

        DynamicLabelList& getIbFaces()
        {
            return ibFacesPtr_();
        }

        const DynamicLabelList& getIbFaces() const
        {
            return ibFacesPtr_();
        }

        DynamicList<bool>& getIbFacesFlips()
        {
            return ibFaceFlipsPtr_();
        }

        const DynamicList<bool>& getIbFacesFlips() const
        {
            return ibFaceFlipsPtr_();
        }
    };
}
#endif
