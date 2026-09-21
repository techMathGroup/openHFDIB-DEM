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

#ifndef interpolationInfo_H
#define interpolationInfo_H

#include "dictionary.H"
#include "fvCFD.H"
// #include "fvMesh.H"
#include "triSurfaceMesh.H"
#include "geomModel.H"

namespace Foam
{
    class interpolationInfo
    {
        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        protected:

        // mesh reference
        const Foam::fvMesh& mesh_;
        // geomModel reference
        std::shared_ptr<geomModel>& geomModel_;

        autoPtr<List<point>> ibPointsPtr_;                              //cellIbPoint mapping
        autoPtr<List<vector>> ibNormalsPtr_;                            //cellIbNormal mapping
        autoPtr<labelList> ibCellLabelsPtr_;                            //ibPointCell mapping
        autoPtr<HashTable<label, label, Hash<label>>> cellToIbPointPtr_;
        scalar charCellSize_;

        protected:

        //-- reset the base arrays shared by all interpolation types
        void resetIntpInfoBase(label newSize)
        {
            ibPointsPtr_.reset(new List<point>(newSize, vector::zero));
            ibNormalsPtr_.reset(new List<vector>(newSize, vector::zero));
            ibCellLabelsPtr_.reset(new labelList(newSize, -1));
            cellToIbPointPtr_.reset(
                new HashTable<label, label, Hash<label>>(2*newSize));
        }

        //-- record the surfCell-ibPoint correspondence
        void setIbCellLabels(const DynamicLabelList& cSurfCells)
        {
            labelList& ibCellLabels = ibCellLabelsPtr_();
            HashTable<label, label, Hash<label>>& cellToIbPoint
                = cellToIbPointPtr_();

            forAll(cSurfCells, cellI)
            {
                ibCellLabels[cellI] = cSurfCells[cellI];

                if (!cellToIbPoint.insert(cSurfCells[cellI], cellI))
                {
                    FatalErrorInFunction
                        << "surface cell " << cSurfCells[cellI]
                        << " appears twice in the surface cell list"
                        << abort(FatalError);
                }
            }
        }

        public:

        // Constructors ------------------------------------------------
        interpolationInfo
        (
            const Foam::fvMesh& mesh,
            std::shared_ptr<geomModel>& gModel
        );

        // Destructors -------------------------------------------------
        virtual ~interpolationInfo();

        List<point>& getIbPoints()
        {
            //~ if(ibPointsPtr_->empty())
            //~ {
            //~     setIntpInfo();
            //~ }
            return ibPointsPtr_();
        }

        List<vector>& getIbNormals()
        {
            return ibNormalsPtr_();
        }

        // get ibPoint index from the mesh cell label
        label findIbPoint(label cellI) const
        {
            if (!cellToIbPointPtr_.valid())
            {
                FatalErrorInFunction
                    << "interpolation info was cleared, setIntpInfo needed"
                    << abort(FatalError);
            }

            if (!cellToIbPointPtr_().found(cellI))
            {
                FatalErrorInFunction
                    << "no ibPoint for cell " << cellI
                    << " - surface cells changed since setIntpInfo"
                    << abort(FatalError);
            }

            return cellToIbPointPtr_()[cellI];
        }

        virtual void setIntpInfo() = 0;

        const DynamicLabelList& getSurfCells() const
        {
            return geomModel_->getSurfaceCellList()[Pstream::myProcNo()];
        }

        virtual void clearIntpInfo()
        {
            ibPointsPtr_.clear();
            ibNormalsPtr_.clear();
            ibCellLabelsPtr_.clear();
            cellToIbPointPtr_.clear();
        }

        void setCharCellSize(scalar charCellSize)
        {
            charCellSize_ = charCellSize;
        }
    };
}
#endif
