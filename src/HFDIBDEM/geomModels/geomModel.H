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
    geomModel.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef geomModel_H
#define geomModel_H

#include <vector>
#include <memory>

#include "dictionary.H"
#include "fvCFD.H"
#include "triSurfaceMesh.H"
#include "triangleFuncs.H"

#include "processorPolyPatch.H"

#include "contactType.H"
#include "subVolume.H"
#include "verletPoint.H"

#include "outputHFDIBDEM.H"
using namespace outputHFDIBDEM;

namespace Foam
{
    // custom field filtration
    template<class Type>
    Foam::tmp<Foam::Field<Type>>
    filterField
    (
        const Field<Type>& field,
        const labelList& cellIDs
    )
    {
        if (isNull(cellIDs))
        {
            return field;
        }
        else
        {
            return tmp<Field<Type>>(new Field<Type>(field, cellIDs));
        }
    }

    class geomModel
    {

        private:

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

        //contact type based on geometry
        const contactType contactType_;

        protected:

        const  fvMesh&   mesh_;

        // body surface cells
        List<DynamicLabelList>  surfCells_;

        // body internal cells
        List<DynamicLabelList>  intCells_;

        // body halo cells
        List<DynamicLabelList>  haloCells_;

        // label mapping from surface to halo cells
        List<List<DynamicLabelList>>  surfToHaloLabels_;

        //list to find the owner of the IB
        DynamicLabelList ibPartialVolume_;
        label owner_;

        label cellToStartInCreateIB_;
        scalar  thrSurf_;

        scalar  intSpan_;
        bool    sdBasedLambda_;
        boundBox curMeshBounds_;
        word overlayLambda_;

        //cachedNeighbours_ may cost problems in OF.com
        // autoPtr<HashTable<const labelList&, label, Hash<label>>> cachedNeighbours_;
        autoPtr<HashTable<labelList, label, Hash<label>>> cachedNeighbours_;
        // HashTable<labelList, label, Hash<label>> cachedNeighbours_;
        autoPtr<HashTable<DynamicList<labelList>, label, Hash<label>>> procNeighbours_;

        //body mass
        scalar M_;
        scalar M0_;                                                     //initial mass (at addition)
        label nCells_;

        //location of center of mass
        vector CoM_;

        //moment of intertia tensor
        symmTensor  I_;

        //Boundbox of the particle
        std::shared_ptr<boundBox> bBox_;

        //characteristic diameter of the body
        scalar  dC_;

        //body density
        dimensionedScalar rhoS_;

        label   refineBuffers_; 

        public:

        // Constructors ------------------------------------------------
        geomModel
        (
            const  fvMesh&   mesh,
            const contactType cType,
            scalar  thrSurf = 1
        );

        // Destructors -------------------------------------------------
        virtual ~geomModel();

        virtual void createImmersedBody
        (
            volScalarField& body,
            Field<label>& octreeField,
            List<labelList>& cellPoints
        )
        {
            InfoH << iB_Info << "1-- geomModelMessage-- "
                << "notImplemented"<<endl;
        }

        virtual std::shared_ptr<geomModel> getCopy() = 0;

        virtual List<DynamicLabelList>& getInternalCellList()
        {
            return intCells_;
        };

        virtual List<DynamicLabelList>& getSurfaceCellList()
        {
            return surfCells_;
        };

        virtual List<DynamicLabelList>& getHaloCellList()
        {
            return haloCells_;
        };

        virtual List<List<DynamicLabelList>>& getSurfToHaloLabels()
        {
            return surfToHaloLabels_;
        };

        virtual void getReferencedLists(
            List<DynamicLabelList>& intLists,
            List<DynamicLabelList>& surfLists,
            List<DynamicLabelList>& haloLists,
            DynamicVectorList& referenceCoM
        )
        {
            intLists.resize(1, intCells_[Pstream::myProcNo()]);
            surfLists.resize(1, surfCells_[Pstream::myProcNo()]);
            haloLists.resize(1, haloCells_[Pstream::myProcNo()]);
            referenceCoM.resize(1, getCoM());
        };

        virtual void getReferencedInternalCellList(
            List<DynamicLabelList>& intLists
        )
        {
            intLists.resize(1, intCells_[Pstream::myProcNo()]);
        };

        virtual void getReferencedSurfaceCellList(
            List<DynamicLabelList>& surfLists
        )
        {
            surfLists.resize(1, surfCells_[Pstream::myProcNo()]);
        };

        virtual void getReferencedHaloCellList(
            List<DynamicLabelList>& haloLists
        )
        {
            haloLists.resize(1, haloCells_[Pstream::myProcNo()]);
        };

        virtual void getReferencedCoMList(
            DynamicVectorList& referenceCoM
        )
        {
            referenceCoM.resize(1, getCoM());
        };

        virtual vector addModelReturnRandomPosition(
            const bool allActiveCellsInMesh,
            const boundBox  cellZoneBounds,
            Random&         randGen
        )
        {
            InfoH << iB_Info << "2-- geomModelMessage-- "
                << "notImplemented"<<endl;
            return vector::zero;
        }

        virtual void bodyMovePoints(vector translVec)
        {
            InfoH << iB_Info << "3-- geomModelMessage-- "
                << "notImplemented"<<endl;
        }

        virtual void bodyScalePoints(scalar scaleFac)
        {
            InfoH << iB_Info << "4-- geomModelMessage-- "
                << "notImplemented"<<endl;
        }

        virtual void bodyRotatePoints
        (
            scalar rotAngle,
            vector axisOfRot
        )
        {
            InfoH << iB_Info << "5-- geomModelMessage-- "
                << "notImplemented"<<endl;
        }

        virtual void resetBody()
        {
            InfoH << iB_Info << "6-- geomModelMessage-- "
                << "notImplemented"<<endl;
        }

        virtual vector getCoM()
        {
            InfoH << iB_Info << "7-- geomModelMessage-- "
                << "notImplemented"<<endl;
            return vector::zero;
        }

        virtual boundBox getBounds()
        {
            InfoH << iB_Info << "8-- geomModelMessage-- "
                << "notImplemented"<<endl;
            boundBox bound;
            return bound;
        }

        virtual void synchronPos(label owner = -1)
        {
            InfoH << iB_Info << "9-- geomModelMessage-- "
                << "notImplemented"<<endl;
        }

        virtual void writeBodyInfo(label bodyId, IOdictionary& dict, word path)
        {
            InfoH << iB_Info << "10-- geomModelMessage-- "
                << "notImplemented"<<endl;
        }

        virtual boolList pointInside(pointField pointI)
        {
            InfoH << iB_Info << "11-- geomModelMessage-- "
                << "notImplemented"<<endl;
            return {false};
        }

        virtual bool pointInside(point pointI)
        {
            InfoH << iB_Info << "12-- geomModelMessage-- "
                << "notImplemented"<<endl;
            return false;
        }

        virtual volumeType getVolumeType(subVolume& sv, bool cIb)
        {
            InfoH << iB_Info << "13-- geomModelMessage-- "
                << "notImplemented"<<endl;
            return volumeType::UNKNOWN;
        }

        virtual scalar getDC()
        {
            return dC_;
        }

        virtual void getClosestPointAndNormal
        (
            const point& startPoint,
            const vector& span,
            point& closestPoint,
            vector& normal
        )
        {
            InfoH << iB_Info << "14-- geomModelMessage-- "
                << "notImplemented"<<endl;
        }

        virtual bool limitFinalSubVolume
        (
            const subVolume& sv,
            bool cIb,
            boundBox& limBBox
        )
        {
            InfoH << iB_Info << "15-- geomModelMessage-- "
                << "notImplemented"<<endl;
            return false;
        }

        virtual void setCoM()
        {
            InfoH << iB_Info << "16-- geomModelMessage-- "
                << "notImplemented"<<endl;
        }

        constexpr contactType getcType()
        {return contactType_;}

        virtual void setOwner()
        {
            Pstream::gatherList(ibPartialVolume_, 0);
            Pstream::broadcast(ibPartialVolume_, 0);                    //OF.com scatter -> broadcast
            for (label i = 0; i < ibPartialVolume_.size(); i++)
            {
                if (ibPartialVolume_[i] == max(ibPartialVolume_))
                {
                    //set owner of the IB which will move this IB
                    owner_ = i;
                    break;
                }
            }
        }

        virtual label getOwner()
        {
            return owner_;
        }

        void setIntSpan(scalar intSpan){intSpan_ = intSpan;}

        void setSdBasedLambda(bool sdBasedLambda){sdBasedLambda_ = sdBasedLambda;}

        void setOverlayLambda(word overlayLambda){overlayLambda_ = overlayLambda;}

        void resetHashTable(){cachedNeighbours_.clear();}

        virtual void setRhoS(dimensionedScalar rhoS){rhoS_ = rhoS;}
        dimensionedScalar getRhoS(){return rhoS_;}

        virtual void calculateGeometricalProperties(
            volScalarField& body
        );

        virtual void calculateGeometricalPropertiesParallel(
            volScalarField& body
        );



        //evaluate center of mass and inertia tensor
        void addToMAndI
        (
            volScalarField& body,
            DynamicLabelList& labelCellLst
        );

        void computeBodyCharPars();

        virtual scalar& getM0()
        {
            return M0_;
        }

        void setM0(scalar M0)
        {
            M0_ = M0;
        }

        scalar& getM()
        {
            return M_;
        }

        void setM(scalar M)
        {
            M_ = M;
        }


        label& getNCells()
        {
            return nCells_;
        }

        void setNCells(label nCells)
        {
            nCells_ = nCells;
        }

        const dimensionedScalar getRhoS()    const
        {
            return rhoS_;
        }

        virtual vector getLVec(const point& toPoint)
        {
            return toPoint - getCoM();
        }

        const symmTensor& getI()    const
        {
            return I_;
        }

        void setI(const symmTensor& I)
        {
            I_ = I;
        }

        virtual void resetBody(volScalarField& body);

        virtual bool isCluster() const
        {
            return false;
        }

        virtual bool checkCluster()
        {
            InfoH << iB_Info << "17-- geomModelMessage-- "
                << "not cluster model"<<endl;
                return false;
        }

        bool isBBoxInMesh();

        DynamicLabelList getPotentSurfCells
        (
            volScalarField& body,
            HashTable<bool, label, Hash<label>>& cellInside,
            List<labelList>& cellPoints
        );

        void correctSurfCells
        (
            volScalarField& body,
            DynamicLabelList& potentSurfCells,
            HashTable<bool, label, Hash<label>>& cellInside,
            List<labelList>& cellPoints
        );

        void findHaloCells
        (
            volScalarField& body
        );

        virtual List<std::shared_ptr<boundBox>> getBBoxes();

        virtual pointField getBodyPoints()
        {
            InfoH << iB_Info << "18-- geomModelMessage-- "
                << "not cluster model"<<endl;
            pointField emptyField;
            return emptyField;
        }

        virtual  void setBodyPosition(pointField pos)
        {
            InfoH << iB_Info << "22-- geomModelMessage-- "
                << "notImplemented"<<endl;
        }

        void setRefineBuffers(label refineBuffers)
        {
            refineBuffers_ = refineBuffers;
        }

        label getRefineBuffers()
        {
            return refineBuffers_;
        }

    };
}
#endif
