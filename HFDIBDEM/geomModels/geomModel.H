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

#include "dictionary.H"
#include "fvCFD.H"
// #include "fvMesh.H"
#include "triSurfaceMesh.H"

#include "contactType.H"

#include "outputHFDIBDEM.H"
using namespace outputHFDIBDEM;

#include "clockTime.H"

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

        //list to find the owner of the IB
        DynamicLabelList ibPartialVolume_;
        label owner_;

        label cellToStartInCreateIB_;
        scalar  thrSurf_;

        scalar  intSpan_;
        bool    sdBasedLambda_;
        boundBox curMeshBounds_;
        autoPtr<HashTable<const labelList&, label, Hash<label>>> cachedNeighbours_;

        //body mass
        scalar M_;
        scalar M0_;                                                     //initial mass (at addition)

        //location of center of mass
        vector CoM_;

        //moment of intertia tensor
        symmTensor  I_;

        //characteristic diameter of the body
        scalar  dC_;

        //body density
        dimensionedScalar rhoS_;

        public:

        static scalar t1;
        static scalar t2;
        static scalar t3;
        static scalar t4;
        static scalar t5;
        static scalar t6;
        static scalar t7;

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
            List<DynamicLabelList>& surfCells,
            List<DynamicLabelList>& intCells,
            List<labelList>& cellPoints
        )
        {
            InfoH << iB_Info << "-- geomModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual bool canAddBody
        (
            const volScalarField& body
        )
        {
            InfoH << iB_Info << "-- geomModelMessage-- "
                << "notIplemented"<<endl;
            return false;
        }

        virtual geomModel* getGeomModel() = 0;

        virtual vector addModelReturnRandomPosition(
            const bool allActiveCellsInMesh,
            const boundBox  cellZoneBounds,
            Random&         randGen
        )
        {
            InfoH << iB_Info << "-- geomModelMessage-- "
                << "notIplemented"<<endl;
            return vector::zero;
        }

        virtual void bodyMovePoints(vector translVec)
        {
            InfoH << iB_Info << "-- geomModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual void bodyScalePoints(scalar scaleFac)
        {
            InfoH << iB_Info << "-- geomModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual void bodyRotatePoints
        (
            scalar rotAngle,
            vector axisOfRot
        )
        {
            InfoH << iB_Info << "-- geomModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual void resetBody()
        {
            InfoH << iB_Info << "-- geomModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual vector getCoM()
        {
            InfoH << iB_Info << "-- geomModelMessage-- "
                << "notIplemented"<<endl;
            return vector::zero;
        }

        virtual boundBox getBounds()
        {
            InfoH << iB_Info << "-- geomModelMessage-- "
                << "notIplemented"<<endl;
            boundBox bound;
            return bound;
        }

        virtual void synchronPos()
        {
            InfoH << iB_Info << "-- geomModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual void writeBodyInfo(label bodyId, IOdictionary& dict, word path)
        {
            InfoH << iB_Info << "-- geomModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual boolList pointInside(pointField pointI)
        {
            InfoH << iB_Info << "-- geomModelMessage-- "
                << "notIplemented"<<endl;
            return {false};
        }

        virtual bool pointInside(point pointI)
        {
            InfoH << iB_Info << "-- geomModelMessage-- "
                << "notIplemented"<<endl;
            return false;
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
            InfoH << iB_Info << "-- geomModelMessage-- "
                << "notIplemented"<<endl;
        }

        virtual void setCoM()
        {
            InfoH << iB_Info << "-- geomModelMessage-- "
                << "notIplemented"<<endl;
        }

        constexpr contactType getcType()
        {return contactType_;}

        label getOwner()
        {
            Pstream::gatherList(ibPartialVolume_, 0);
            Pstream::scatter(ibPartialVolume_, 0);
            for (label i = 0; i < ibPartialVolume_.size(); i++)
            {
                if (ibPartialVolume_[i] == max(ibPartialVolume_))
                {
                    //set owner of the IB which will move this IB
                    owner_ = i;
                    break;
                }
            }
            return owner_;
        }

        void setIntSpan(scalar intSpan){intSpan_ = intSpan;}

        void setSdBasedLambda(bool sdBasedLambda){sdBasedLambda_ = sdBasedLambda;}

        void resetHashTable(){cachedNeighbours_.clear();}

        virtual void setRhoS(dimensionedScalar rhoS){rhoS_ = rhoS;}
        dimensionedScalar getRhoS(){return rhoS_;}

        void calculateGeometricalProperties(
            volScalarField& body,
            List<DynamicLabelList>& surfCells,
            List<DynamicLabelList>& intCells
        );

        //evaluate center of mass and inertia tensor
        void addToMAndI
        (
            volScalarField& body,
            DynamicLabelList& labelCellLst
        );

        void computeBodyCharPars(List<DynamicLabelList>& surfCells);

        virtual scalar& getM0()
        {
            return M0_;
        }

        scalar& getM()
        {
            return M_;
        }

        const dimensionedScalar getRhoS()    const
        {
            return rhoS_;
        }

        const vector& getCoM()  const
        {
            return CoM_;
        }

        const symmTensor& getI()    const
        {
            return I_;
        }

        bool isBBoxInMesh();

        DynamicLabelList getPotentSurfCells
        (
            volScalarField& body,
            List<DynamicLabelList>& intCells,
            HashTable<bool, label, Hash<label>>& cellInside,
            List<labelList>& cellPoints
        );

        void correctSurfCells
        (
            volScalarField& body,
            List<DynamicLabelList>& intCells,
            List<DynamicLabelList>& surfCells,
            DynamicLabelList& potentSurfCells,
            HashTable<bool, label, Hash<label>>& cellInside,
            List<labelList>& cellPoints
        );
    };
}
#endif
