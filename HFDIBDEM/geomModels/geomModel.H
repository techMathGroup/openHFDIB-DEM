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
    addModels.C

Contributors
    Martin Isoz (2019-*), Martin Šourek (2019-*), 
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef geomModel_H
#define geomModel_H

#include "dictionary.H"
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "triSurfaceMesh.H"

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
        
        protected:
            
        const  dynamicFvMesh&   mesh_; 
        //list to find the owner of the IB
        DynamicLabelList ibPartialVolume_;
        label owner_;
        
        label cellToStartInCreateIB_;
        scalar  thrSurf_;
        Vector<label> geometricD_;//geometric directions
        
        scalar  intSpan_;        
        bool    sdBasedLambda_;
        boundBox curMeshBounds_;
        HashTable<pointField,label,Hash<label>> lastIbPoints_;
        
        public:
        
        // Constructors ------------------------------------------------
        geomModel
        (
            const  dynamicFvMesh&   mesh,
            scalar  thrSurf = 1,
            Vector<label> geometricD = {1,1,1}
        );
        
        // Destructors -------------------------------------------------
        virtual ~geomModel();
        
        virtual void createImmersedBody
        (
            volScalarField& body,
            Field<label>& octreeField,
            List<DynamicLabelList>& surfCells,
            List<DynamicLabelList>& intCells,
            List<pointField>& cellPoints
        )
        {Info << "-- shapeModelMessage-- " << "notIplemented"<<endl;}    
        
        virtual bool canAddBody
        (
            const volScalarField& body
        )
        {Info << "-- shapeModelMessage-- " << "notIplemented"<<endl; return false;} 
        
        virtual geomModel* getGeomModel()
        {
            Info << "-- shapeModelMessage-- " << "notIplemented"<<endl;
            geomModel* bodyGeomModel = new geomModel(mesh_);
            return bodyGeomModel;            
        } 
        
        virtual vector addModelReturnRandomPosition(
            const bool allActiveCellsInMesh,
            const boundBox  cellZoneBounds,       
            Random&         randGen
        )
        {Info << "-- shapeModelMessage-- " << "notIplemented"<<endl; return vector::zero;}  
        
        virtual void bodyMovePoints(vector translVec)        
        {Info << "-- shapeModelMessage-- " << "notIplemented"<<endl;} 
        
        virtual void bodyScalePoints(scalar scaleFac)        
        {Info << "-- shapeModelMessage-- " << "notIplemented"<<endl;}
        
        virtual void bodyRotatePoints
        (
            scalar rotAngle,
            vector axisOfRot
        )      
        {Info << "-- shapeModelMessage-- " << "notIplemented"<<endl;}
            
        virtual void resetBody()
        {Info << "-- shapeModelMessage-- " << "notIplemented"<<endl;}
        
        virtual vector getCoM()
        {Info << "-- shapeModelMessage-- " << "notIplemented"<<endl; return vector::zero;}
        
        virtual boundBox getBounds()
        {Info << "-- shapeModelMessage-- " << "notIplemented"<<endl; 
            boundBox bound;
            return bound;}
        
         virtual void recordHistory()
        {Info << "-- shapeModelMessage-- " << "notIplemented"<<endl;}
        
         virtual void returnHistory()
        {Info << "-- shapeModelMessage-- " << "notIplemented"<<endl;}
        
        virtual void synchronPos()
        {Info << "-- shapeModelMessage-- " << "notIplemented"<<endl;}
        
        virtual void writeBodyInfo(label bodyId, IOdictionary& dict, word path)
        {Info << "-- shapeModelMessage-- " << "notIplemented"<<endl;}
        
        virtual boolList pointInside(pointField pointI)
        {Info << "-- shapeModelMessage-- " << "notIplemented"<<endl; return {false};}
        
        virtual vector getClosestNormal(point startPoint, vector span)
        {Info << "-- shapeModelMessage-- " << "notIplemented"<<endl; return vector::zero;}
                
        label getOwner()
        {    
            Pstream::gatherList(ibPartialVolume_, 0);
            Pstream::scatter(ibPartialVolume_, 0);
            for (label i = 0; i < ibPartialVolume_.size(); i++)
            {
                if (ibPartialVolume_[i] == max(ibPartialVolume_))
                {
                    //Set owner of the IB which will move this IB
                    owner_ = i;
                    break;
                }
            }
            return owner_;
        }
        
        void setIntSpan(scalar intSpan){intSpan_ = intSpan;}
        
        void setSdBasedLambda(bool sdBasedLambda){sdBasedLambda_ = sdBasedLambda;}
        
        void resetHashTable(){lastIbPoints_.clear();}
    };
}
#endif
