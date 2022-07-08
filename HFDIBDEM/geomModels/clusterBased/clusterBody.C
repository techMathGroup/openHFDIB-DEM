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

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "clusterBody.H"

using namespace Foam;

//---------------------------------------------------------------------------//
// create immersed body for convex body
void clusterBody::createImmersedBody
(
    volScalarField& body,
    Field<label>& octreeField,
    List<labelList>& cellPoints
)
{
    forAll(ibGeomModelList, ibI)
    {
        ibGeomModelList[ibI].createImmersedBody(
            body,
            octreeField,
            cellPoints
        );
    }
}
//---------------------------------------------------------------------------//
void clusterBody::getReferencedLists
(
    List<DynamicLabelList>& intLists,
    List<DynamicLabelList>& surfLists,
    DynamicVectorList& referenceCoM
)
{
    intLists.resize(ibGeomModelList.size());
    surfLists.resize(ibGeomModelList.size());
    referenceCoM.resize(ibGeomModelList.size());

    forAll(ibGeomModelList, ibI)
    {
        intLists[ibI] = 
            ibGeomModelList[ibI].getInternalCellList()[Pstream::myProcNo()];
        
        surfLists[ibI] = 
            ibGeomModelList[ibI].getSurfaceCellList()[Pstream::myProcNo()];

        referenceCoM[ibI] = 
            ibGeomModelList[ibI].getCoM();
    }
}
//---------------------------------------------------------------------------//
void clusterBody::updateSurfList()
{
    surfCells_.clear();
    surfCells_.setSize(Pstream::nProcs());

    forAll(ibGeomModelList, ibI)
    {
        surfCells_[Pstream::myProcNo()].append(
            ibGeomModelList[ibI].getSurfaceCellList()[Pstream::myProcNo()]
        );
    }
}
//---------------------------------------------------------------------------//
void clusterBody::updateIntList()
{
    intCells_.clear();
    intCells_.setSize(Pstream::nProcs());

    forAll(ibGeomModelList, ibI)
    {
        intCells_[Pstream::myProcNo()].append(
            ibGeomModelList[ibI].getInternalCellList()[Pstream::myProcNo()]
        );
    }
}
//---------------------------------------------------------------------------//
void clusterBody::calculateGeometricalProperties
(
    volScalarField& body
)
{
    M_      = scalar(0);
    I_      = symmTensor::zero;

    forAll(ibGeomModelList, ibI)
    {
        ibGeomModelList[ibI].calculateGeometricalProperties(body);
        M_ += ibGeomModelList[ibI].getM();
        I_ += ibGeomModelList[ibI].getI();
    }
}
//---------------------------------------------------------------------------//
void clusterBody::bodyMovePoints
(
    vector translVec
)
{
    forAll(ibGeomModelList, ibI)
    {
        if(ibGeomModelList[ibI].getOwner() == Pstream::myProcNo())
        {
            ibGeomModelList[ibI].bodyMovePoints(translVec);
        }
    }
}
//---------------------------------------------------------------------------//
void clusterBody::bodyScalePoints
(
    scalar scaleFac
)
{
    forAll(ibGeomModelList, ibI)
    {
        ibGeomModelList[ibI].bodyScalePoints(scaleFac);
    }
}
//---------------------------------------------------------------------------//
void clusterBody::bodyRotatePoints
(
    scalar rotAngle,
    vector axisOfRot
)
{
    forAll(ibGeomModelList, ibI)
    {
        if(ibGeomModelList[ibI].getOwner() == Pstream::myProcNo())
        {
            ibGeomModelList[ibI].bodyRotatePoints(rotAngle, axisOfRot);
        }
    }
}
//---------------------------------------------------------------------------//
vector clusterBody::getCoM()
{    
    return ibGeomModelList[0].getCoM();
}
//---------------------------------------------------------------------------//
boundBox clusterBody::getBounds()
{
    DynamicPointList allBounds;
    forAll(ibGeomModelList, ibI)
    {
        boundBox bBoxI = ibGeomModelList[ibI].getBounds();
        allBounds.append(bBoxI.min());
        allBounds.append(bBoxI.max());
    }
    
    return boundBox(allBounds);
}
//---------------------------------------------------------------------------//
void clusterBody::synchronPos()
{
    forAll(ibGeomModelList, ibI)
    {
        ibGeomModelList[ibI].synchronPos();
    }
}
//---------------------------------------------------------------------------//
boolList clusterBody::pointInside(pointField pointF)
{
    boolList inside(pointF.size());

    forAll(pointF,pointI)
    {
        bool pointInside = false;
        forAll(ibGeomModelList, ibI)
        {
            pointInside = ibGeomModelList[ibI].pointInside(pointF[pointI]);
            if(pointInside)
                break;
        }
        inside[pointI] = pointInside;
    }

    return inside;
}
//---------------------------------------------------------------------------//
bool clusterBody::pointInside(point pointI)
{
    bool pointInside = false;
    forAll(ibGeomModelList, ibI)
    {
        pointInside = ibGeomModelList[ibI].pointInside(pointI);
        if(pointInside)
            break;
    }
    return pointInside;
}
//---------------------------------------------------------------------------//
scalar clusterBody::getDC()
{
    scalar dc = ibGeomModelList[0].getDC();
    for(int i = 1; i < ibGeomModelList.size(); ++i)
    {
        if(dc < ibGeomModelList[i].getDC())
        {
            dc = ibGeomModelList[i].getDC();
        }
    }
    return dc;
}
//---------------------------------------------------------------------------//
void clusterBody::getClosestPointAndNormal
(
    const point& startPoint,
    const vector& span,
    point& closestPoint,
    vector& normal
)
{
    List<point> closestPoints(ibGeomModelList.size());
    List<vector> closestNormals(ibGeomModelList.size());

    forAll(ibGeomModelList, ibI)
    {
        ibGeomModelList[ibI].getClosestPointAndNormal(
            startPoint,
            span,
            closestPoints[ibI],
            closestNormals[ibI]
        );
    }

    closestPoint = closestPoints[0];
    normal = closestNormals[0];
    for (int i = 1; i < closestPoints.size(); ++i)
    {
        if(mag(startPoint - closestPoint) 
            > mag(startPoint - closestPoints[i]))
        {
            closestPoint = closestPoints[i];
            normal = closestNormals[i];
        }
    }
}
//---------------------------------------------------------------------------//
label clusterBody::getOwner()
{
    return Pstream::myProcNo();
}
//---------------------------------------------------------------------------//
scalar& clusterBody::getM0()
{
    return ibGeomModelList[0].getM0();
}
//---------------------------------------------------------------------------//
vector clusterBody::getLVec(const point& toPoint)
{
    List<point> closestPoints(ibGeomModelList.size());
    List<vector> closestNormals(ibGeomModelList.size());

    forAll(ibGeomModelList, ibI)
    {
        ibGeomModelList[ibI].getClosestPointAndNormal(
            toPoint,
            (ibGeomModelList[ibI].getCoM() - toPoint),
            closestPoints[ibI],
            closestNormals[ibI]
        );
    }

    vector closestPoint = closestPoints[0];
    label cGModel = 0;

    for(int i = 1; i < closestPoints.size(); ++i)
    {
        if(mag(toPoint - closestPoint) 
            > mag(toPoint - closestPoints[i]))
        {
            closestPoint = closestPoints[i];
            cGModel = i;
        }
    }

    return toPoint - ibGeomModelList[cGModel].getCoM();
}
//---------------------------------------------------------------------------//
void clusterBody::resetBody(volScalarField& body)
{
    forAll(ibGeomModelList, ibI)
    {
        ibGeomModelList[ibI].resetBody(body);
    }
}
//---------------------------------------------------------------------------//
bool clusterBody::shouldBeUnclustered()
{
    int remBodies = 0;
    forAll(ibGeomModelList, ibI)
    {
        if(ibGeomModelList[ibI].getM() > 0)
        {
            ++remBodies;
        }
    }

    if(remBodies == 1)
    {
        return true;
    }
    return false;
}
//---------------------------------------------------------------------------//
autoPtr<geomModel> clusterBody::getRemGeomModel()
{
    forAll(ibGeomModelList, ibI)
    {
        if(ibGeomModelList[ibI].getM() > 0)
        {
            return ibGeomModelList.set(ibI, nullptr);
        }
    }

    return ibGeomModelList.set(0, nullptr);
}
