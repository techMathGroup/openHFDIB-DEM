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
    for(std::shared_ptr<geomModel>& gModel : ibGeomModelList)
    {
        gModel->createImmersedBody(
            body,
            octreeField,
            cellPoints
        );

        Info << "Periodic body created" << " mass: " << gModel->getM() << endl;
        Info << "Periodic body created" << " bbox: " << gModel->getBounds().min() << " " << gModel->getBounds().max() << endl;
        Info << "Periodic body created" << " intList: " << gModel->getInternalCellList()[Pstream::myProcNo()].size() << endl;
    }
}
//---------------------------------------------------------------------------//
void clusterBody::updateSurfList()
{
    surfCells_.clear();
    surfCells_.setSize(Pstream::nProcs());

    for(std::shared_ptr<geomModel>& gModel : ibGeomModelList)
    {
        surfCells_[Pstream::myProcNo()].append(
            gModel->getSurfaceCellList()[Pstream::myProcNo()]
        );
    }
}
//---------------------------------------------------------------------------//
void clusterBody::updateIntList()
{
    intCells_.clear();
    intCells_.setSize(Pstream::nProcs());

    for(std::shared_ptr<geomModel>& gModel : ibGeomModelList)
    {
        intCells_[Pstream::myProcNo()].append(
            gModel->getInternalCellList()[Pstream::myProcNo()]
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

    for(std::shared_ptr<geomModel>& gModel : ibGeomModelList)
    {
        gModel->calculateGeometricalProperties(body);
        M_ += gModel->getM();
        I_ += gModel->getI();
    }
}
//---------------------------------------------------------------------------//
void clusterBody::calculateGeometricalPropertiesParallel
(
    volScalarField& body
)
{
    M_      = scalar(0);
    I_      = symmTensor::zero;

    for(std::shared_ptr<geomModel>& gModel : ibGeomModelList)
    {
        gModel->calculateGeometricalPropertiesParallel(body);

        Info << "Periodic body calculated" << " mass: " << gModel->getM() << endl;
    }
}
//---------------------------------------------------------------------------//
void clusterBody::setMassAndInertia()
{
    M_      = scalar(0);
    I_      = symmTensor::zero;

    for(std::shared_ptr<geomModel>& gModel : ibGeomModelList)
    {
        M_ += gModel->getM();
        I_ += gModel->getI();
    }
}
//---------------------------------------------------------------------------//
vector clusterBody::getCoM()
{
    return ibGeomModelList[0]->getCoM();
}
//---------------------------------------------------------------------------//
boundBox clusterBody::getBounds()
{
    DynamicPointList allBounds;
    for(std::shared_ptr<geomModel>& gModel : ibGeomModelList)
    {
        boundBox bBoxI = gModel->getBounds();
        allBounds.append(bBoxI.min());
        allBounds.append(bBoxI.max());
    }

    return boundBox(allBounds);
}
//---------------------------------------------------------------------------//
void clusterBody::synchronPos(label owner)
{
    for(std::shared_ptr<geomModel>& gModel : ibGeomModelList)
    {
        gModel->synchronPos(owner_);
    }
}
//---------------------------------------------------------------------------//
boolList clusterBody::pointInside(pointField pointF)
{
    boolList inside(pointF.size());

    forAll(pointF,pointI)
    {
        bool pointInside = false;
        for(std::shared_ptr<geomModel>& gModel : ibGeomModelList)
        {
            pointInside = gModel->pointInside(pointF[pointI]);
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
    for(std::shared_ptr<geomModel>& gModel : ibGeomModelList)
    {
        if(gModel->pointInside(pointI))
            return true;
    }
    return false;
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

    for(size_t ibI = 0; ibI < ibGeomModelList.size(); ++ibI)
    {
        ibGeomModelList[ibI]->getClosestPointAndNormal(
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
void clusterBody::resetBody(volScalarField& body)
{
    for(std::shared_ptr<geomModel>& gModel : ibGeomModelList)
    {
        gModel->resetBody(body);
    }
}
//---------------------------------------------------------------------------//
List<std::shared_ptr<boundBox>> clusterBody::getBBoxes()
{
    List<std::shared_ptr<boundBox>> retList;
    for(std::shared_ptr<geomModel>& gModel : ibGeomModelList)
    {
        List<std::shared_ptr<boundBox>> bBoxI = gModel->getBBoxes();
        retList.append(bBoxI);
    }
    return retList;
}
//---------------------------------------------------------------------------//
