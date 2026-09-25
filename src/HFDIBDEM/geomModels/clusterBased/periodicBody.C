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
#include "periodicBody.H"

using namespace Foam;

//---------------------------------------------------------------------------//
void periodicBody::getReferencedLists
(
    List<DynamicLabelList>& intLists,
    List<DynamicLabelList>& surfLists,
    DynamicVectorList& referenceCoM
)
{
    intLists.resize(ibGeomModelList.size());
    surfLists.resize(ibGeomModelList.size());
    referenceCoM.resize(ibGeomModelList.size());

    for(size_t ibI = 0; ibI < ibGeomModelList.size(); ++ibI)
    {
        intLists[ibI] =
            ibGeomModelList[ibI]->getInternalCellList()[Pstream::myProcNo()];

        surfLists[ibI] =
            ibGeomModelList[ibI]->getSurfaceCellList()[Pstream::myProcNo()];

        referenceCoM[ibI] =
            ibGeomModelList[ibI]->getCoM();
    }
}
//---------------------------------------------------------------------------//
void periodicBody::bodyMovePoints
(
    vector translVec
)
{
    for(std::shared_ptr<geomModel>& gModel : ibGeomModelList)
    {
        gModel->bodyMovePoints(translVec);
    }
}
//---------------------------------------------------------------------------//
void periodicBody::bodyScalePoints
(
    scalar scaleFac
)
{
    for(std::shared_ptr<geomModel>& gModel : ibGeomModelList)
    {
        gModel->bodyScalePoints(scaleFac);
    }
}
//---------------------------------------------------------------------------//
void periodicBody::bodyRotatePoints
(
    scalar rotAngle,
    vector axisOfRot
)
{
    for(std::shared_ptr<geomModel>& gModel : ibGeomModelList)
    {
        gModel->bodyRotatePoints(rotAngle, axisOfRot);
    }
}
//---------------------------------------------------------------------------//
vector periodicBody::getCoM()
{
    return ibGeomModelList[0]->getCoM();
}
//---------------------------------------------------------------------------//
scalar periodicBody::getDC()
{
    scalar dc = ibGeomModelList[0]->getDC();
    for(size_t ibI = 0; ibI < ibGeomModelList.size(); ++ibI)
    {
        if(dc < ibGeomModelList[ibI]->getDC())
        {
            dc = ibGeomModelList[ibI]->getDC();
        }
    }
    return dc;
}
//---------------------------------------------------------------------------//
void periodicBody::setOwner()
{
    HashTable<label, label, Hash<label>> frequency;
    for(std::shared_ptr<geomModel>& gModel : ibGeomModelList)
    {
        gModel->setOwner();

        label ownerI = gModel->getOwner();
        if(!frequency.found(ownerI))
        {
            frequency.insert(ownerI, 0);
        }

        frequency[ownerI]++;
    }

    label maxFreq = 0;
    label maxProc = 0;
    for (auto it = frequency.begin(); it != frequency.end(); ++it)
    {
        if (*it > maxFreq)
        {
            maxFreq = *it;
            maxProc = it.key();
        }
    }

    owner_ = maxProc;
}
//---------------------------------------------------------------------------//
label periodicBody::getOwner()
{
    return owner_;
}
//---------------------------------------------------------------------------//
scalar& periodicBody::getM0()
{
    return ibGeomModelList[0]->getM0();
}
//---------------------------------------------------------------------------//
vector periodicBody::getLVec(const point& toPoint)
{
    List<point> closestPoints(ibGeomModelList.size());
    List<vector> closestNormals(ibGeomModelList.size());

    for(size_t ibI = 0; ibI < ibGeomModelList.size(); ++ibI)
    {
        ibGeomModelList[ibI]->getClosestPointAndNormal(
            toPoint,
            (ibGeomModelList[ibI]->getCoM() - toPoint),
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

    return toPoint - ibGeomModelList[cGModel]->getCoM();
}
//---------------------------------------------------------------------------//
bool periodicBody::shouldBeUnclustered()
{
    int remBodies = 0;
    Info << "Periodic body " << ibGeomModelList.size() << " size." << endl;
    for(std::shared_ptr<geomModel>& gModel : ibGeomModelList)
    {
        Info << "Periodic body " << gModel->getM() << " size." << endl;
        if(gModel->getM() > 0)
        {
            ++remBodies;
        }
    }

    return remBodies == 1;
}
//---------------------------------------------------------------------------//
std::shared_ptr<geomModel> periodicBody::getRemGeomModel()
{
    auto iter = find_if(ibGeomModelList.begin(), ibGeomModelList.end(),
        [] (const std::shared_ptr<geomModel>& gP) { return gP->getM() > 0; } );

    if (iter != ibGeomModelList.end())
    {
        return *iter;
    }
    else
    {
        return ibGeomModelList[0];
    }
}
//---------------------------------------------------------------------------//
