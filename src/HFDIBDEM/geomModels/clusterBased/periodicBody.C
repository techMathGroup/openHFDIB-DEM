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
void periodicBody::bodyMovePoints
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
void periodicBody::bodyScalePoints
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
void periodicBody::bodyRotatePoints
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
vector periodicBody::getCoM()
{    
    return ibGeomModelList[0].getCoM();
}
//---------------------------------------------------------------------------//
scalar periodicBody::getDC()
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
void periodicBody::setOwner()
{
    forAll(ibGeomModelList, ibI)
    {
        ibGeomModelList[ibI].setOwner();
    }

    HashTable<label, label, Hash<label>> frequency;
    forAll(ibGeomModelList, ibI)
    {
        label ownerI = ibGeomModelList[ibI].getOwner();
        if(!frequency.found(ownerI))
        {
            frequency.insert(ownerI, 1);
        }
        else
        {
            frequency[ownerI]++;
        }
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
    return ibGeomModelList[0].getM0();
}
//---------------------------------------------------------------------------//
vector periodicBody::getLVec(const point& toPoint)
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
bool periodicBody::shouldBeUnclustered()
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
autoPtr<geomModel> periodicBody::getRemGeomModel()
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
//---------------------------------------------------------------------------//
