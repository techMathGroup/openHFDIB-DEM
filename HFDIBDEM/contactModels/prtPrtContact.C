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
    Martin Isoz (2019-*), Martin Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "prtPrtContact.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace contactModel
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//---------------------------------------------------------------------------//
bool detectPrtPrtContact(
    const dynamicFvMesh&   mesh,
    contactInfo& cInfo,
    contactInfo& tInfo
)
{
    const contactType cT(cInfo.getGeomModel().getcType());
    const contactType tT(tInfo.getGeomModel().getcType());

    if(cT == sphere && tT == sphere)
    {
        return detectPrtPrtContact<sphere,sphere>
        (
            mesh,
            cInfo,
            tInfo
        );
    }
    else
    {
        return detectPrtPrtContact<arbShape,arbShape>(
            mesh,
            cInfo,
            tInfo
        );
    }
}
//---------------------------------------------------------------------------//
template <contactType cT, contactType tT>
bool
detectPrtPrtContact(
    const dynamicFvMesh&   mesh,
    contactInfo& cInfo,
    contactInfo& tInfo
)
{
    List<DynamicLabelList>    commonCells;
    commonCells.setSize(Pstream::nProcs());

    List<DynamicLabelList> cSurfCells(cInfo.getSurfCells());
    List<DynamicLabelList> cIntCells(cInfo.getIntCells());
    List<DynamicLabelList> tSurfCells(tInfo.getSurfCells());

    // iterate over surfCells to find commen cells
    forAll (cSurfCells[Pstream::myProcNo()],cSCellI)
    {
        forAll (tSurfCells[Pstream::myProcNo()],tSCellI)
        {
            if (mag(cSurfCells[Pstream::myProcNo()][cSCellI]-tSurfCells[Pstream::myProcNo()][tSCellI]) < SMALL)
            {
                commonCells[Pstream::myProcNo()].append(cSurfCells[Pstream::myProcNo()][cSCellI]);
            }
        }
    }

    scalar intersectedVolume(0);

    if (commonCells[Pstream::myProcNo()].size() > SMALL)
    {
        const pointField& pp = mesh.points();

        boolList tcenterInsideList = cInfo.getGeomModel().pointInside(mesh.C());
        boolList ccenterInsideList = tInfo.getGeomModel().pointInside(mesh.C());

        // iterate over all surfCells and intCells to evaluate intersected volume
        forAll (cSurfCells[Pstream::myProcNo()],cSCellI)
        {
            const labelList& vertexLabels = mesh.cellPoints()[cSurfCells[Pstream::myProcNo()][cSCellI]];
            const pointField vertexPoints(pp,vertexLabels);
            boolList cvertexesInside = cInfo.getGeomModel().pointInside( vertexPoints );
            boolList tvertexesInside = tInfo.getGeomModel().pointInside( vertexPoints );
            bool ccenterInside(ccenterInsideList[cSurfCells[Pstream::myProcNo()][cSCellI]] );
            bool tcenterInside(tcenterInsideList[cSurfCells[Pstream::myProcNo()][cSCellI]] );
            scalar rVInSize(0.5/(tvertexesInside.size()+1));
            // Note: weight of a single vertex in the cell

            scalar partialVolume(0);
            forAll (tvertexesInside, verIn)
            {
                if (tvertexesInside[verIn]==true && cvertexesInside[verIn]==true)
                {
                    partialVolume += rVInSize;                          //fraction of cell covered
                }
            }

            if (tcenterInside==true && ccenterInside==true)
            {
                partialVolume += 0.5;                                   //fraction of cell covered
            }

            intersectedVolume += mesh.V()[cSurfCells[Pstream::myProcNo()][cSCellI]] * partialVolume;

            if (intersectedVolume > 0) break;
        }

        forAll (cIntCells[Pstream::myProcNo()],cSCellI)
        {
            const labelList& vertexLabels = mesh.cellPoints()[cIntCells[Pstream::myProcNo()][cSCellI]];
            const pointField vertexPoints(pp,vertexLabels);
            boolList cvertexesInside = cInfo.getGeomModel().pointInside( vertexPoints );
            boolList tvertexesInside = tInfo.getGeomModel().pointInside( vertexPoints );
            bool ccenterInside(ccenterInsideList[cIntCells[Pstream::myProcNo()][cSCellI]] );
            bool tcenterInside(tcenterInsideList[cIntCells[Pstream::myProcNo()][cSCellI]] );
            scalar rVInSize(0.5/(tvertexesInside.size()+1));
            // Note: weight of a single vertex in the cell

            scalar partialVolume(0);
            forAll (tvertexesInside, verIn)
            {
                if (tvertexesInside[verIn]==true && cvertexesInside[verIn]==true)
                {
                    partialVolume += rVInSize;                          //fraction of cell covered
                }
            }

            if (tcenterInside==true && ccenterInside==true)
            {
                partialVolume += 0.5;                                   //fraction of cell covered
            }

            intersectedVolume += mesh.V()[cIntCells[Pstream::myProcNo()][cSCellI]] * partialVolume;

            if (intersectedVolume > 0) break;
        }
    }

    reduce(intersectedVolume, sumOp<scalar>());

    if (intersectedVolume > 0)
    {
        return true;
    }

    return false;
}
//---------------------------------------------------------------------------//
template <>
bool detectPrtPrtContact <sphere,sphere>
(
    const dynamicFvMesh&   mesh,
    contactInfo& cInfo,
    contactInfo& tInfo
)
{
    vector cCenter(cInfo.getGeomModel().getCoM());
    vector tCenter(tInfo.getGeomModel().getCoM());
    scalar cRadius(cInfo.getGeomModel().getDC() / 2);
    scalar tRadius(tInfo.getGeomModel().getDC() / 2);

    if(mag(cCenter-tCenter) < (cRadius + tRadius))
    {
        return true;
    }
    return false;
}
//---------------------------------------------------------------------------//
template <contactType cT, contactType tT>
void getPrtPrtContactVars(
    const dynamicFvMesh&   mesh,
    contactInfo& cInfo,
    contactInfo& tInfo,
    vector geometricD,
    prtPrtContactVars& prtPrtCntVars
)
{
    List<DynamicLabelList>    commonCells;
    commonCells.setSize(Pstream::nProcs());

    List<DynamicLabelList> cSurfCells(cInfo.getSurfCells());
    List<DynamicLabelList> cIntCells(cInfo.getIntCells());
    List<DynamicLabelList> tSurfCells(tInfo.getSurfCells());

    // iterate over surfCells to find commen cells
    forAll (cSurfCells[Pstream::myProcNo()],cSCellI)
    {
        forAll (tSurfCells[Pstream::myProcNo()],tSCellI)
        {
            if (mag(cSurfCells[Pstream::myProcNo()][cSCellI]-tSurfCells[Pstream::myProcNo()][tSCellI]) < SMALL)
            {
                commonCells[Pstream::myProcNo()].append(cSurfCells[Pstream::myProcNo()][cSCellI]);
            }
        }
    }

    label numOfComCells(0);
    vector contactCenter(vector::zero);
    // if there are any common cells check if the surfaces are intersected
    if (commonCells[Pstream::myProcNo()].size() > SMALL)
    {
        // evaluate center of the contact area
        forAll (commonCells[Pstream::myProcNo()], cCell)
        {
            contactCenter += mesh.C()[commonCells[Pstream::myProcNo()][cCell]];
        }
        numOfComCells = commonCells[Pstream::myProcNo()].size();
    }

    sumReduce(contactCenter, numOfComCells);

    if (numOfComCells <= 0)
    {
        prtPrtCntVars.contactCenter_ = vector::zero;
        prtPrtCntVars.contactVolume_ = 0;
        prtPrtCntVars.contactNormal_ = vector::zero;
        prtPrtCntVars.contactArea_ = 0;
        return;
    }

    contactCenter /= numOfComCells;

    vector cCoM(cInfo.getGeomModel().getCoM());                         //center of mass (current)
    vector tCoM(tInfo.getGeomModel().getCoM());                         //center of mass (neighbor)
    scalar tDC(tInfo.getGeomModel().getDC());                           //characteristic diameter

    vector normalVector(tInfo.getGeomModel().getClosestNormal(contactCenter,vector::one * tDC));
    scalar contactArea(0);
    bool case3D(true);
    // check if the case is 3D
    forAll (geometricD, direction)
    {
        if (geometricD[direction] == -1)
        {
            case3D = false;
            break;
        }
    }
    // use edge cells to find contact area (better precision then surfcells)
    DynamicVectorList edgePoints;
    DynamicLabelList edgeCells;

    scalar intersectedVolume(0);

    if (commonCells[Pstream::myProcNo()].size() > SMALL)
    {
        const pointField& pp = mesh.points();

        boolList tcenterInsideList = tInfo.getGeomModel().pointInside( mesh.C());
        boolList ccenterInsideList = cInfo.getGeomModel().pointInside( mesh.C());

        forAll (cSurfCells[Pstream::myProcNo()],cSCellI)
        {
            bool partiallyInT(false);

            const labelList& vertexLabels = mesh.cellPoints()[cSurfCells[Pstream::myProcNo()][cSCellI]];
            const pointField vertexPoints(pp,vertexLabels);
            boolList tvertexesInside = tInfo.getGeomModel().pointInside( vertexPoints );
            boolList cvertexesInside = cInfo.getGeomModel().pointInside( vertexPoints );
            bool tcenterInside(tcenterInsideList[cSurfCells[Pstream::myProcNo()][cSCellI]] );
            bool ccenterInside(ccenterInsideList[cSurfCells[Pstream::myProcNo()][cSCellI]] );
            scalar rVInSize(0.5/tvertexesInside.size());
            // Note: weight of a single vertex in the cell

            DynamicVectorList edgePointsI;

            scalar partialVolume(0);
            forAll (tvertexesInside, verIn)
            {
                if (tvertexesInside[verIn]==true && cvertexesInside[verIn]==true)
                {
                    partialVolume += rVInSize;                          //fraction of cell covered
                    partiallyInT = true;
                    edgePointsI.append(vertexPoints[verIn]);
                }
            }

            if (tcenterInside==true && ccenterInside==true)
            {
                partialVolume += 0.5;                                   //fraction of cell covered
                partiallyInT = true;
                for(label i = 0; i < tvertexesInside.size(); i++)
                    edgePointsI.append(mesh.C()[cSurfCells[Pstream::myProcNo()][cSCellI]] );
            }
            // cells is edge cell when the cell is surfcell in both proccessors
            if (partialVolume + SMALL < 1 && partiallyInT)
            {
                edgeCells.append(cSurfCells[Pstream::myProcNo()][cSCellI]);
                vector edgePoint(vector::zero);
                forAll(edgePointsI,pointI)
                {
                    edgePoint += edgePointsI[pointI];
                }
                edgePoint /= edgePointsI.size();
                edgePoints.append(edgePoint);
            }


            intersectedVolume += mesh.V()[cSurfCells[Pstream::myProcNo()][cSCellI]] * partialVolume;
        }
        // calculate remaining intersected volume
        forAll (cIntCells[Pstream::myProcNo()],cSCellI)
        {
            const labelList& vertexLabels = mesh.cellPoints()[cIntCells[Pstream::myProcNo()][cSCellI]];
            const pointField vertexPoints(pp,vertexLabels);
            boolList tvertexesInside = tInfo.getGeomModel().pointInside( vertexPoints );
            boolList cvertexesInside = cInfo.getGeomModel().pointInside( vertexPoints );
            bool tcenterInside(tcenterInsideList[cIntCells[Pstream::myProcNo()][cSCellI]] );
            bool ccenterInside(ccenterInsideList[cIntCells[Pstream::myProcNo()][cSCellI]] );
            scalar rVInSize(0.5/tvertexesInside.size());
            // Note: weight of a single vertex in the cell

            scalar partialVolume(0);
            forAll (tvertexesInside, verIn)
            {
                if (tvertexesInside[verIn]==true && cvertexesInside[verIn]==true)
                {
                    partialVolume += rVInSize;                          //fraction of cell covered
                }
            }

            if (tcenterInside==true && ccenterInside==true)
            {
                partialVolume += 0.5;                                   //fraction of cell covered
            }

            intersectedVolume += mesh.V()[cIntCells[Pstream::myProcNo()][cSCellI]] * partialVolume;
        }
    }

    reduce(intersectedVolume, sumOp<scalar>());

    if (intersectedVolume > 0)
    {
        if (case3D)
        {
            Tuple2<scalar,vector> returnTuple = get3DcontactInfo(mesh, edgeCells, edgePoints, normalVector, contactCenter, cInfo.getGeomModel().getOwner());
            contactArea = returnTuple.first();
            normalVector = returnTuple.second();
        }
        else
        {
            Tuple2<scalar,vector> returnTuple = get2DcontactInfo(mesh, commonCells[Pstream::myProcNo()], normalVector, contactCenter, geometricD);
            contactArea = returnTuple.first();
            normalVector = returnTuple.second();
        }
    }

    if (intersectedVolume > 0 && contactArea > 0)
    {
        prtPrtCntVars.contactCenter_ = contactCenter;
        prtPrtCntVars.contactVolume_ = intersectedVolume;
        prtPrtCntVars.contactNormal_ = normalVector;
        prtPrtCntVars.contactArea_ = contactArea;
    }
    else
    {
        prtPrtCntVars.contactCenter_ = vector::zero;
        prtPrtCntVars.contactVolume_ = 0;
        prtPrtCntVars.contactNormal_ = vector::zero;
        prtPrtCntVars.contactArea_ = 0;
    }
}
//---------------------------------------------------------------------------//
Tuple2<scalar,vector> get3DcontactInfo(
    const dynamicFvMesh&   mesh,
    DynamicLabelList commonCells,
    DynamicVectorList edgePoints,
    vector normalVector,
    vector contactCenter,
    label owner
)
{
    //Collect edge positions from all processors
    List<DynamicPointList> commCellsPositionsProc;
    commCellsPositionsProc.setSize(Pstream::nProcs());
    DynamicPointList commCellsPositions;
    forAll (commonCells, cCell)
    {
        commCellsPositionsProc[Pstream::myProcNo()].append(mesh.C()[commonCells[cCell]]);
    }

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UOPstream send(proci, pBufs);
            send << commCellsPositionsProc[Pstream::myProcNo()];
        }
    }

    pBufs.finishedSends();
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recv(proci, pBufs);
            DynamicPointList commCellsPositionsi (recv);
            commCellsPositionsProc[proci].append(commCellsPositionsi);
        }
    }

    for (label i = 0; i < commCellsPositionsProc[owner].size(); i++)
    {
        commCellsPositions.append(commCellsPositionsProc[owner][i]);
    }

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != owner)
        {
            for (label i = 0; i < commCellsPositionsProc[proci].size(); i++)
            {
                commCellsPositions.append(commCellsPositionsProc[proci][i]);
            }
        }
    }

    pBufs.clear();

    // collect edge Points from all processors
    List<DynamicPointList> commPointsPositionsProc;
    commPointsPositionsProc.setSize(Pstream::nProcs());
    DynamicPointList commPointsPositions;
    forAll (edgePoints, cCell)
    {
        commPointsPositionsProc[Pstream::myProcNo()].append(edgePoints[cCell]);
    }

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UOPstream send(proci, pBufs);
            send << commPointsPositionsProc[Pstream::myProcNo()];
        }
    }

    pBufs.finishedSends();
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recv(proci, pBufs);
            DynamicPointList commPointsPositionsi (recv);
            commPointsPositionsProc[proci].append(commPointsPositionsi);
        }
    }

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        for (label i = 0; i < commCellsPositionsProc[proci].size(); i++)
        {
            commPointsPositions.append(commPointsPositionsProc[proci][i]);
        }
    }

    scalar area(0.0);
    vector normalVec(vector::zero);

    if (commPointsPositions.size() >= 3)
    {
        bool normOk(false);
        vector center(vector::zero);

        forAll (commPointsPositions,cell)
        {
            center += commPointsPositions[cell];
        }
        center /= commPointsPositions.size();

        scalar xx(0);
        scalar xy(0);
        scalar xz(0);
        scalar yy(0);
        scalar yz(0);
        scalar zz(0);

        forAll (commPointsPositions,cell)
        {
            vector subPoint(commPointsPositions[cell] - center);
            if(subPoint != vector::zero)
                subPoint = subPoint/mag(subPoint);
            xx += subPoint[0] * subPoint[0];
            xy += subPoint[0] * subPoint[1];
            xz += subPoint[0] * subPoint[2];
            yy += subPoint[1] * subPoint[1];
            yz += subPoint[1] * subPoint[2];
            zz += subPoint[2] * subPoint[2];
        }

        xx /= commPointsPositions.size();
        xy /= commPointsPositions.size();
        xz /= commPointsPositions.size();
        yy /= commPointsPositions.size();
        yz /= commPointsPositions.size();
        zz /= commPointsPositions.size();

        vector weightedDir(vector::zero);


        scalar detX(yy*zz-yz*yz);
        vector axisDirX(detX,xz*yz-xy*zz,xy*yz-xz*yy);
        scalar weightX(detX*detX);
        if((weightedDir & axisDirX) < 0.0)
            weightX = -weightX;
        weightedDir += axisDirX * weightX;

        scalar detY(xx*zz-xz*xz);
        vector axisDirY(xz*yz-xy*zz,detY,xy*xz-yz*xx);
        scalar weightY(detY*detY);
        if((weightedDir & axisDirY) < 0.0)
            weightY = -weightY;
        weightedDir += axisDirY * weightY;

        scalar detZ(xx*yy-xy*xy);
        vector axisDirZ(xy*yz-xz*yy,xy*xz-yz*xx,detZ);
        scalar weightZ(detZ*detZ);
        if((weightedDir & axisDirZ) < 0.0)
            weightZ = -weightZ;
        weightedDir += axisDirZ * weightZ;

        if(mag(weightedDir) > SMALL)
        {
            normOk = true;
            normalVec = weightedDir/mag(weightedDir);
        }
        if (!normOk || mag(normalVec) < 1)
            normalVec = normalVector;

        // create best fitting plane
        plane bestFitPlane(contactCenter, normalVec);
        normalVec = bestFitPlane.normal();
        DynamicPointList commCellsPosInPlane;
        forAll (commCellsPositions,cell)
            commCellsPosInPlane.append(bestFitPlane.nearestPoint(commCellsPositions[cell]));

        vector q1(1.0, 0.0, 0.0);
        vector q2(0.0, 1.0, 0.0);
        if (abs(q1 & bestFitPlane.normal()) > abs(q2 & bestFitPlane.normal()))
            q1 = q2;

        vector u(bestFitPlane.normal() ^ q1);
        vector v(bestFitPlane.normal() ^ u);

        DynamicList<plane> clockwisePlanes;
        List<scalar> helpList(6);
        helpList[0] = 0.0;
        helpList[1] = 1.0;
        helpList[2] = 0.0;
        helpList[3] = -1.0;
        helpList[4] = 0.0;
        helpList[5] = 1.0;

        // loop over parts of plane to find and sort points
        DynamicVectorList commCellsInSections;
        for (label i = 0; i < 4; i++)
        {
            scalar uStep(helpList[i + 1] -  helpList[i]);
            scalar vStep(helpList[i + 2] -  helpList[i + 1]);
            DynamicPointList pointsInSection;
            for (scalar j = 0.0; j < 3.0; j += 1.0)
            {
                plane uPlane(contactCenter, u*(helpList[i] + uStep*j/4.0) + v*(helpList[i + 1] + vStep*j/4.0));
                plane vPlane(contactCenter, u*(helpList[i] + uStep*(j+1)/4.0) + v*(helpList[i + 1] + vStep*(j+1)/4.0));

                forAll (commCellsPosInPlane, celli)
                {
                    if (uPlane.sideOfPlane(commCellsPosInPlane[celli]) == 0 && vPlane.sideOfPlane(commCellsPosInPlane[celli]) == 1)
                        pointsInSection.append(commCellsPosInPlane[celli]);
                }

                if (pointsInSection.size() > SMALL)
                {
                    vector average(vector::zero);
                    forAll (pointsInSection, pointI)
                        average += pointsInSection[pointI];
                    average /= pointsInSection.size();

                    commCellsInSections.append(average);
                }
            }
        }
        // calculate contact area
        for (label i = 0; i + 1 < commCellsInSections.size(); i++)
        {
            vector AC(commCellsInSections[i] - contactCenter);
            vector BC(commCellsInSections[i + 1] - contactCenter);
            vector crossPr(AC ^ BC);
            area += mag(crossPr)/2;
        }

        if (commCellsInSections.size() > 2)
        {
            vector AC(commCellsInSections[commCellsInSections.size() - 1] - contactCenter);
            vector BC(commCellsInSections[0] - contactCenter);
            vector crossPr(AC ^ BC);
            area += mag(crossPr)/2;
        }
    }

    if (normalVec == vector::zero)
        normalVec = normalVector;

    reduce(area, sumOp<scalar>());
    area /= Pstream::nProcs();
    reduce(normalVec, sumOp<vector>());
    normalVec /= Pstream::nProcs();

    Tuple2<scalar,vector> returnValue(area,normalVec);

    return returnValue;
}
//---------------------------------------------------------------------------//
Tuple2<scalar,vector>  get2DcontactInfo(
    const dynamicFvMesh&   mesh,
    DynamicLabelList commonCells,
    vector normalVector,
    vector contactCenter,
    vector geometricD
)
{
    DynamicPointList commCellsPositions;

    forAll (commonCells, cCell)
    {
        commCellsPositions.append(mesh.C()[commonCells[cCell]]);
    }

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UOPstream send(proci, pBufs);
            send << commCellsPositions;
        }
    }

    pBufs.finishedSends();

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recv(proci, pBufs);
            DynamicPointList commCellsPositionsi (recv);
            commCellsPositions.append(commCellsPositionsi);
        }
    }

    // evaluate contact area
    scalar contactArea(0);

    scalar highestDistance(0);

    for (label i = 1; i < commCellsPositions.size(); i++)
    {
        if (mag(commCellsPositions[i-1] - commCellsPositions[i]) > highestDistance)
        {
            highestDistance = mag(commCellsPositions[i-1] - commCellsPositions[i]);
        }
    }

    reduce(highestDistance, maxOp<scalar>());
    if (commonCells.size() > 0)
        contactArea = sqrt(mag(mesh.Sf()[commonCells[0]])) * highestDistance;
    reduce(contactArea, maxOp<scalar>());

    vector geomDirections(geometricD);
    forAll (geomDirections, direction)
    {
        if (geomDirections[direction] == 1)
            geomDirections[direction] = 0;
        if (geomDirections[direction] == -1)
            geomDirections[direction] = 1;
    }

    vector normalVector2(vector::zero);

    forAll (commonCells, cCell)
    {
        vector tempor((mesh.C()[commonCells[cCell]] - contactCenter) ^ geomDirections);
        if ((tempor & normalVector) < 0)
            tempor *= -1;
        normalVector2 += tempor;
    }

    label numOfComCells(commonCells.size());
    sumReduce(normalVector2, numOfComCells);

    normalVector2 /= numOfComCells;

    Tuple2<scalar,vector> returnValue(contactArea,normalVector);

    if (mag(normalVector2) == 0) return returnValue;

    normalVector2 = normalVector2/mag(normalVector2);

    // assure that the normal vector points out of the body
    if ((normalVector2 & normalVector) < 0)
        normalVector2 *= -1;

    returnValue.second() = normalVector2;

    return returnValue;
}
//---------------------------------------------------------------------------//
template<>
void getPrtPrtContactVars<sphere,sphere>
(
    const dynamicFvMesh&   mesh,
    contactInfo& cInfo,
    contactInfo& tInfo,
    vector geometricD,
    prtPrtContactVars& prtPrtCntVars
)
{
    scalar pi = Foam::constant::mathematical::pi;

    bool case3D = true;
    label emptyDim = 0;
    // check if the case is 3D
    forAll (geometricD, direction)
    {
        if (geometricD[direction] == -1)
        {
            case3D = false;
            emptyDim = direction;
            break;
        }
    }

    vector cCenter(cInfo.getGeomModel().getCoM());
    vector tCenter(tInfo.getGeomModel().getCoM());
    scalar cRadius(cInfo.getGeomModel().getDC() / 2);
    scalar tRadius(tInfo.getGeomModel().getDC() / 2);

    Info << "Case 3D: " << case3D << endl;
    Info << "cCenter: " << cCenter << " cRadius " << cRadius << endl;
    Info << "tCenter: " << tCenter << " tRadius " << tRadius << endl;

    vector centerDir = cCenter - tCenter;
    Info << "centerDir " << centerDir << endl;
    scalar d = mag(centerDir);
    Info << "d " << d << endl;

    if(mag(centerDir) < SMALL || d > (cRadius + tRadius))
    {
        prtPrtCntVars.contactCenter_ = vector::zero;
        prtPrtCntVars.contactVolume_ = 0;
        prtPrtCntVars.contactNormal_ = vector::zero;
        prtPrtCntVars.contactArea_ = 0;
        return;
    }

    scalar xLength = (sqr(d) - sqr(cRadius) + sqr(tRadius))/(2*d);
    Info << "xLength " << xLength << endl;
    vector contactCenter = tCenter + (centerDir/d)*xLength;
    Info << "contactCenter " << contactCenter << endl;
    if(sqr(xLength) > sqr(tRadius))
    {
        prtPrtCntVars.contactCenter_ = contactCenter;
        prtPrtCntVars.contactVolume_ = (4/3)*pi*pow(tRadius,3);
        prtPrtCntVars.contactNormal_ = centerDir/d;
        prtPrtCntVars.contactArea_ = pi*sqr(tRadius);
        return;
    }
    scalar contactRad = sqrt(sqr(tRadius) - sqr(xLength));
    Info << "contactRad " << contactRad << endl;

    if(case3D)
    {
        scalar cSphCapH = cRadius - d + xLength;
        Info << "cSphCapH " << cSphCapH << endl;
        scalar tSphCapH = tRadius - xLength;
        Info << "tSphCapH " << tSphCapH << endl;

        scalar cSphCapV = (pi*sqr(cSphCapH)*(3*cRadius - cSphCapH)) / 3;
        Info << "cSphCapV " << cSphCapV  << endl;
        scalar tSphCapV = (pi*sqr(tSphCapH)*(3*tRadius - tSphCapH)) / 3;
        Info << "tSphCapV " << tSphCapV  << endl;

        prtPrtCntVars.contactCenter_ = contactCenter;
        prtPrtCntVars.contactVolume_ = cSphCapV + tSphCapV;
        prtPrtCntVars.contactNormal_ = centerDir/d;
        prtPrtCntVars.contactArea_ = pi*sqr(contactRad);
    }
    else
    {
        boundBox meshBounds = mesh.bounds();
        scalar emptyLength = meshBounds.max()[emptyDim]
                            - meshBounds.min()[emptyDim];

        scalar cCirSeg = sqr(cRadius)*acos((d - xLength)/cRadius)
                         - (d - xLength)*sqrt(sqr(cRadius) - sqr(d - xLength));
        scalar tCirSeg = sqr(tRadius)*acos(xLength/tRadius)
                         - xLength*sqrt(sqr(tRadius) - sqr(xLength));

        prtPrtCntVars.contactCenter_ = contactCenter;
        prtPrtCntVars.contactVolume_ = (cCirSeg + tCirSeg)*emptyLength;
        prtPrtCntVars.contactNormal_ = centerDir/d;
        prtPrtCntVars.contactArea_ = 2*contactRad*emptyLength;
    }
}
//---------------------------------------------------------------------------//
void solvePrtContact(
    const dynamicFvMesh&   mesh,
    contactInfo& cInfo,
    contactVars* cVars,
    contactInfo& tInfo,
    contactVars* tVars,
    vector geometricD,
    scalar deltaT,
    Tuple2<Tuple2<vector,vector>,Tuple2<vector,vector>>& outVars
)
{
    const contactType cT(cInfo.getGeomModel().getcType());
    const contactType tT(tInfo.getGeomModel().getcType());

    prtPrtContactVars prtPrtCntVars;

    if(cT == sphere && tT == sphere)
    {
        getPrtPrtContactVars<sphere,sphere>(
            mesh,
            cInfo,
            tInfo,
            geometricD,
            prtPrtCntVars
        );

        prtPrtContactVars prtPrtCntVarsTest;
        getPrtPrtContactVars<arbShape,arbShape>(
            mesh,
            cInfo,
            tInfo,
            geometricD,
            prtPrtCntVarsTest
        );
    }
    else
    {
        getPrtPrtContactVars<arbShape,arbShape>(
            mesh,
            cInfo,
            tInfo,
            geometricD,
            prtPrtCntVars
        );
    }

    if (prtPrtCntVars.contactVolume_ < SMALL)
    {
        return;
    }
    // get data from the current body necessary for computation of
    // contact forces
    vector cCoM(cInfo.getGeomModel().getCoM());                         //center of mass
    vector cVel(cVars->Vel_);                                           //linear velocity
    scalar cKN(cInfo.getkN());                                          //elastic normal stiffness
    scalar cGammaN(cInfo.getgammaN());                                  //normal viscosity
    scalar cKt(cInfo.getkt());                                          //elastic tangential stiffness
    scalar cGammat(cInfo.getgammat());                                  //tangential viscosity
    scalar cmu(cInfo.getmu());                                          //firction coef
    scalar cadhN(cInfo.getAdhN());                                      //adhesive force
    scalar cadhEqui(cInfo.getAdhEqui());                                //adhesive force
    scalar cM(cInfo.getGeomModel().getM0());                            //mass
    scalar cRhoS(cVars->rhoS_.value());                                 //viscosity
    vector cAxis(cVars->Axis_);                                         //Axis
    scalar cOmega(cVars->omega_);                                       //Omega

    // get data from the neighbor body
    vector tCoM(tInfo.getGeomModel().getCoM());                         //center of mass
    vector tVel(tVars->Vel_);                                           //linear velocity
    scalar tKN(tInfo.getkN());                                          //elastic normal stiffness
    scalar tGammaN(tInfo.getgammaN());                                  //normal viscosity
    scalar tKt(tInfo.getkt());                                          //elastic tangential stiffness
    scalar tGammat(tInfo.getgammat());                                  //tangential viscosity
    scalar tmu(tInfo.getmu());                                          //firction coef
    scalar tadhN(tInfo.getAdhN());                                      //adhesive force
    scalar tadhEqui(tInfo.getAdhEqui());                                //adhesive force
    scalar tM(tInfo.getGeomModel().getM0());                            //mass
    scalar tRhoS(tVars->rhoS_.value());                                 //viscosity
    vector tAxis(tVars->Axis_);                                         //Axis
    scalar tOmega(tVars->omega_);                                       //Omega


    vector  FN(vector::zero);                                           //placeholder for normal force
    vector  Ft(vector::zero);                                           //placeholder for tangential force
    Info << "-- Detected particle-particle contact " << cVars->bodyId_ << " : " << tVars->bodyId_ << endl;
    Info << "-- Particle-particle contact center " << prtPrtCntVars.contactCenter_ << endl;
    Info << "-- Particle-particle contact normal " << prtPrtCntVars.contactNormal_ << endl;
    Info << "-- Particle-particle contact volume " << prtPrtCntVars.contactVolume_ << endl;
    Info << "-- Particle-particle contact area "   << prtPrtCntVars.contactArea_ << endl;

    // compute mean model parameters
    scalar aKN((cKN*tKN)/(cKN+tKN+SMALL));
    scalar aGammaN((cGammaN*tGammaN)/(cGammaN+tGammaN+SMALL));
    scalar aKt((cKt*tKt)/(cKt+tKt+SMALL));
    scalar aGammat((cGammat*tGammat)/(cGammat+tGammat+SMALL));
    scalar amu((cmu*tmu)/(cmu+tmu+SMALL));
    scalar aadhN((cadhN*tadhN)/(cadhN+tadhN+SMALL));

    vector cLVec(prtPrtCntVars.contactCenter_-cCoM);
    vector tLVec(prtPrtCntVars.contactCenter_-tCoM);

    // compute normal to movement and relative velocity
    vector nVec(prtPrtCntVars.contactNormal_);

    vector cplanarVec       =  cLVec- cAxis*((cLVec)&cAxis);
    vector cVeli(-(cplanarVec^cAxis)*cOmega + cVel);

    vector tplanarVec       =  tLVec- tAxis*((tLVec)&tAxis);
    vector tVeli(-(tplanarVec^tAxis)*tOmega + tVel);

    scalar Vn(-(cVeli - tVeli) & nVec);
    scalar Lc(4*mag(cLVec)*mag(tLVec)/(mag(cLVec)+mag(tLVec)));

    scalar reduceM(cM*tM/(cM+tM));

    // compute the normal force
    FN = (aKN*prtPrtCntVars.contactVolume_/(Lc+SMALL) + aGammaN*sqrt(aKN*reduceM/pow(Lc+SMALL,3))*(prtPrtCntVars.contactArea_ * Vn))*nVec;
    // compute adhesive force
    scalar FAc(aadhN*prtPrtCntVars.contactArea_);
    scalar FAeq(min(aKN*((cadhEqui*cM)/(cRhoS + SMALL))/(Lc+SMALL), aKN*((tadhEqui*tM)/(tRhoS + SMALL))/(Lc+SMALL)));
    scalar partMul(max(prtPrtCntVars.contactVolume_ * cRhoS / (cM+SMALL) / (cadhEqui+SMALL), prtPrtCntVars.contactVolume_ * tRhoS / (tM+SMALL) / (tadhEqui+SMALL)));
    if(partMul > 1)
    {
        partMul = 1;
    }
    vector FA((FAeq * partMul  + FAc * (1-partMul)) * nVec);

    vector cFtLast(vector::zero);
    vector tFtLast(vector::zero);

    // find history of tangential force between these two particles
    DynamicList<Tuple2<label,Tuple2<label,vector>>>& chistoryFt(cInfo.getHistoryhistoryFt());
    DynamicList<Tuple2<label,Tuple2<label,vector>>>& thistoryFt(tInfo.getHistoryhistoryFt());

    bool cFtLastFinded(false);
    forAll (chistoryFt,cFti)
    {
        if (chistoryFt[cFti].first() == tVars->bodyId_)
        {
            cFtLastFinded = true;
            cFtLast = chistoryFt[cFti].second().second();
            break;
        }
    }

    bool tFtLastFinded(false);
    forAll (thistoryFt,tFti)
    {
        if (thistoryFt[tFti].first() == cVars->bodyId_)
        {
            tFtLastFinded = true;
            tFtLast = thistoryFt[tFti].second().second();
            break;
        }
    }

    // magnitude of last Ft should be same but only oposite
    vector FtLast((cFtLast - tFtLast)/2);
    // project last Ft into a new direction
    vector FtLastP(FtLast - (FtLast & nVec) * nVec);
    // scale projected Ft to have same magnitude as FtLast
    vector FtLastr(mag(FtLast) * (FtLastP/(mag(FtLastP)+SMALL)));
    // compute relative tangential velocity
    vector cVeliNorm((cVeli & nVec)*nVec);
    cVeli -= cVeliNorm;

    vector tVeliNorm((cVeli & nVec)*nVec);
    tVeli -= tVeliNorm;

    vector Vt(cVeli - tVeli);
    // compute tangential force
    vector Ftdi(- aGammat*sqrt(aKN*reduceM*Lc)*Vt);
    Ft = (FtLastr - aKt*Vt*deltaT + Ftdi);

    if (mag(Ft) > amu * mag(FN))
    {
        Ft *= amu * mag(FN) / mag(Ft);
    }

    FN -= FA;

    vector F(FN+Ft);

    // update history of tangential force
    if (cFtLastFinded)
    {
        forAll (chistoryFt,cFti)
        {
            if (chistoryFt[cFti].first() == tVars->bodyId_)
            {
                Tuple2<label,vector> help(1,Ft);
                chistoryFt[cFti].second() = help;
                break;
            }
        }
    }
    else
    {
        Tuple2<label,vector> help(1,Ft);
        Tuple2<label,Tuple2<label,vector>> help2(tVars->bodyId_, help);

        chistoryFt.append(help2);
    }

    if (tFtLastFinded)
    {
        forAll (thistoryFt,tFti)
        {
            if (thistoryFt[tFti].first() == cVars->bodyId_)
            {
                Tuple2<label,vector> help(1,-Ft);
                thistoryFt[tFti].second() = help;
                break;
            }
        }
    }
    else
    {
        Tuple2<label,vector> help(1,-Ft);
        Tuple2<label,Tuple2<label,vector>> help2(cVars->bodyId_, help);

        thistoryFt.append(help2);
    }

    // add the computed force to the affected bodies
    vector cTN(cLVec ^  F);
    vector tTN(tLVec ^ -F);
    outVars.first().first() = F;
    outVars.first().second() = cTN;
    outVars.second().first() = -F;
    outVars.second().second() = tTN;
}
//---------------------------------------------------------------------------//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace contactModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
