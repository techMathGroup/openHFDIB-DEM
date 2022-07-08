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
#include "prtContact.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace contactModel
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//---------------------------------------------------------------------------//
bool detectPrtPrtContact(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass
)
{
    if
    (
        cClass.getGeomModel().getcType() == sphere
        &&
        tClass.getGeomModel().getcType() == sphere
    )
    {
        return detectPrtPrtContact_Sphere
        (
            mesh,
            cClass,
            tClass
        );
    }
    else if
    (
        cClass.getGeomModel().getcType() == cluster
        ||
        tClass.getGeomModel().getcType() == cluster
    )
    {
        return detectPrtPrtContact_Cluster
        (
            mesh,
            cClass,
            tClass
        );
    }
    else
    {
        return detectPrtPrtContact_ArbShape
        (
            mesh,
            cClass,
            tClass
        );
    }
}
//---------------------------------------------------------------------------//
bool detectPrtPrtContact_ArbShape(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass
)
{
    List<DynamicLabelList>    commonCells;
    commonCells.setSize(Pstream::nProcs());

    List<DynamicLabelList> cSurfCells(cClass.getSurfCells());
    List<DynamicLabelList> cIntCells(cClass.getIntCells());
    List<DynamicLabelList> tSurfCells(tClass.getSurfCells());

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

        boolList tcenterInsideList = cClass.getGeomModel().pointInside(mesh.C());
        boolList ccenterInsideList = tClass.getGeomModel().pointInside(mesh.C());

        // iterate over all surfCells and intCells to evaluate intersected volume
        forAll (cSurfCells[Pstream::myProcNo()],cSCellI)
        {
            const labelList& vertexLabels = mesh.cellPoints()[cSurfCells[Pstream::myProcNo()][cSCellI]];
            const pointField vertexPoints(pp,vertexLabels);
            boolList cvertexesInside = cClass.getGeomModel().pointInside( vertexPoints );
            boolList tvertexesInside = tClass.getGeomModel().pointInside( vertexPoints );
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
            boolList cvertexesInside = cClass.getGeomModel().pointInside( vertexPoints );
            boolList tvertexesInside = tClass.getGeomModel().pointInside( vertexPoints );
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
bool detectPrtPrtContact_Sphere
(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass
)
{
    if
    (
        mag(cClass.getGeomModel().getCoM()-tClass.getGeomModel().getCoM())
        <
        ((cClass.getGeomModel().getDC() / 2) + (tClass.getGeomModel().getDC() / 2))
    )
    {
        Info << "There is prt prt contact" << endl;
        return true;
    }
    Info << "There is NOT prt prt contact" << endl;
    return false;
}
//---------------------------------------------------------------------------//
bool detectPrtPrtContact_Cluster
(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass
)
{
    PtrList<geomModel> cBodies(0);
    PtrList<geomModel> tBodies(0);

    if(cClass.getGeomModel().isCluster())
    {
        clusterBody& cCluster = dynamic_cast<clusterBody&>(cClass.getGeomModel());
        PtrList<geomModel>& cBodiesR = cCluster.getClusterBodies();
        forAll(cBodiesR, cIbI)
        {
            cBodies.append(cBodiesR[cIbI].getGeomModel());
        }
    }
    else
    {
        cBodies.append(cClass.getGeomModel().getGeomModel());
    }

    if(tClass.getGeomModel().isCluster())
    {
        clusterBody& tCluster = dynamic_cast<clusterBody&>(tClass.getGeomModel());
        PtrList<geomModel>& tBodiesR = tCluster.getClusterBodies();
        forAll(tBodiesR, tIbI)
        {
            tBodies.append(tBodiesR[tIbI].getGeomModel());
        }
    }
    else
    {
        tBodies.append(tClass.getGeomModel().getGeomModel());
    }

    forAll(cBodies, cIbI)
    {
        forAll(tBodies, tIbI)
        {
            autoPtr<geomModel> cGeomModel(cBodies[cIbI].getGeomModel());
            autoPtr<ibContactClass> cIbClassI(new ibContactClass(
                cGeomModel,
                cClass.getMatInfo()
            ));

            autoPtr<geomModel> tGeomModel(tBodies[tIbI].getGeomModel());
            autoPtr<ibContactClass> tIbClassI(new ibContactClass(
                tGeomModel,
                tClass.getMatInfo()
            ));

            if (detectPrtPrtContact(
                mesh,
                cIbClassI(),
                tIbClassI()
            ))
            {
                return true;
            }
        }
    }
    return false;
}
//---------------------------------------------------------------------------//
void getPrtContactVars_ArbShape(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass,
    prtContactVars& prtCntVars
)
{
    List<DynamicLabelList>    commonCells;
    commonCells.setSize(Pstream::nProcs());

    List<DynamicLabelList> cSurfCells(cClass.getSurfCells());
    List<DynamicLabelList> cIntCells(cClass.getIntCells());
    List<DynamicLabelList> tSurfCells(tClass.getSurfCells());

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
        prtCntVars.contactCenter_ = vector::zero;
        prtCntVars.contactVolume_ = 0;
        prtCntVars.contactNormal_ = vector::zero;
        prtCntVars.contactArea_ = 0;
        return;
    }

    contactCenter /= numOfComCells;

    scalar tDC(tClass.getGeomModel().getDC());                           //characteristic diameter

    point closestPoint = vector::zero;
    vector normalVector = vector::zero;
    tClass.getGeomModel().getClosestPointAndNormal
    (
        contactCenter,
        vector::one * tDC,
        closestPoint,
        normalVector
    );

    scalar contactArea(0);
    // use edge cells to find contact area (better precision then surfcells)
    DynamicVectorList edgePoints;
    DynamicLabelList edgeCells;

    scalar intersectedVolume(0);

    if (commonCells[Pstream::myProcNo()].size() > SMALL)
    {
        const pointField& pp = mesh.points();

        boolList tcenterInsideList = tClass.getGeomModel().pointInside( mesh.C());
        boolList ccenterInsideList = cClass.getGeomModel().pointInside( mesh.C());

        forAll (cSurfCells[Pstream::myProcNo()],cSCellI)
        {
            bool partiallyInT(false);

            const labelList& vertexLabels = mesh.cellPoints()[cSurfCells[Pstream::myProcNo()][cSCellI]];
            const pointField vertexPoints(pp,vertexLabels);
            boolList tvertexesInside = tClass.getGeomModel().pointInside( vertexPoints );
            boolList cvertexesInside = cClass.getGeomModel().pointInside( vertexPoints );
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
            boolList tvertexesInside = tClass.getGeomModel().pointInside( vertexPoints );
            boolList cvertexesInside = cClass.getGeomModel().pointInside( vertexPoints );
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
            Tuple2<scalar,vector> returnTuple = get3DcontactVars(mesh, edgeCells, edgePoints, normalVector, contactCenter, cClass.getGeomModel().getOwner());
            contactArea = returnTuple.first();
            normalVector = returnTuple.second();
        }
        else
        {
            Tuple2<scalar,vector> returnTuple = get2DcontactVars(mesh, commonCells[Pstream::myProcNo()], normalVector, contactCenter);
            contactArea = returnTuple.first();
            normalVector = returnTuple.second();
        }
    }

    if (intersectedVolume > 0 && contactArea > 0)
    {
        prtCntVars.contactCenter_ = contactCenter;
        prtCntVars.contactVolume_ = intersectedVolume;
        prtCntVars.contactNormal_ = normalVector;
        prtCntVars.contactArea_ = contactArea;
    }
    else
    {
        prtCntVars.contactCenter_ = vector::zero;
        prtCntVars.contactVolume_ = 0;
        prtCntVars.contactNormal_ = vector::zero;
        prtCntVars.contactArea_ = 0;
    }
}
//---------------------------------------------------------------------------//
Tuple2<scalar,vector> get3DcontactVars(
    const fvMesh&   mesh,
    DynamicLabelList& commonCells,
    DynamicVectorList& edgePoints,
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
Tuple2<scalar,vector>  get2DcontactVars(
    const fvMesh&   mesh,
    DynamicLabelList& commonCells,
    vector normalVector,
    vector contactCenter
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

    vector normalVector2(vector::zero);

    forAll (commonCells, cCell)
    {
        vector tempor((mesh.C()[commonCells[cCell]] - contactCenter) ^ emptyDir);
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
void getPrtContactVars_Sphere
(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass,
    prtContactVars& prtCntVars
)
{
    scalar cRadius = cClass.getGeomModel().getDC() / 2;
    scalar tRadius = tClass.getGeomModel().getDC() / 2;
    vector tCenter = tClass.getGeomModel().getCoM();

    vector centerDir = cClass.getGeomModel().getCoM()
                        - tCenter;

    scalar d = mag(centerDir);

    if(mag(centerDir) < SMALL || d > (cRadius + tRadius))
    {
        prtCntVars.contactCenter_ = vector::zero;
        prtCntVars.contactVolume_ = 0;
        prtCntVars.contactNormal_ = vector::zero;
        prtCntVars.contactArea_ = 0;
        return;
    }

    scalar xLength = (sqr(d) - sqr(cRadius) + sqr(tRadius))/(2*d);

    if(sqr(xLength) > sqr(tRadius))
    {
        prtCntVars.contactCenter_ = tCenter + (centerDir/d)*xLength;
        prtCntVars.contactVolume_ = (4/3)*Foam::constant::mathematical::pi
            *pow(tRadius,3);
        prtCntVars.contactNormal_ = centerDir/d;
        prtCntVars.contactArea_ = Foam::constant::mathematical::pi
            *sqr(tRadius);
        return;
    }

    if(case3D)
    {

        scalar cSphCapV = (Foam::constant::mathematical::pi
                            *sqr(cRadius - d + xLength)
                            *(3*cRadius - (cRadius - d + xLength))) / 3;
        scalar tSphCapV = (Foam::constant::mathematical::pi
                            *sqr(tRadius - xLength)
                            *(3*tRadius - (tRadius - xLength))) / 3;

        prtCntVars.contactCenter_ = tCenter + (centerDir/d)*xLength;
        prtCntVars.contactVolume_ = cSphCapV + tSphCapV;
        prtCntVars.contactNormal_ = centerDir/d;
        prtCntVars.contactArea_ = Foam::constant::mathematical::pi
                                    *sqr(sqrt(sqr(tRadius) - sqr(xLength)));
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

        prtCntVars.contactCenter_ = tCenter + (centerDir/d)*xLength;
        prtCntVars.contactVolume_ = (cCirSeg + tCirSeg)*emptyLength;
        prtCntVars.contactNormal_ = centerDir/d;
        prtCntVars.contactArea_ = 2*sqrt(sqr(tRadius)
                                        - sqr(xLength))*emptyLength;
    }
}
//---------------------------------------------------------------------------//
void getPrtContactVars_Cluster
(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass,
    prtContactVars& prtCntVars
)
{
    prtCntVars.contactCenter_ = vector::zero;
    prtCntVars.contactVolume_ = 0;
    prtCntVars.contactNormal_ = vector::zero;
    prtCntVars.contactArea_ = 0;

    PtrList<geomModel> cBodies(0);
    PtrList<geomModel> tBodies(0);

    if(cClass.getGeomModel().isCluster())
    {
        clusterBody& cCluster = dynamic_cast<clusterBody&>(cClass.getGeomModel());
        PtrList<geomModel>& cBodiesR = cCluster.getClusterBodies();
        forAll(cBodiesR, cIbI)
        {
            cBodies.append(cBodiesR[cIbI].getGeomModel());
        }
    }
    else
    {
        cBodies.append(cClass.getGeomModel().getGeomModel());
    }

    if(tClass.getGeomModel().isCluster())
    {
        clusterBody& tCluster = dynamic_cast<clusterBody&>(tClass.getGeomModel());
        PtrList<geomModel>& tBodiesR = tCluster.getClusterBodies();
        forAll(tBodiesR, tIbI)
        {
            tBodies.append(tBodiesR[tIbI].getGeomModel());
        }
    }
    else
    {
        tBodies.append(tClass.getGeomModel().getGeomModel());
    }

    forAll(cBodies, cIbI)
    {
        forAll(tBodies, tIbI)
        {
            autoPtr<geomModel> cGeomModel(cBodies[cIbI].getGeomModel());
            autoPtr<ibContactClass> cIbClassI(new ibContactClass(
                cGeomModel,
                cClass.getMatInfo()
            ));

            autoPtr<geomModel> tGeomModel(tBodies[tIbI].getGeomModel());
            autoPtr<ibContactClass> tIbClassI(new ibContactClass(
                tGeomModel,
                tClass.getMatInfo()
            ));

            prtContactVars prtCntVarsI;

            getPrtContactVars(
                mesh,
                cIbClassI(),
                tIbClassI(),
                prtCntVarsI
            );

            if (prtCntVarsI.contactVolume_ > prtCntVars.contactVolume_)
            {
                prtCntVars = prtCntVarsI;
            }
        }
    }
}
//---------------------------------------------------------------------------//
void getPrtContactVars
(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass,
    prtContactVars& prtCntVars
)
{
    if
    (
        cClass.getGeomModel().getcType() == sphere
        &&
        tClass.getGeomModel().getcType() == sphere
    )
    {
        getPrtContactVars_Sphere(
            mesh,
            cClass,
            tClass,
            prtCntVars
        );
    }
    else if
    (
        cClass.getGeomModel().getcType() == cluster
        ||
        tClass.getGeomModel().getcType() == cluster
    )
    {
        getPrtContactVars_Cluster(
            mesh,
            cClass,
            tClass,
            prtCntVars
        );
    }
    else
    {
        getPrtContactVars_ArbShape(
            mesh,
            cClass,
            tClass,
            prtCntVars
        );
    }
}
//---------------------------------------------------------------------------//
bool solvePrtContact(
    const fvMesh&   mesh,
    prtContactInfo& cInfo,
    scalar deltaT
)
{
    getPrtContactVars(
        mesh,
        cInfo.getcClass(),
        cInfo.gettClass(),
        cInfo.getprtCntVars()
    );

    if (!(cInfo.getprtCntVars().contactVolume_ > 0))
    {
        return false;
    }

    InfoH << DEM_Info << "-- Detected Particle-particle contact: -- body "
            << cInfo.getcPair().first() << " & -- body "
            << cInfo.getcPair().second() << endl;
    InfoH << DEM_Info << "-- Particle-particle contact cBody pos: "
            << cInfo.getcClass().getGeomModel().getCoM() << " & tBody pos: "
            << cInfo.gettClass().getGeomModel().getCoM() << endl;
    InfoH << DEM_Info << "-- Particle-particle contact center "
            << cInfo.getprtCntVars().contactCenter_ << endl;
    InfoH << DEM_Info << "-- Particle-particle contact normal "
            << cInfo.getprtCntVars().contactNormal_ << endl;
    InfoH << DEM_Info << "-- Particle-particle contact volume "
            << cInfo.getprtCntVars().contactVolume_ << endl;
    InfoH << DEM_Info << "-- Particle-particle contact area "
            << cInfo.getprtCntVars().contactArea_ << endl;

    cInfo.evalVariables();

    // compute the normal force
    vector F = cInfo.getFNe();
    InfoH << DEM_Info << "-- Particle-particle contact FNe " << F << endl;

    vector FNd = cInfo.getFNd();
    InfoH << DEM_Info << "-- Particle-particle contact FNd " << FNd << endl;

    F += FNd;
    InfoH << DEM_Info << "-- Particle-particle contact FN " << F << endl;

    vector Ft = cInfo.getFt(deltaT);
    InfoH << DEM_Info << "-- Particle-particle contact Ft " << Ft << endl;

    if (mag(Ft) > cInfo.getMu() * mag(F))
    {
        Ft *= cInfo.getMu() * mag(F) / mag(Ft);
    }
    InfoH << DEM_Info << "-- Particle-particle contact Ft clamped " << Ft << endl;
    F += Ft;    

    vector FA = cInfo.getFA();
    InfoH << DEM_Info << "-- Particle-particle contact FA " << FA << endl;
    F -= FA;

    InfoH << DEM_Info << "-- Resolved Particle-particle contact: -- body "
            << cInfo.getcPair().first() << " & -- body "
            << cInfo.getcPair().second() << endl;

    // add the computed force to the affected bodies
    cInfo.getOutForce().first().F = F;
    cInfo.getOutForce().first().T = cInfo.getcLVec() ^  F;
    cInfo.getOutForce().second().F = -F;
    cInfo.getOutForce().second().T = cInfo.gettLVec() ^ -F;
    return true;
}
//---------------------------------------------------------------------------//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace contactModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
