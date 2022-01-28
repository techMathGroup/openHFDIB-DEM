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
#include "lineInt.H"

#define ORDER 2

using namespace Foam;

//---------------------------------------------------------------------------//
lineInt::lineInt(dictionary& interpDict):
interpDict_(interpDict)
{}
lineInt::~lineInt()
{}
//---------------------------------------------------------------------------//
void lineInt::ibInterpolate
(
    interpolationInfo& intpInfo,
    volVectorField& Ui,
    vectorField ibPointsVal,
    const Foam::fvMesh& mesh
)
{
    lineIntInfo& lsInfo
        = dynamic_cast<lineIntInfo&>(intpInfo);

    correctVelocity
    (
        lsInfo,
        Ui,
        ibPointsVal,
        mesh
    );
}
//---------------------------------------------------------------------------//
void lineInt::correctVelocity
(
    lineIntInfo& intpInfo,
    volVectorField& Ui,
    vectorField& ibPointsVal,
    const Foam::fvMesh& mesh
)
{
    const DynamicLabelList& cSurfCells = intpInfo.getSurfCells();

    List<point>& ibPoints = intpInfo.getIbPoints();
    List<List<intPoint>>& intPoints = intpInfo.getIntPoints();

    getCurVelocity(intPoints);
    List<label> intOrder = getIntOrder(intPoints);

    forAll(intPoints, ibp)
    {
        label cellI = cSurfCells[ibp];

        switch(intOrder[ibp])
        {
            case 0:
            {
                Ui[cellI] = ibPointsVal[ibp];
                break;
            }

            case 1:
            {
                vector VP1 = intPoints[ibp][0].iVel_ - ibPointsVal[ibp];

                // distance between interpolation points
                scalar deltaR = mag(intPoints[ibp][0].iPoint_ - ibPoints[ibp]);

                // cell center to surface distance
                scalar ds = mag(mesh.C()[cellI] - ibPoints[ibp]);

                vector linCoeff = VP1/(deltaR+SMALL);

                Ui[cellI] = linCoeff*ds + ibPointsVal[ibp];
                break;
            }

            case 2:
            {
                vector VP1 =  intPoints[ibp][0].iVel_ - ibPointsVal[ibp];
                vector VP2 =  intPoints[ibp][1].iVel_ - ibPointsVal[ibp];

                // distance between interpolation points
                scalar deltaR1 = mag(intPoints[ibp][1].iPoint_
                        - intPoints[ibp][0].iPoint_);
                scalar deltaR2 = mag(intPoints[ibp][0].iPoint_ - ibPoints[ibp]);

                // cell center to surface distance
                scalar ds = mag(mesh.C()[cellI] - ibPoints[ibp]);

                vector quadCoeff = (VP2 - VP1)*deltaR1 - VP1*deltaR2;
                quadCoeff       /= (deltaR1*deltaR2*(deltaR1 + deltaR2)+SMALL);

                vector linCoeff  = (VP1-VP2)*Foam::pow(deltaR1,2.0);
                linCoeff        += 2.0*VP1*deltaR1*deltaR2;
                linCoeff        += VP1*Foam::pow(deltaR2,2.0);
                linCoeff        /= (deltaR1*deltaR2*(deltaR1 + deltaR2)+SMALL);

                Ui[cellI] = quadCoeff*ds*ds + linCoeff*ds + ibPointsVal[ibp];
                break;
            }
        }
    }
}
//---------------------------------------------------------------------------//
void lineInt::getCurVelocity
(
    List<List<intPoint>>& intPoints
)
{
    List<DynamicPointList> intPointToSync(Pstream::nProcs());
    List<DynamicLabelList> cellLabelToSync(Pstream::nProcs());
    List<DynamicList<Tuple2<label,label>>> indexesToSync(Pstream::nProcs());

    forAll(intPoints, ibCellI)
    {
        forAll(intPoints[ibCellI], iPoint)
        {
            intPoint& curIPoint = intPoints[ibCellI][iPoint];

            if(curIPoint.iProc_ == Pstream::myProcNo())
            {
                curIPoint.iVel_ =  interpV_->interpolate(
                    curIPoint.iPoint_,
                    curIPoint.iCell_
                );
            }
            else
            {
                if(curIPoint.iProc_ != -1)
                {
                    intPointToSync[curIPoint.iProc_].append(curIPoint.iPoint_);
                    cellLabelToSync[curIPoint.iProc_].append(curIPoint.iCell_);
                    indexesToSync[curIPoint.iProc_].append(
                        Tuple2<label,label>(ibCellI, iPoint)
                    );
                }
            }
        }
    }

    PstreamBuffers pBufsIntP(Pstream::commsTypes::nonBlocking);
    PstreamBuffers pBufsIntC(Pstream::commsTypes::nonBlocking);

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if(proci != Pstream::myProcNo())
        {
            UOPstream sendIntP(proci, pBufsIntP);
            UOPstream sendIntC(proci, pBufsIntC);

            sendIntP << intPointToSync[proci];
            sendIntC << cellLabelToSync[proci];
        }
    }

    pBufsIntP.finishedSends();
    pBufsIntC.finishedSends();

    List<DynamicPointList> intPRecv(Pstream::nProcs());
    List<DynamicLabelList> intCRecv(Pstream::nProcs());

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recvIntP(proci, pBufsIntP);
            UIPstream recvIntC(proci, pBufsIntC);

            DynamicPointList recIntP (recvIntP);
            DynamicLabelList recIntC (recvIntC);

            intPRecv[proci] = recIntP;
            intCRecv[proci] = recIntC;
        }
    }

    pBufsIntP.clear();
    pBufsIntC.clear();

    List<DynamicVectorList> intVelRtrn(Pstream::nProcs());

    forAll(intPRecv, proci)
    {
        forAll(intPRecv[proci], pi)
        {
            intVelRtrn[proci].append(
                interpV_->interpolate(
                    intPRecv[proci][pi],
                    intCRecv[proci][pi]
            ));
        }
    }

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if(proci != Pstream::myProcNo())
        {
            UOPstream sendIntV(proci, pBufsIntP);

            sendIntV << intVelRtrn[proci];
        }
    }

    pBufsIntP.finishedSends();

    List<DynamicVectorList> intVelRcv(Pstream::nProcs());

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recvIntV(proci, pBufsIntP);

            DynamicVectorList recIntV (recvIntV);

            intVelRcv[proci] = recIntV;
        }
    }

    pBufsIntP.clear();

    forAll(intVelRcv, proci)
    {
        forAll(intVelRcv[proci], pi)
        {
            Tuple2<label, label> cInd = indexesToSync[proci][pi];
            intPoints[cInd.first()][cInd.second()].iVel_ = intVelRcv[proci][pi];
        }
    }
}
//---------------------------------------------------------------------------//
List<label> lineInt::getIntOrder
(
    List<List<intPoint>>& intPoints
)
{
    List<label> intOrderToRtn(intPoints.size(), 0);

    forAll(intPoints, ibp)
    {
        forAll(intPoints[ibp], ip)
        {
            if(intPoints[ibp][ip].iProc_ != -1)
            {
                ++intOrderToRtn[ibp];
            }
        }
    }

    return intOrderToRtn;
}
//---------------------------------------------------------------------------//
