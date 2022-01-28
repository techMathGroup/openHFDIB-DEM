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
#include "interpolationInfo.H"

using namespace Foam;

autoPtr<vectorField> interpolationInfo::surfNorm_;

//---------------------------------------------------------------------------//
interpolationInfo::interpolationInfo
(
    const Foam::fvMesh& mesh,
    const volScalarField& body,
    List<DynamicLabelList>& surfCells,
    bool sdBasedLambda,
    scalar intSpan
)
:
mesh_(mesh),
body_(body),
surfCells_(surfCells),
sdBasedLambda_(sdBasedLambda),
intSpan_(intSpan)
{}

interpolationInfo::~interpolationInfo()
{}
//---------------------------------------------------------------------------//
point interpolationInfo::getIbPoint
(
    label scell,
    vectorField& surfNorm,
    const volScalarField& body
)
{
    if (sdBasedLambda_)
    {
        scalar minMaxBody(max(min(body[scell],1.0-SMALL),SMALL));
        scalar intDist = Foam::atanh(2.0*minMaxBody - 1.0)
                            *Foam::pow(mesh_.V()[scell],0.333)/intSpan_;
        return (mesh_.C()[scell] + surfNorm[scell]*intDist);
    }
    else
    {
        labelList cellNb(mesh_.cellCells()[scell]); //list of neighbours
        scalar intDist = 0;
        forAll (cellNb,nbCellI)
        {
            intDist += mag(mesh_.C()[cellNb[nbCellI]] - mesh_.C()[scell]);
        }
        intDist /= scalar(cellNb.size());
        return (mesh_.C()[scell] - surfNorm[scell]*(0.5-body[scell])*intDist);
    }
}
//---------------------------------------------------------------------------//
