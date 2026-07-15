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
#include "leastSquaresInt.H"
#include "scalarMatrices.H"

using namespace Foam;

//---------------------------------------------------------------------------//
leastSquaresInt::leastSquaresInt
(
    const fvMesh& mesh,
    scalar distFactor,
    scalar radiusFactor,
    scalar angleFactor,
    scalar maxCCRows
)
:
mesh_(mesh),
distFactor_(distFactor),
radiusFactor_(radiusFactor),
angleFactor_(angleFactor),
maxCCRows_(maxCCRows)
{}
leastSquaresInt::~leastSquaresInt()
{}
//---------------------------------------------------------------------------//
void leastSquaresInt::ibInterpolate
(
    interpolationInfo& intpInfo,
    volVectorField& Ui,
    vectorField ibPointsVal,
    const Foam::fvMesh& mesh
)
{
    leastSquaresIntInfo& lsInfo
        = dynamic_cast<leastSquaresIntInfo&>(intpInfo);

    imposeDirichletCondition
    (
        lsInfo,
        Ui,
        ibPointsVal,
        mesh
    );
}
//---------------------------------------------------------------------------//
void leastSquaresInt::imposeDirichletCondition
(
    leastSquaresIntInfo& intpInfo,
    volVectorField& Ui,
    vectorField& ibPointsVal,
    const Foam::fvMesh& mesh
)
{
    const DynamicLabelList& ibCells = intpInfo.getSurfCells();
    const List<point>& ibPoints = intpInfo.getIbPoints();
    const List<vector>& ibNormals = intpInfo.getIbNormals();
    const labelListList& cellCells = intpInfo.getCellCells();
    const PtrList<scalarRectangularMatrix>& cInvDirMats =
        intpInfo.getInvDirMats();

    autoPtr<vectorField> tpolyU(new vectorField(Ui.internalField(),ibCells));
    vectorField& polyU = tpolyU();

    const vectorField& C = mesh_.cellCentres();

    label nCoeffs = 5;
    if(case3D)
    {
        nCoeffs += 4;
    }

    label counter = 0;
    scalarField error(ibCells.size(), 0);
    scalar maxError = 0;

    do
    {
        counter++;
        scalar maxMagPolyU = 0;
        forAll(ibCells, cellI)
        {
            label curCell = ibCells[cellI];
            const labelList& curCells = cellCells[cellI];
            const scalarRectangularMatrix& curInvMatrix = cInvDirMats[cellI];

            if(curCells.size() < nCoeffs)
            {
                Info << "Not eneough interp Points" << endl;
            }
            else
            {
                vectorField coeffs(nCoeffs, vector::zero);
                vectorField source(curCells.size(), vector::zero);

                for(label i = 0; i < curCells.size(); i++)
                {
                    source[i] = Ui[curCells[i]] - ibPointsVal[cellI];
                }

                for(label i = 0; i < nCoeffs; i++)
                {
                    for(label j = 0; j < source.size(); j++)
                    {
                        coeffs[i] += curInvMatrix[i][j]*source[j];
                    }
                }

                vector oldPolyU = polyU[cellI];
                vector R = C[curCell] - ibPoints[cellI];
                if((R & ibNormals[cellI]) < 0)
                    R = vector::zero;

                if(case3D)
                {
                    polyU[cellI] = ibPointsVal[cellI]
                        + coeffs[0]*R.x()
                        + coeffs[1]*R.y()
                        + coeffs[2]*R.x()*R.y()
                        + coeffs[3]*sqr(R.x())
                        + coeffs[4]*sqr(R.y())
                        + coeffs[5]*R.z()
                        + coeffs[6]*R.x()*R.z()
                        + coeffs[7]*R.y()*R.z()
                        + coeffs[8]*sqr(R.z());
                }
                else
                {
                    labelList A (2,0.0);
                    label validDim = 0;
                    for(label dim = 0; dim < 3; dim++)
                    {
                        if(dim != emptyDim)
                        {
                            A[validDim++] = dim;
                        }
                    }
                    polyU[cellI] = ibPointsVal[cellI]
                        + coeffs[0]*R[A[0]]
                        + coeffs[1]*R[A[1]]
                        + coeffs[2]*R[A[0]]*R[A[1]]
                        + coeffs[3]*sqr(R[A[0]])
                        + coeffs[4]*sqr(R[A[1]]);
                }

                error[cellI] = mag(polyU[cellI] - oldPolyU);
            }
        }

        forAll(ibCells, ci)
        {
            Ui[ibCells[ci]] = polyU[ci];
        }

        if(!polyU.empty())
        {
            maxMagPolyU = max(mag(polyU));
        }
        else
        {
            maxMagPolyU = 0;
        }

        error /= maxMagPolyU + SMALL;

        if(!polyU.empty())
        {
            maxError = max(error);
        }
        else
        {
            maxError = 0;
        }

    } while(maxError > 0.01 && counter < 5);
}
//---------------------------------------------------------------------------//
void leastSquaresInt::adjustPhi
(
    leastSquaresIntInfo& intpInfo,
    surfaceScalarField& phi
)
{
    //scalarField& phi = phi.internalField();

    scalar fluxIn = 0;
    scalar fluxOut = 0;
    //scalar fixedFlux = 0;

    const labelList& ibFaces = intpInfo.getIbFaces();
    const boolList& ibFaceFlips = intpInfo.getIbFacesFlips();

    forAll (ibFaces, faceI)
    {
        const label curFace = ibFaces[faceI];
        const bool curFlip = ibFaceFlips[faceI];

        const scalar curFlux = phi[curFace];

        if (!curFlip)
        {
            // Face points out of the live cell
            if (curFlux >= 0)
            {
                // Flux out of the live cell
                fluxOut += curFlux;
            }
            else
            {
                // Flux into the live cell
                fluxIn -= curFlux;
            }
        }
        else
        {
            // Face points into the live cell: flip it
            if (curFlux >= 0)
            {
                // Flux into of the live cell
                fluxIn += curFlux;
            }
            else
            {
                // Flux out the live cell
                fluxOut -= curFlux;
            }
        }
    }

    scalar imbalance = fluxIn - fluxOut;

    Info<< " fluxIn = " << fluxIn << " fluxOut = " << fluxOut
        << " imbalance = " << imbalance
        << endl;

    scalar massCorr = 1.0;

    if (mag(imbalance) > SMALL)
    {
        // Scaling required: scale to match the smaller of two
        // fluxes
        if (fluxIn > fluxOut)
        {
            // Scale down incoming flux
            // Note change of sign: imbalance is negative
            massCorr = 1 - imbalance/(fluxIn + SMALL);

            Info<< "Scaling down incoming flux with factor = "
                << massCorr << endl;

            scalar newFluxIn = 0;

            // Visit all incoming flux faces and re-scale the flux
            forAll (ibFaces, faceI)
            {
                const label curFace = ibFaces[faceI];
                const bool curFlip = ibFaceFlips[faceI];

                if (mesh_.isInternalFace(curFace))
                {
                    // Take reference to current flux
                    scalar& curFlux = phi[curFace];

                    if (!curFlip)
                    {
                        // Face points out of the live cell
                        if (curFlux < 0)
                        {
                            // Flux out of the live cell
                            curFlux *= massCorr;
                            newFluxIn -= curFlux;
                        }
                    }
                    else
                    {
                        // Face points into the live cell: flip it
                        if (curFlux >= 0)
                        {
                            // Flux out the live cell
                            curFlux *= massCorr;
                            newFluxIn += curFlux;
                        }
                    }
                }
            }
        }
        else
        {
            // Scale down outgoing flux
            massCorr = 1 + imbalance/(fluxOut + SMALL);

            Info<< "Scaling down outgoing flux with factor = "
                << massCorr << endl;

            scalar newFluxOut = 0;

            // Visit all outgoing flux faces and re-scale the flux
            forAll (ibFaces, faceI)
            {
                const label curFace = ibFaces[faceI];
                const bool curFlip = ibFaceFlips[faceI];

                if (mesh_.isInternalFace(curFace))
                {
                    // Take reference to current flux
                    scalar& curFlux = phi[curFace];

                    if (!curFlip)
                    {
                        // Face points out of the live cell
                        if (curFlux >= 0)
                        {
                            // Flux out of the live cell
                            curFlux *= massCorr;
                            newFluxOut += curFlux;
                        }
                    }
                    else
                    {
                        // Face points into the live cell: flip it
                        if (curFlux < 0)
                        {
                            // Flux out the live cell
                            curFlux *= massCorr;
                            newFluxOut -= curFlux;
                        }
                    }
                }
            }
        }
    }
}
//---------------------------------------------------------------------------//
