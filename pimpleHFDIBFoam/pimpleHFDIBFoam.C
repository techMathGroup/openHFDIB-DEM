/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pimpleHFDIBDyMFoam

Description
    Transient solver for incompressible flow.

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

    The code is prepared to use the HFDIB method and refine the mesh around
    the solids via fvdynamicRefineMesh
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "openHFDIBDEM.H"
#include "singlePhaseTransportModel.H"
#include "kinematicMomentumTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

#include "triSurfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "createFvOptions.H"
    #include "createMRF.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    #include "readDynMeshDict.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    volScalarField surface("surf",lambda);
    scalar thrSurf(readScalar(HFDIBDEMDict.lookup("surfaceThreshold")));

    bool isFirstTime(true);
    openHFDIBDEM  HFDIBDEM(mesh);

    while (runTime.run())
    {

        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        word startTime(runTime.timeName());

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (isFirstTime)
        {
            Info << "\nInitializing HFDIBDEM\n" << endl;

            HFDIBDEM.initialize(lambda,U,refineF,maxRefinementLevel,startTime);

            if (maxRefinementLevel > 0 && startTime == "0")
            {
                #include "initialMeshRefinement.H"
            }

            isFirstTime = false;
        }

        HFDIBDEM.preUpdateBodies(lambda,f);

        surface = lambda;
        forAll (Ui,cellI)
        {
            if (lambda[cellI]>thrSurf)
            {
                surface[cellI] =1.0;
            }
            else
            {
                surface[cellI] =0.0;
                Ui[cellI] *=0;
            }
        }

        surface.correctBoundaryConditions();
        Ui.correctBoundaryConditions();

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }

                    lambda *= 0;
                    HFDIBDEM.recreateBodies(lambda,refineF);

                    surface = lambda;
                    forAll (Ui,cellI)
                    {
                        if (lambda[cellI]>thrSurf)
                        {
                            surface[cellI] =1.0;
                        }
                        else
                        {
                            surface[cellI] =0.0;
                            Ui[cellI] *=0;
                        }
                    }

                    surface.correctBoundaryConditions();
                    Ui.correctBoundaryConditions();
                }
            }
            f.storePrevIter();

            #include "UEqnExtF.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqnExtF.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        Info << "updating HFDIBDEM" << endl;
        HFDIBDEM.postUpdateBodies(lambda,f);

        HFDIBDEM.addRemoveBodies(lambda,U,refineF);

        HFDIBDEM.moveBodies(lambda,refineF);

        HFDIBDEM.correctContact(lambda,refineF);
        Info << "updated HFDIBDEM" << endl;


        runTime.write();

        if(runTime.outputTime())
        {
            HFDIBDEM.writeBodiesInfo();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }

    Info<< "End\n" << endl;

    return 0;
};


// ************************************************************************* //
