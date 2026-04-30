/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    pimpleHFDIBFoam

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids,
    with optional mesh motion and mesh topology changes.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"

#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

#include "triSurfaceMesh.H"
#include "openHFDIBDEM.H"
#include "clockTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    #include "readDynMeshDict.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nInitializing HFDIBDEM\n" << endl;
    openHFDIBDEM  HFDIBDEM(mesh);
    HFDIBDEM.initialize(lambda,U,refineF,maxRefinementLevel,runTime.timeName());
    #include "initialMeshRefinement.H"

    if(HFDIBDEM.getRecordFirstTime())
    {
        HFDIBDEM.setRecordFirstTime(false);
        HFDIBDEM.writeBodiesInfo();
    }

    Info<< "\nStarting time loop\n" << endl;

    scalar CFDTime_(0.0);
    scalar DEMTime_(0.0);
    scalar suplTime_(0.0);

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

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        clockTime createBodiesTime; // OS time efficiency testing
        HFDIBDEM.createBodies(lambda,refineF);
        HFDIBDEM.updateBodiesRhoF(rho.value());
        suplTime_ += createBodiesTime.timeIncrement(); // OS time efficiency testing

        clockTime preUpdateBodiesTime; // OS time efficiency testing
        HFDIBDEM.preUpdateBodies(lambda);
        suplTime_ += preUpdateBodiesTime.timeIncrement(); // OS time efficiency testing
        
        // --- pre-compute gradient of lambda field (force updates)
        volVectorField gradLambda(fvc::grad(lambda));
        gradLambda.correctBoundaryConditions();        
        
        // --- construct surface field where the momentum source should
        //     be switched on
        forAll(surface, sI)
        {
            if (lambda[sI] > thrSurf)
                surface[sI] = 1;
            else
                surface[sI] = 0;
        }
        surface.correctBoundaryConditions();

        clockTime pimpleRunClockTime; // OS time efficiency testing
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
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

                    lambda *= 0.0;

                    HFDIBDEM.recreateBodies(lambda,refineF);
                    
                    volVectorField gradLambda(fvc::grad(lambda));                    
                    forAll(surface, sI)
                    {
                        if (lambda[sI] > thrSurf)
                            surface[sI] = 1;
                        else
                            surface[sI] = 0;
                    }
                    gradLambda.correctBoundaryConditions();
                    surface.correctBoundaryConditions();
                }

            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }
        CFDTime_ += pimpleRunClockTime.timeIncrement();
        Info << "updating HFDIBDEM" << endl;
        clockTime postUpdateBodiesTime;
        
        fDragVisc = (f - fvc::grad(p))*rho;
        fDragPress= -gradLambda*p*rho;
        
        fDragPress.correctBoundaryConditions();
        fDragVisc.correctBoundaryConditions();
        
        for (label pass=0; pass<=fDragSmoothingIter; pass++)
        {
            fDragPress = fvc::average(fvc::interpolate(fDragPress));
            fDragVisc  = fvc::average(fvc::interpolate(fDragVisc));
            fDragPress.correctBoundaryConditions();
            fDragVisc.correctBoundaryConditions();
        }
        
        HFDIBDEM.postUpdateBodies(lambda,gradLambda,fDragPress,fDragVisc);
        suplTime_ += postUpdateBodiesTime.timeIncrement();


        clockTime addRemoveTime;
        HFDIBDEM.addRemoveBodies(lambda,U,refineF);
        HFDIBDEM.updateBodiesRhoF(rho.value());
        suplTime_ += addRemoveTime.timeIncrement();

        clockTime updateDEMTime;
        HFDIBDEM.updateDEM(lambda,refineF);
        DEMTime_ += updateDEMTime.timeIncrement();
        Info << "updated HFDIBDEM" << endl;


        runTime.write();

        clockTime writeBodiesInfoTime;
        if(runTime.outputTime())
        {
            HFDIBDEM.writeBodiesInfo();
        }
        suplTime_ += writeBodiesInfoTime.timeIncrement();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        Info<< " CFDTime_                 = " << CFDTime_             << " s \n" <<
               " Solver suplementary time = " << suplTime_            << " s \n" << 
               " DEMTime_                 = " << DEMTime_             << " s \n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
