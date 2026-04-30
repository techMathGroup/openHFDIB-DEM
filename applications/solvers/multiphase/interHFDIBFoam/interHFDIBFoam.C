/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    interFoam

Group
    grpMultiphaseSolvers

Description
    Solver for two incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"

// hfdib-dem inclusions
#include "triSurfaceMesh.H"
#include "openHFDIBDEM.H"
#include "clockTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using VOF phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"
    

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }
    
    // hfdib-dem inclusions
    #include "readDynMeshDict.H"
    
    // hfdib-dem code modification
    Info << "\nInitializing HFDIBDEM\n" << endl;
    openHFDIBDEM  HFDIBDEM(mesh);
    HFDIBDEM.initialize(lambda,U,refineF,maxRefinementLevel,runTime.timeName());
    #include "initialMeshRefinement.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

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
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        // hfdib-dem code modification
        HFDIBDEM.createBodies(lambda,refineF);
        //~ HFDIBDEM.updateBodiesRhoF(rho);
        HFDIBDEM.updateBodiesRhoF(alpha1,lambda,rho1.value(),rho2.value());
        HFDIBDEM.preUpdateBodies(lambda);

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }
            
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }
        
        // hfdib-dem code modification
        // --- compute viscous forces and update coupling
        volVectorField gradLambda(fvc::grad(lambda));
        //~ scalar omega1(0.5);                                             //formulation weight
        //~ scalar omega1(0.5);                                             //formulation weight
        //~ scalar omega1(0.0);                                             //formulation weight
        //~ scalar omega1(1.0);                                             //formulation weight
        
        //~ fDragPress = fvc::grad(p);
        //~ fDragPress = -fvc::grad(lambda)*p;
        //~ fDragPress = -0.0*gradLambda*p;
        //~ fDragPress = -fvc::ddt(U);
        //~ fDragVisc  = -fvc::div(turbulence->devReff());
        //~ fDragVisc  = -gradLambda & turbulence->devReff();
        //~ fDragVisc  = f;
        //~ fDragPress = -0.5*omega1*(fvc::ddt(U) - f);
        //~ fDragVisc  = fDragPress;
        //~ fDragPress+= (1.0-omega1)*fvc::grad(p);
        //~ fDragVisc += (1.0-omega1)*fvc::div(turbulence->devRhoReff());      //this sign might actually be correct
        //~ fDragPress*= 0.0;
        //~ fDragVisc *= 0.0;
        
        //~ fDragPress = 0.5*f/rho;
        //~ fDragVisc  = fDragPress;
        
        fDragPress = -gradLambda*p;
        
        volTensorField gradU = fvc::grad(U);
        volTensorField tau = -mixture.mu()*(gradU + gradU.T());
        fDragVisc = -gradLambda & tau;
        
        //~ fDragPress /= rho;
        //~ fDragVisc  /= rho;
        for (label pass=0; pass<=fDragSmoothingIter; pass++)
        {
            fDragPress = fvc::average(fvc::interpolate(fDragPress));
            fDragVisc  = fvc::average(fvc::interpolate(fDragVisc));
            fDragPress.correctBoundaryConditions();
            fDragVisc.correctBoundaryConditions();
        }
        
        HFDIBDEM.postUpdateBodies(lambda,gradLambda,fDragPress,fDragVisc);
        HFDIBDEM.addRemoveBodies(lambda,U,refineF);
        HFDIBDEM.updateBodiesRhoF(rho);
        HFDIBDEM.updateDEM(lambda,refineF);
        Info << "updated HFDIBDEM" << endl;

        runTime.write();
        
        if(runTime.outputTime())
        {
            HFDIBDEM.writeBodiesInfo();
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
