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
    
    The code is prepared to use HFDIB method and refine the mesh around
    the solids via fvdynamicRefineMesh
    
    Note: the general code structure is the same as in pisoHFDIBFoam
          by Federico Municchi but the HFDIB library itself was slightly
          changed
    Note: although the code should theoretically work with moving meshes
          I only borrow the mesh refinement functionality. the 
          initialMeshRefinement.H is completely useles for the rest of
          dynamic meshes types


Changelog:
* 20190107
    - hoping to correct the problems with hanging pointers
    - the solver still crashes with deltaT -> 0
    - the solver still crashes when the consistent version of PIMPLE is
      used
    - if less than 2 nCorrectors are used, the solver crashes with
      "sizes of addressing and field are different" (WTF?)
    - nCorrectors = 2, nNonOrthogonalCorrectors = 1 => crash with
      "corrupted size vs. prev_size: 0x0000000006d96630"
    - nCorrectors = 2, nNonOrthogonalCorrectors = 2 => runs (CN)
                                                    => crashes (Euler)
      "free(): invalid pointer: 0x0000000006fcb0c0"
    - nCorrectors = 2, nNonOrthogonalCorrectors = 3 => runs but unstable
    - the problem with pointers does not seem to be related to temporary
      fields. neither to temporary matrices
    - manifestation and type of the pointer problem depends on the
      default gradScheme
* 20190110
    - updated the HFDIB class to remove dynamic lists (vectors), which
      solved a lot of problems with mesh refinement
    - forced re-initialization of lambda and f fields on the mesh update
    - !! I need to implement a refinement of the stl surface up to the
      refinementLevel !! (load stl, refine) x refinementLevel
      Note: this will probably force me to load (somehow) the maximal
            mesh refinement level from the mesh object
      Q:    would it be possible to move this whole mesh update and
            stl reloading into the immersed body initialization?
            (I want to keep this class as simple as possible)
* 20190111
    - managed to somehow update the mesh to level 2 (but I needed to
      hardcode the level!?) -> no I don't (see below)
    - CREATE a readDynamicMeshDict.H in which you will get the dynamic
      mesh type and if needed, also the maximum refinement level
      (seems to be the easiest way around)
    - I run the initial mesh-stl ping-pong maxRefinement+1 times to
      enable the final mesh update (split-points)
    - now, with change in the mesh, the solver seems to crash once more
    - created readDynMeshDict.H file in which I explore the dynamicMeshDict
      and if needed, read the maxRefinement and use it later on
    - I still need to solve the bug with solver crashing when I do not
      output each time step at the begining of the simulation
      -> it did not crush (but it may be due to shear luck)
    - I still need to make the solver run in parallel
    - more vigorous testing will be needed (I can test the stl files from
      Martin's (Š) BC thesis and maybe put together an article on
      comparison of this method with standard OpenFOAM solvers (framework)
      -> look into the case handling complexity, runTime and memory
         requirements

* 20190122
    - the solver seems to work OK (fingers crossed) in serial mode
    - I need to make it run in parallel
    - added forced update of U and Ui if pressure needs a reference
      (now solver does not crash when ran on closed domain, but it is
      far from ideal) => !! STILL CRASHES ON CLOSED DOMAIN !!


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "openHFDIB.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H" // OF 4.x

#include "triSurfaceMesh.H"  // for fd read from stl file

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createControls.H"
    #include "createFields.H"
    #include "createUf.H"
    #include "createFvOptions.H"
    #include "createMRF.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    
    #include "readDynMeshDict.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    volScalarField surface("surf",lambda);
    
    bool isFirstTime(true);
    openHFDIB  HFDIB(mesh);

    while (runTime.run())
    {
        
        #include "readControls.H"
        #include "CourantNo.H"

        #include "setDeltaT.H"
        
        runTime++;
        
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        if (isFirstTime)
        {
            Info << "\nInitializing HFDIB\n" << endl;
            
            if (maxRefinementLevel > 0)
            {
                #include "initialMeshRefinementV2.H"
            }
        }
        else
        {
            Info << "Updating mesh" << endl;
            mesh.update();
            Info << "Updated mesh" << endl;
        }

        // Calculate absolute flux from the mapped surface velocity
        Info << "updating phi" << endl;
        phi = mesh.Sf() & Uf;
        
        Info << "updated phi" << endl;

        if (mesh.changing())
        {
            HFDIB.forceReinit(lambda,f);
            // If there is a change in the mesh topology, I need to
            // reinitialize all the immersed bodies
            // Note: what does happen if more than one immersed body is
            //       created?
            // Note: I should forceReinit even if correctPhi=false
        }
        if (mesh.changing() && correctPhi)
        {
            #include "correctPhi.H"
            Info << "Corrected Phi" << endl;
        }

        // Make the flux relative to the mesh motion
        fvc::makeRelative(phi, U);
        
        if (mesh.changing() && checkMeshCourantNo)
        {
            #include "meshCourantNo.H"
            Info << "Computed Co" << endl;
        }
        
        if (isFirstTime)
        {
            HFDIB.initialize(); 
            isFirstTime = false;
        }
        
        Info << "trying to update HFDIB" << endl;
        HFDIB.update(lambda,f);
        surface = lambda;
        Info << "updated HFDIB" << endl;
        
        forAll(Ui,cellI)
        {
            //~ Ui[cellI].component(0)=0.0;
            //~ if( lambda[cellI]>0.0001 ) surface[cellI] =1.0;
            //~ else surface[cellI] =0.0;
            if( lambda[cellI]>0.0001 )
            {
                surface[cellI] =1.0;
            }
            else 
            {
                surface[cellI] =0.0;
                Ui[cellI] *=0;
            }
        }
        // Note: well, this is an ugly code snippet. I initialize 
        //       surface with lambda, but I do not update its boundary
        //       conditions, ever -> this might be giving me the
        //       segFault
        // Note: even changing the position of surface field
        //       initialization will change the problem, I should still
        //       repair this to take into account boundaries

        
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqnWOTMP.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqnNotCons.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }


        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
            
    }

    Info<< "End\n" << endl;

    return 0;
};


// ************************************************************************* //
