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
    HFDIBDEMFoam

Description
    Transient solver for pure DEM simulations.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "openHFDIBDEM.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "fvcSmooth.H"

#include "triSurfaceMesh.H"
#include "clockTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "createFvOptions.H"
    #include "createMRF.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    #include "readDynMeshDict.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nInitializing HFDIBDEM\n" << endl;
    openHFDIBDEM  HFDIBDEM(mesh);
    HFDIBDEM.initialize(lambda,U,refineF,maxRefinementLevel,runTime.timeName());
    HFDIBDEM.setSolverInfo();
    #include "initialMeshRefinement.H"

    Info<< "\nStarting time loop\n" << endl;

    if(HFDIBDEM.getRecordFirstTime())
    {
        HFDIBDEM.setRecordFirstTime(false);
        HFDIBDEM.writeFirtsTimeBodiesInfo();
    }

    while (runTime.run())
    {

        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        HFDIBDEM.createBodies(lambda,refineF);
        
        HFDIBDEM.preUpdateBodies(lambda,f);

        mesh.update();

        if (mesh.changing())
        {
            lambda *= 0;
            HFDIBDEM.recreateBodies(lambda,refineF);
        }

        Info << "updating HFDIBDEM" << endl;

        HFDIBDEM.postUpdateBodies(lambda,f);

        HFDIBDEM.addRemoveBodies(lambda,U,refineF);

        HFDIBDEM.updateDEM(lambda,refineF);

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
