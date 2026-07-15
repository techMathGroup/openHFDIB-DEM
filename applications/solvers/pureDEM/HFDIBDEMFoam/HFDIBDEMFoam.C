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
    // OS time efficiency testing

    // scalar preUpdateTime_(0.0);
    // scalar postUpdateTime_(0.0);
    // scalar addRemoveTime_(0.0);
    // scalar updateDEMTime_(0.0);
    // scalar writeBodiesInfoTime_(0.0);
    // scalar meshUpdateTime_(0.0);
    // scalar meshChangingTime_(0.0);
    scalar createBodiesTime_(0.0);

    // OS time efficiency testing
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

        clockTime createBodiesTime; // OS time efficiency testing
        HFDIBDEM.createBodies(lambda,refineF);
        createBodiesTime_ += createBodiesTime.timeIncrement(); // OS time efficiency testing

        // clockTime preUpdateBodiesTime; // OS time efficiency testing
        HFDIBDEM.preUpdateBodies(lambda);
        // preUpdateTime_ += preUpdateBodiesTime.timeIncrement(); // OS time efficiency testing

        // clockTime meshUpdateTime; // OS time efficiency testing
        mesh.update();
        // meshUpdateTime_ += meshUpdateTime.timeIncrement(); // OS time efficiency testing

        // clockTime meshChangingTime; // OS time efficiency testing
        if (mesh.changing())
        {
            lambda *= 0.;
            HFDIBDEM.recreateBodies(lambda,refineF);
        }
        // meshChangingTime_ += meshChangingTime.timeIncrement(); // OS time efficiency testing

        Info << "updating HFDIBDEM" << endl;

        // clockTime postUpdateBodiesTime;
        volVectorField gradLambda(fvc::grad(lambda));
        HFDIBDEM.postUpdateBodies(lambda,gradLambda,f,f);               //MI: here, we should clean up interfaces
        // postUpdateTime_ += postUpdateBodiesTime.timeIncrement();

        // clockTime addRemoveTime;
        HFDIBDEM.addRemoveBodies(lambda,U,refineF);
        // addRemoveTime_ += addRemoveTime.timeIncrement();

        // clockTime updateDEMTime;
        HFDIBDEM.updateDEM(lambda,refineF);
        // updateDEMTime_ += updateDEMTime.timeIncrement();

        Info << "updated HFDIBDEM" << endl;

        runTime.write();

        // clockTime writeBodiesInfoTime;
        if(runTime.outputTime())
        {
            HFDIBDEM.writeBodiesInfo();
        }
        // writeBodiesInfoTime_ += writeBodiesInfoTime.timeIncrement();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        Info << "createBodiesTime    = " << createBodiesTime_    << " s " << endl;

    // Info<< "preUpdateTime       = " << preUpdateTime_       << " s \n"
    //     << "createBodiesTime    = " << createBodiesTime_    << " s \n"
    //     << "meshUpdateTime      = " << meshUpdateTime_      << " s \n"
    //     << "meshChangingTime    = " << meshChangingTime_       << " s \n"
    //     << "postUpdateTime      = " << postUpdateTime_      << " s \n"
    //     << "addRemoveTime       = " << addRemoveTime_       << " s \n"
    //     << "updateDEMTime       = " << updateDEMTime_       << " s \n"
    //     << "   wallContactTime_         = " << HFDIBDEM.wallContactTime_        << " s \n"
    //     << "   wallContactParallelTime_ = " << HFDIBDEM.wallContactParallelTime_<< " s \n"
    //     << "   wallContactSCTime_       = " << HFDIBDEM.wallContactReduceTime_  << " s \n"
    //     << "   prtContactTime_          = " << HFDIBDEM.prtContactTime_         << " s \n"
    //     << "   prtContactParallelTime_  = " << HFDIBDEM.prtContactParallelTime_ << " s \n"
    //     << "   prtContactSCTime_        = " << HFDIBDEM.prtContactReduceTime_   << " s \n"
    //     << "   demItegrationTime_       = " << HFDIBDEM.demItegrationTime_      << " s \n"
    //     << "writeBodiesInfoTime = " << writeBodiesInfoTime_ << " s \n" << endl;

    }

    // Info<< "preUpdateTime       = " << preUpdateTime_       << " s \n"
    //     << "createBodiesTime    = " << createBodiesTime_    << " s \n"
    //     << "meshUpdateTime      = " << meshUpdateTime_      << " s \n"
    //     << "meshChangingTime    = " << meshChangingTime_       << " s \n"
    //     << "postUpdateTime      = " << postUpdateTime_      << " s \n"
    //     << "addRemoveTime       = " << addRemoveTime_       << " s \n"
    //     << "updateDEMTime       = " << updateDEMTime_       << " s \n"
    //     << "   wallContactTime_         = " << HFDIBDEM.wallContactTime_        << " s \n"
    //     << "   wallContactParallelTime_ = " << HFDIBDEM.wallContactParallelTime_<< " s \n"
    //     << "   wallContactSCTime_       = " << HFDIBDEM.wallContactReduceTime_  << " s \n"
    //     << "   prtContactTime_          = " << HFDIBDEM.prtContactTime_         << " s \n"
    //     << "   prtContactParallelTime_  = " << HFDIBDEM.prtContactParallelTime_ << " s \n"
    //     << "   prtContactSCTime_        = " << HFDIBDEM.prtContactReduceTime_   << " s \n"
    //     << "   demItegrationTime_       = " << HFDIBDEM.demItegrationTime_      << " s \n"
    //     << "writeBodiesInfoTime = " << writeBodiesInfoTime_ << " s \n" << endl;
    Info<< "End\n" << endl;

    return 0;
};


// ************************************************************************* //
