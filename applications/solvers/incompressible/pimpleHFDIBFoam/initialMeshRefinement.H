Info << "Running initial mesh refinement for maxRefinementLevel: " << maxRefinementLevel << endl;
for (label level=0;level<maxRefinementLevel+1;level++)
{
    
    mesh.update();

    Info << "mesh.changing(): " << mesh.changing() << endl;
    // Info << refineF << endl;

    if (mesh.changing())
    {
        Info << "Mesh changed, updating immersed bodies" << endl;
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
    }
}
