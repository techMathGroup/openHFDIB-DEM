for (label level=0;level<maxRefinementLevel+1;level++)
{
    if (isA<dynamicRefineFvMesh>(mesh))
    {
        advector.surf().reconstruct();
    }

    mesh.update();

    if (mesh.changing())
    {
        gh = (g & mesh.C()) - ghRef;
        ghf = (g & mesh.Cf()) - ghRef;

        if (isA<dynamicRefineFvMesh>(mesh))
        {
            advector.surf().mapAlphaField();
            alpha2 = 1.0 - alpha1;
            alpha2.correctBoundaryConditions();
            rho == alpha1*rho1 + alpha2*rho2;
            rho.correctBoundaryConditions();
            rho.oldTime() = rho;
            alpha2.oldTime() = alpha2;
        }

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

        lambda *= 0.0;
        HFDIBDEM.recreateBodies(lambda,refineF);
    }
}
