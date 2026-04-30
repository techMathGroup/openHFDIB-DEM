{
    // Temporarily making U relative to mesh motion
    if (mesh.moving())
    {
        U -= fvc::reconstruct(mesh.phi());
    }

    // Updating alpha1
    alpha1.storePrevIter();
    #include "alphaSuSp.H"
    advector.advect(Sp, Su);

    // -- correction to allow for free phase passing through solids
    alpha1 = (1.0-lambda)*alpha1 + lambda*alpha1.prevIter();
    alpha1.correctBoundaryConditions();

    // Making U absolute again after advection step
    if (mesh.moving())
    {
        U += fvc::reconstruct(mesh.phi());
    }

    // -- original rhoPhi update, nicely consistent with isoAdvector
    #include "rhofs.H"
    rhoPhi = advector.getRhoPhi(rho1f, rho2f);

    // -- updated rhoPhi update, consistent with corrected alpha1
    alpha2 = 1.0 - alpha1;
    // rhoPhi = fvc::interpolate(alpha1*rho1 + alpha2*rho2)*phi;           //prevents entry into liquid
    scalar blend = 0.99;                                                //blending increases stability (slightly)
    // rhoPhi = blend * rhoPhi + (1.0 - blend) * fvc::interpolate(alpha1*rho1 + alpha2*rho2)*phi;
    surfaceScalarField lambdaF(fvc::interpolate(lambda));
    rhoPhi =
      (1.0 - lambdaF)*rhoPhi                                            //outside of solids, use isoAdvector rhoPhi
    + lambdaF*(blend * rhoPhi                                           //inside solids, blend isoAdvector rhoPhi with fvc interpolation
    + (1.0 - blend)
    * fvc::interpolate(alpha1*rho1 + alpha2*rho2)*phi);
    // Note (MI): this overly complex blending is UNDER CONSTRUCTION,
    //            - for 2D water entry case, it increases stability
    //            - was NOT TESTED FOR ANYTHING ELSE
    // Note (MI): motivation is to increase flux consistency (I destroy 
    //            it with the alpha1 blending above) while keeping as
    //            much as possible of the isoAdvector rhoPhi

    mixture.correct();
}

scalar domainFraction = 0;
if (porosityEnabled)
{
    const volScalarField& porosity = tporosity.cref();
    rhoPhi *= scalar(1)/fvc::interpolate(porosity);
    domainFraction = alpha1.weightedAverage(mesh.Vsc()*porosity).value();
}
else
{
    domainFraction = alpha1.weightedAverage(mesh.Vsc()).value();
}

Info<< "Phase-1 volume fraction = "
    << domainFraction
    << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
    << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
    << endl;
