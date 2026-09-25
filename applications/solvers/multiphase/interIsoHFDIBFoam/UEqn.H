MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevRhoReff(rho, U)
 ==
    fvOptions(rho, U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

// hfdib-dem modifications
volScalarField UA(UEqn.A());
HFDIBDEM.interpolateIB(U,Ui,lambda);
Ui.correctBoundaryConditions();

// surfaceScalarField lambdaF(fvc::interpolate(lambda));
// surfaceScalarField lambdaF(fvc::interpolate(surface));
volVectorField ghf_gradRho(
    fvc::reconstruct(
        ghf*fvc::snGrad(rho) * mesh.magSf() * (1.0 - surfaceF)
    )
);
ghf_gradRho.correctBoundaryConditions();

// const volVectorField neg_gradP_rgh(
volVectorField neg_gradP_rgh(
    fvc::reconstruct
    (
        (
            mixture.surfaceTensionForce()
        //   - ghf*fvc::snGrad(rho)
          - fvc::snGrad(p_rgh)
        ) * mesh.magSf()
    )
);
neg_gradP_rgh -= ghf_gradRho;
neg_gradP_rgh.correctBoundaryConditions();

// float-stepper test implementation
scalar fCoeff = 0.8;
f.storePrevIter();
// f = surface*(UEqn.A()*Ui - UEqn.H() - neg_gradP_rgh);
f = (1.0-fCoeff)*f + fCoeff*surface*(UEqn.A()*Ui - UEqn.H() - neg_gradP_rgh);
f.correctBoundaryConditions();  

// surfaceScalarField phiSolid(fvc::flux(Ui));
// surfaceScalarField rhoPhiSolid = 
//     fvc::interpolate(rhoSolidExcess)*phiSolid;
// surfaceScalarField rhoPhiSolid = 
//     fvc::interpolate(rhoSolidExcess)*phi;
// f+= fvc::div(rhoPhiSolid, U, "div(phiSolid,U)") + fvc::ddt(rhoS, U);
// f.correctBoundaryConditions();

// volVectorField solidConvCorr
// (
//     IOobject
//     (
//         "solidConvCorr",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     fvc::div(rhoPhiSolid, U, "div(phiSolid,U)")
// );

// volVectorField solidDdtCorr
// (
//     IOobject
//     (
//         "solidDdtCorr",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     // fvc::ddt(rhoSolidExcess, U)
//     (rhoSolidExcess*fvc::ddt(lambda) - lambda*fvc::ddt(rho))*U
// );

// volVectorField solidMomentumCorr
// (
//     IOobject
//     (
//         "solidMomentumCorr",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     0*solidConvCorr + 0*solidDdtCorr
// );

if (pimple.momentumPredictor())
{
    for (label i = 0; i < fDirectForcingIter; ++i)
    {
        if (max(mag(surface*(Ui - U))).value() < tolUEqn)
        {
            Info << "HFDIBDEM: U converged to Ui within max tolerance " << tolUEqn << endl;
            break;
        }
        
        // Info<< "convection:         " << max(mag(fvc::div(rhoPhi,U))).value() << endl;
        // Info<< "solid momentum:     " << max(mag(solidMomentumCorr)).value() << endl;
        // Info<< "solid convection:   " << max(mag(solidConvCorr)).value() << endl;
        // Info<< "solid ddt:          " << max(mag(solidDdtCorr)).value() << endl;

        Info << "Looping over UEqn, iteration: " << i << endl;
        solve(UEqn == neg_gradP_rgh + f);// - solidMomentumCorr);

        fvOptions.correct(U);

        HFDIBDEM.interpolateIB(U,Ui,lambda);
        Ui.correctBoundaryConditions();

        // solidConvCorr = fvc::div(rhoPhiSolid, U, "div(phiSolid,U)");
        // solidDdtCorr = fvc::ddt(rhoS, U);
        // solidDdtCorr = (rhoSolidExcess*fvc::ddt(lambda) - lambda*fvc::ddt(rho))*U;
        // solidMomentumCorr = solidConvCorr + solidDdtCorr;

        // phiSolid = fvc::flux(Ui);
        // rhoPhiSolid =
        //         fvc::interpolate(rhoSolidExcess)*phiSolid;
        // rhoPhiSolid =
        //         fvc::interpolate(rhoSolidExcess)*phi;

        f = (1.0-fCoeff)*f + fCoeff*surface*(UEqn.A()*Ui - UEqn.H() - neg_gradP_rgh);
        // f+= fvc::div(rhoPhiSolid, Ui, "div(phiSolid,U)") + fvc::ddt(rhoS, Ui);
        // f+= fvc::div(rhoPhiSolid, U, "div(phiSolid,U)") + fvc::ddt(rhoS, U);
        f.correctBoundaryConditions(); 

        if (i == fDirectForcingIter - 1)
        {
            Info << "Loop over UEqn did not converge within " << fDirectForcingIter << " iterations" << endl;
        }
    }
}

UEqn -= f;// - solidMomentumCorr;
