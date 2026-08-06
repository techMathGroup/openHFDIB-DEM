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
f.storePrevIter();
scalar fCoeff = 0.8;
// f = (1.0-fCoeff)*f + fCoeff*surface*(UEqn.A()*Ui - UEqn.H() - neg_gradP_rgh);
f = surface*(UEqn.A()*Ui - UEqn.H() - neg_gradP_rgh);
f.correctBoundaryConditions(); 
// f.storePrevIter();
// for (label i = 0; i < 3; ++i)
// {
//     // sync the flow with the body movement
//     HFDIBDEM.updateFSCoupling(lambda,f);
//     HFDIBDEM.interpolateIB(U,Ui,lambda);
//     Ui.correctBoundaryConditions();

//     // f = surface*(UEqn.A()*Ui - UEqn.H() - neg_gradP_rgh);
//     f = (1.0-fCoeff)*f + fCoeff*surface*(UEqn.A()*Ui - UEqn.H() - neg_gradP_rgh);
//     f.correctBoundaryConditions(); 
// }

// // sync the flow with the body movement
// HFDIBDEM.updateFSCoupling(lambda,f);
// HFDIBDEM.interpolateIB(U,Ui,lambda);
// Ui.correctBoundaryConditions();

if (pimple.momentumPredictor())
{
    for (label i = 0; i < fDirectForcingIter; ++i)
    {
        if (max(mag(surface*(Ui - U))).value() < tolUEqn)
        {
            Info << "HFDIBDEM: U converged to Ui within max tolerance " << tolUEqn << endl;
            break;
        }

        Info << "Looping over UEqn, iteration: " << i << endl;
        solve(UEqn == neg_gradP_rgh + f);

        fvOptions.correct(U);

        // HFDIBDEM.updateFSCoupling(lambda,f);
        // HFDIBDEM.interpolateIB(U,Ui,lambda);
        // Ui.correctBoundaryConditions();

        f = (1.0-fCoeff)*f + fCoeff*surface*(UEqn.A()*Ui - UEqn.H() - neg_gradP_rgh);
        f.correctBoundaryConditions(); 

        if (i == fDirectForcingIter - 1)
        {
            Info << "Loop over UEqn did not converge within " << fDirectForcingIter << " iterations" << endl;
        }
    }
}

UEqn -= f;// - solidMomentumCorr;
