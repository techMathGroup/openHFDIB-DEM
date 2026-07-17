// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

// Base operator
tmp<fvVectorMatrix> tUEqnBase
(
    fvm::ddt(U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
);

// current convection operator
tmp<fvVectorMatrix> tConv
(
    fvm::div(phi, U)
);

// helper for assembling the complete momentum equation
auto assembleUEqn =
[&]() -> tmp<fvVectorMatrix>
{
    tmp<fvVectorMatrix> tEqn(new fvVectorMatrix(tUEqnBase()));

    fvVectorMatrix& eqn = tEqn.ref();

    eqn += tConv();

    eqn.relax();

    fvOptions.constrain(eqn);

    return tEqn;
};

const volVectorField gradP(fvc::grad(p));
const volVectorField maskedPressTerm(surface*gradP);
// Forcing update definition
auto updateForcing = [&](const fvVectorMatrix& eqn)
{
    const volScalarField& UEqnA = eqn.A();
    const volScalarField maskedUEqnA(surface*UEqnA);

    f = Ui*maskedUEqnA;
    f += maskedPressTerm;
    f -= surface*eqn.H();
    f.correctBoundaryConditions();
};

// initial momentum equation
tmp<fvVectorMatrix> tUEqn = assembleUEqn();
updateForcing(tUEqn());

if (pimple.momentumPredictor())
{
    for (label i = 0; i < fDirectForcingIter; ++i)
    {
        if (max(mag(surface*(Ui - U))).value() < tolUEqn)
        {
            Info<< "HFDIBDEM: U converged to Ui within max tolerance "
                << tolUEqn << endl;
            break;
        }

        Info<< "Looping over UEqn, iteration: " << i << endl;
        fvVectorMatrix& UEqnIter = tUEqn.ref();
        solve(UEqnIter == -gradP + f);

        fvOptions.correct(U);

        // Use the updated movement to change Ui
        HFDIBDEM.interpolateIB(U, Ui, lambda);
        Ui.correctBoundaryConditions();

        // Rebuild only convection
        phi = fvc::flux(U);
        MRF.makeRelative(phi);
        fvc::makeRelative(phi, U);
        tConv = fvm::div(phi, U);

        // Assemble fresh momentum equation
        tUEqn = assembleUEqn();
        updateForcing(tUEqn());

        if (i == fDirectForcingIter - 1)
        {
            Info<< "Loop over UEqn did not converge to "
                << tolUEqn
                << " within "
                << fDirectForcingIter
                << " iterations"
                << endl;
        }
    }
}

// Assemble final momentum equation for the pressure corrector
tConv = fvm::div(phi, U);
tUEqn = assembleUEqn();
fvVectorMatrix& UEqn = tUEqn.ref();
UEqn -= f;