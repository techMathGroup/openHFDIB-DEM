// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)

);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

const volVectorField gradP(fvc::grad(p));
const volScalarField& UEqnA = UEqn.A();
const volScalarField maskedUEqnA(surface*UEqnA);
const volVectorField maskedPressTerm(surface*gradP);

// Forcing update definition
auto updateForcing = [&]()
{
    f = Ui*maskedUEqnA;
    f += maskedPressTerm;
    f -= surface*UEqn.H();
    f.correctBoundaryConditions();
};

// Direct forcing term computation
updateForcing();

if (pimple.momentumPredictor())
{
    for (label i = 0; i < fDirectForcingIter; i++)
    {
        if (max(mag(surface*(Ui - U))).value() < tolUEqn)
        {
            Info << "HFDIBDEM: U converged to Ui within max tolerance " << tolUEqn << endl;
            break;
        }

        Info << "Looping over UEqn, iteration: " << i << endl;
        solve(UEqn == -gradP + f);

        fvOptions.correct(U);

        // Direct forcing term computation
        updateForcing();

        // Use the updated movement to change Ui
        HFDIBDEM.interpolateIB(U,Ui,lambda);
        Ui.correctBoundaryConditions();

        if (i == fDirectForcingIter - 1)
        {
            Info << "Loop over UEqn did not converge to " << tolUEqn << " within " << fDirectForcingIter << " iterations" << endl;
        }
    }
}

// Asign the current f as RHS for momentum equation
UEqn -= f;
