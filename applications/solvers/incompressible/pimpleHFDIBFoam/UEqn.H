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

volVectorField gradP(fvc::grad(p));

// Direct forcing term computation
f = surface*(UEqn.A()*Ui - UEqn.H() + gradP);
// f = lambda*(UEqn.A()*Ui - UEqn.H() + gradP);
f.correctBoundaryConditions(); 

if (pimple.momentumPredictor())
{
    // solve(UEqn == -fvc::grad(p) + f);

    // fvOptions.correct(U);
    for (int i = 0; i < fDirectForcingIter; i++)
    {
        Info << "Looping over UEqn, iteration: " << i << endl;
        solve(UEqn == -fvc::grad(p) + f);

        fvOptions.correct(U);

        if (max(mag(surface*(Ui - U))).value() < tolUEqn)
        {
            Info << "HFDIBDEM: U converged to Ui within max tolerance " << tolUEqn << endl;
            break;
        }


        // f += 0.8*surface*(Ui - U)/runTime.deltaT();                     //MI: in this loop, correct the force incrementaly
        f = surface*(UEqn.A()*Ui - UEqn.H() + gradP);
        // f = lambda*(UEqn.A()*Ui - UEqn.H() + gradP);
        f.correctBoundaryConditions(); 

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
