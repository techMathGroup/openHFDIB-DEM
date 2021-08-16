// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevSigma(U)
 ==
    fvOptions(U)
  + g
);


UEqn.relax();

fvOptions.constrain(UEqn);

volScalarField UA(UEqn.A());
UA.correctBoundaryConditions();
HFDIBDEM.updateFSCoupling(lambda,f);                                    //updates V_el, omega_ and Axis_
HFDIBDEM.interpolateIB(U,Ui,lambda);
Ui.correctBoundaryConditions();

f = surface*(UEqn.A()*Ui - UEqn.H() + fvc::grad(p));
f.correctBoundaryConditions();    

if (pimple.momentumPredictor())
{    
    solve(UEqn == -fvc::grad(p) + f);
    fvOptions.correct(U);
}

