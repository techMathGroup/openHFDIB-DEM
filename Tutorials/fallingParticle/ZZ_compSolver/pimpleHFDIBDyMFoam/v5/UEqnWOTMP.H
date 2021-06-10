// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
  + g
);

UEqn.relax();

fvOptions.constrain(UEqn);

//if (pimple.momentumPredictor())
//{
    volScalarField UA(UEqn.A());
    UA.correctBoundaryConditions();
        
    // Note: is UA really r? it seems just like UA
    
    HFDIB.interpolateIB(U,Ui,lambda);
    
    f = surface*(UA*Ui - UEqn.H() + fvc::grad(p));
        
    solve(UEqn == -fvc::grad(p) + f);
    
    fvOptions.correct(U);
//}
