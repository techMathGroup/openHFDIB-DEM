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

volVectorField neg_gradP_rgh(
    fvc::reconstruct
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
          - fvc::snGrad(p_rgh)
        ) * mesh.magSf()
    )
);

f = surface*(UEqn.A()*Ui - UEqn.H() - neg_gradP_rgh);
//~ f = surface*(UEqn.A()*Ui - UEqn.H() + fvc::grad(p));
f.correctBoundaryConditions();  

//~ if (pimple.momentumPredictor())
//~ {
    //~ solve
    //~ (
        //~ UEqn
     //~ ==
        //~ fvc::reconstruct
        //~ (
            //~ (
                //~ mixture.surfaceTensionForce()
              //~ - ghf*fvc::snGrad(rho)
              //~ - fvc::snGrad(p_rgh)
            //~ ) * mesh.magSf()
        //~ )
    //~ );

    //~ fvOptions.correct(U);
//~ }

for (int i = 0; i < fDirectForcingIter; i++)
{
    Info << "Looping over UEqn, iteration: " << i << endl;
    solve(UEqn == neg_gradP_rgh + f);

    fvOptions.correct(U);
    
    HFDIBDEM.interpolateIB(U,Ui,lambda);
    
    if (max(mag(surface*(Ui - U))).value() < tolUEqn)
    {
        Info << "HFDIBDEM: U converged to Ui within max tolerance " << tolUEqn << endl;
        break;
    }
    
    Ui.correctBoundaryConditions();
    f += 0.8*surface*rho*(Ui - U)/runTime.deltaT(); 
    f.correctBoundaryConditions(); 
    
    if (i == fDirectForcingIter - 1)
    {
        Info << "Loop over UEqn did not converge within " << fDirectForcingIter << " iterations" << endl;
    }
}

UEqn -= lambda*f;
