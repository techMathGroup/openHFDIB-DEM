volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
  + MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi, Uf))
);

MRF.makeRelative(phiHbyA);

if (p.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p);
    // enforce update of U and so on
//     HFDIBDEM.updateFSCoupling(lambda,f);
    HFDIBDEM.interpolateIB(U,Ui,lambda);
    f = surface*(UA*Ui - UEqn.H() + fvc::grad(p));
    Ui.correctBoundaryConditions();
    f.correctBoundaryConditions();
    fvc::makeAbsolute(phiHbyA, U);
}

//~ tmp<volScalarField> rAU(rAU);

//~ if (pimple.consistent())
//~ {
    //~ rAU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
    //~ phiHbyA +=
        //~ fvc::interpolate(rAU - rAU)*fvc::snGrad(p)*mesh.magSf();
    //~ HbyA -= (rAU - rAU)*fvc::grad(p);
//~ }

// Note: don't I actually need this later?
surfaceScalarField phiFbyA
(
    fvc::interpolate(rAU)*fvc::flux(f)
);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA+phiFbyA, rAU, MRF);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    // Pressure corrector

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAU, p) ==  fvc::div(phiHbyA)
                                 + fvc::div(phiFbyA)
    );
    // Note: replaced rAU by rAU for consistent PIMPLE
    // Note: is this correct? shouldn't I have some sort of flux form
    //       on the right hand side?

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux() + phiFbyA;
    }
    
    p.correctBoundaryConditions();
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();
// Explicitly relax the forcing term for momentum corrector
f.relax();

U = HbyA - rAU*fvc::grad(p) + rAU*f;//with relaxation of f
//~ U = HbyA - rAU*fvc::grad(p) + fvc::reconstruct(phiFbyA);//consistent with pEqn

U.correctBoundaryConditions();

fvOptions.correct(U);

// Correct Uf if the mesh is moving
fvc::correctUf(Uf, U, phi);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);

//Update imposed field
UA = 1.0/rAU;
UA.correctBoundaryConditions();
// HFDIBDEM.updateFSCoupling(lambda,f);
HFDIBDEM.interpolateIB(U,Ui,lambda);
f = surface*(UA*Ui - UEqn.H() + fvc::grad(p) );
// Note: This part of the code does not look so clean either. If I do
//       not use rAU in consistent, then how did rAU changed?

// post update correction of boundary conditions
f.correctBoundaryConditions();

