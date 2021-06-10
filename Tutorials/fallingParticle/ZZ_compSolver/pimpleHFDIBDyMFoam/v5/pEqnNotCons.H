volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA
//~ (
    //~ "phiHbyA",
    //~ (fvc::interpolate(HbyA) & mesh.Sf())
  //~ + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
//~ );
(
    "phiHbyA",
    fvc::flux(HbyA)
  //~ + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
  + fvc::interpolate(rAU)*fvc::ddtCorr(U, Uf)
);

MRF.makeRelative(phiHbyA);

if (p.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p);
    // enforce update of U and so on (it seems that I need it)
    HFDIB.interpolateIB(U,Ui,lambda);
    f = surface*(UA*Ui - UEqn.H() + fvc::grad(p));
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

//~ if (pimple.nCorrPISO() <= 1)
//~ {
    //~ tUEqn.clear();
//~ }
// Note: don't I actually need this later?

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAU, MRF);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    // Pressure corrector

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAU, p) ==  fvc::div(phiHbyA)
                                 + fvc::div(lambda*rAU*f)
    );
    // Note: replaced rAU by rAU for consistent PIMPLE
    // Note: is this correct? shouldn't I have some sort of flux form
    //       on the right hand side?

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

//~ Info << "After cont Eq" << endl;

// Explicitly relax pressure for momentum corrector
p.relax();

//~ Info << "p relaxed" << endl;

U = HbyA - rAU*fvc::grad(p) + rAU*f;

//~ Info << "U updated " << endl;

U.correctBoundaryConditions();

//~ Info << "U BC corrected " << endl;

fvOptions.correct(U);

{
    Uf = fvc::interpolate(U);
    surfaceVectorField n(mesh.Sf()/mesh.magSf());
    Uf += n*(phi/mesh.magSf() - (n & Uf));
}

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);

//~ Info << "Before HFDIB update " << endl;

//Update imposed field
UA = 1.0/rAU;

//~ Info << "Past the dangerous statement" << endl;
UA.correctBoundaryConditions();
HFDIB.interpolateIB(U,Ui,lambda);
f = surface*(UA*Ui - UEqn.H() +fvc::grad(p) );
// Note: This part of the code does not look so clean either. If I do
//       not use rAU in consistent, then how did rAU changed?

