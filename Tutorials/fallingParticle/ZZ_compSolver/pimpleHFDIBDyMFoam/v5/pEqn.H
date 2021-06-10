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
  + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
);

MRF.makeRelative(phiHbyA);

if (p.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p);
    fvc::makeAbsolute(phiHbyA, U);
}

tmp<volScalarField> rAtU(rAU);

if (pimple.consistent())
{
    rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p);
}

//~ if (pimple.nCorrPISO() <= 1)
//~ {
    //~ tUEqn.clear();
//~ }
// Note: don't I actually need this later?

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAtU(), MRF);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    // Pressure corrector

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p) ==  fvc::div(phiHbyA)
                                 //~ + fvc::div(lambda*rAtU()*f)
        //~ fvm::laplacian(rAtU, p) ==  fvc::div(phiHbyA)
                                 + fvc::div(lambda*rAtU*f)
    );
    // Note: replaced rAU by rAtU for consistent PIMPLE
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

//~ U = HbyA - rAtU()*fvc::grad(p) + rAtU()*f;
U = HbyA - rAtU*fvc::grad(p) + rAtU*f;
//~ U = HbyA - rAtU*fvc::grad(p);
U.correctBoundaryConditions();
// Note: I believe that there is a bug in the original pimpleDyMFoam
//       (rAtU -> rAtU(), check higher oF versions)

Info << "U updated " << endl;

fvOptions.correct(U);

{
    Uf = fvc::interpolate(U);
    surfaceVectorField n(mesh.Sf()/mesh.magSf());
    Uf += n*(phi/mesh.magSf() - (n & Uf));
}

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);

Info << "Before HFDIB update " << endl;

//Update imposed field
//~ UA = 1.0/rAtU();
UA = 1.0/rAtU;

Info << "Past the dangerous statement" << endl;
UA.correctBoundaryConditions();
HFDIB.interpolateIB(U,Ui,lambda);
f = surface*(UA*Ui - UEqn.H() +fvc::grad(p) );
// Note: This part of the code does not look so clean either. If I do
//       not use rAtU() in consistent, then how did rAU changed?

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

