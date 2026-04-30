#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// hfdib-dem-related field creation
InfoH<< "Reading field Ui\n" << endl;
volVectorField Ui
(
    IOobject
    (
        "Ui",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U
);

InfoH<< "Reading field lambda\n" << endl;
volScalarField lambda
(
    IOobject
    (
        "lambda",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimless,0)
);

InfoH<< "Reading field f\n" << endl;
volVectorField f
(
    IOobject
    (
        "f",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



InfoH<< "Reading field refineF\n" << endl;
volScalarField refineF
(
    IOobject
    (
        "refineF",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimless,0)
);

InfoH<< "Reading field surface\n" << endl;
volScalarField surface
(
    IOobject
    (
        "surface",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimless,0)
);

#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

// hfdib-dem-related field creation
// Note: this is kind of important to keep rhoF_*F in updateCoupling
//       in immersedBody
InfoH<< "Reading field fDragPress\n" << endl;
volVectorField fDragPress
(
    IOobject
    (
        "fDragPress",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    f
);

InfoH<< "Reading field fDragVisc\n" << endl;
volVectorField fDragVisc
(
    IOobject
    (
        "fDragVisc",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    f
);


typedef incompressibleInterPhaseTransportModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    > transportModelType;

autoPtr<transportModelType> turbulence
(
    new transportModelType
    (
        rho, U, phi, rhoPhi, mixture
    )
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvOptions.H"
#include "createPorosity.H"

isoAdvection advector(alpha1, phi, U);

// hfdib-dem-related dictionary loading
IOdictionary HFDIBDEMDict
(
    IOobject
    (
        "HFDIBDEMDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

scalar thrSurf(SMALL);                                                  //surface threshold for identification 
if(HFDIBDEMDict.found("surfaceThreshold"))
{
    thrSurf = scalar(readScalar(HFDIBDEMDict.lookup("surfaceThreshold")));
}
Info << " -- forcing term will be active in cells with lambda > " << thrSurf << " times" << endl;
label fDragSmoothingIter(5);
if(HFDIBDEMDict.found("nDragForceSmoothingIterations"))
{
    fDragSmoothingIter = label(readLabel(HFDIBDEMDict.lookup("nDragForceSmoothingIterations")));
}
Info << " -- will smooth the particle-fluid interaction forces " << fDragSmoothingIter << " times" << endl;
label fDirectForcingIter(3);
if(HFDIBDEMDict.found("nDirectForcingIterations"))
{
    fDirectForcingIter = label(readLabel(HFDIBDEMDict.lookup("nDirectForcingIterations")));
}
scalar tolUEqn(SMALL);                                                  //surface threshold for identification 
if(HFDIBDEMDict.found("immersedUTolerance"))
{
    tolUEqn = scalar(readScalar(HFDIBDEMDict.lookup("immersedUTolerance")));
    Info << " -- found tolerance for immersed U: " << tolUEqn << endl;
    Info << " -- UEqn will loop until the tolerance is met or for: " << fDirectForcingIter << " iterations" << endl;
}
else
{
    Info << " -- no residual control for immmersed U tolerance set, UEqn will loop for: " << fDirectForcingIter << "iterations" << endl;
}
