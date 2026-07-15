#include "createRDeltaT.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
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


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


// Q (MI, 20250611): do we want to even keep this here? for fluid, we
//                   should include this through fvOptions and for
//                   solids, we might want to load this when initializing
//                   the DEM part of the code
Info<< "\nReading g" << endl;
uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);



singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)      //OF.com: momentumTransportModel -> turbulenceModel
);



#include "createMRF.H"
#include "createFvOptions.H"

// auxiliary variables from HFDIBDEMDict
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

// auxiliary variables from transportProperties
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const dimensionedScalar rho
(
    transportProperties.getOrDefault<dimensionedScalar>("rho", 1.0)
);
if (not transportProperties.found("rho"))
{
    Info << " -- no entry for fluid density found in transportProperties, defaulting to 1.0" << endl;
}

// hfdib-dem additional forces (fluid->solid interaction)
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
    f*rho
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
    f*rho
);

