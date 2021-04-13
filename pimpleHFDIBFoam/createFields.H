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

Info<< "Creating field Ui\n" << endl;
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

Info<< "Reading field lambda\n" << endl;
volScalarField lambda
(
    IOobject
    (
        "lambda",
        runTime.timeName(),
        mesh,
        //~ IOobject::MUST_READ,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimless,0) // initialising to 0
);
// Note: OK, this is not a very nice coding, I read the field and do not
//       use the values at all... a more suitable definition of lambda
//       would be appropriate (MI)

Info<< "Reading field f\n" << endl;
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

Info<< "Creating field refineF\n" << endl;
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
    dimensionedScalar("zero",dimless,0) // initialising to 0
);


#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

#include "readGravitationalAcceleration.H"

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

dimensionedScalar rho = dimensionedScalar(
  transportProperties.lookup("rho")
);

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);


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

scalar thrSurf(readScalar(HFDIBDEMDict.lookup("surfaceThreshold")));
