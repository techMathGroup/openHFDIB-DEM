Info<< "\nReading g" << endl;
//~ const meshObjects::gravity& g = meshObjects::gravity::New(runTime);

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
