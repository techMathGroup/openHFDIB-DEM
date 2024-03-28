// assume static mesh
label maxRefinementLevel = 0;

// load dynamicMeshDict
IOdictionary dynamicMeshDict
(
    IOobject
    (
        "dynamicMeshDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// get the dynamicFvMeshType
const word dynamicFvMeshTypeName(dynamicMeshDict.lookup("dynamicFvMesh"));
Info<< "Just to confirm, working with " << dynamicFvMeshTypeName << endl;

const word stMesh("staticFvMesh");
const word rfMesh("dynamicRefineFvMesh");

if (dynamicFvMeshTypeName!=stMesh)
{
    if (dynamicFvMeshTypeName==rfMesh)
    {
        maxRefinementLevel = readLabel(dynamicMeshDict.subDict("dynamicRefineFvMeshCoeffs").lookup("maxRefinement"));
        
        Info << "Will try to refine mesh based on the bodies location" << endl;
    }
}
