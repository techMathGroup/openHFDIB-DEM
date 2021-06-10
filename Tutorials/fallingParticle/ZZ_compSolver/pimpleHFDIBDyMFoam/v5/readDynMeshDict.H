// Note: in this file, I need to read the necessary data from the
//       dynamicMeshDict (if it is present?)
label maxRefinementLevel = 0;
//~ scalar maxRefinementLevel = 0;
//~ bool  staticMesh         = true;
// Note: assume staticFvMesh

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
    //~ staticMesh=false;
    if (dynamicFvMeshTypeName==rfMesh)
    {
        maxRefinementLevel = readLabel(dynamicMeshDict.subDict("dynamicRefineFvMeshCoeffs").lookup("maxRefinement"));
        
        Info << "Will try to refine mesh based on the bodies location" << endl;
    }
}
