Info << "Running initial mesh refinement for maxRefinementLevel: " << maxRefinementLevel << endl;
for (label level=0;level<maxRefinementLevel+1;level++)
{
    
    mesh.update();

    if (mesh.changing())
    {
        Info << "Recreating immersed bodies after mesh update" << endl;
        lambda *= 0.0;
        HFDIBDEM.recreateBodies(lambda,refineF);
    }
}
