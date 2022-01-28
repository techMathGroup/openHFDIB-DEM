for (label level=0;level<maxRefinementLevel+1;level++)
{
    // update mesh
    Info << "Updating the mesh: initial refinment" << endl;
    mesh.update();

    if (mesh.changing())
    {                
        lambda *= 0;
        HFDIBDEM.recreateBodies(lambda,refineF);
    }
}
