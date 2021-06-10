for (label level=0;level<maxRefinementLevel+1;level++)
{
    // read the objects (stl files) and initialize IB
    openHFDIB  HFDIBT(mesh);
    
    HFDIBT.initialize();
    
    HFDIBT.update(lambda,f);
    // Note: from this, I have the lambda field necessary for mesh
    //       refinement.
    // Note: I am trying to use different HFDIB object to update lambda
    //       values for mesh refinement
    
    // update mesh
    Info << "Updating the mesh" << endl;
    mesh.update();
    
    lambda *= scalar(0);
    f      *= scalar(0);
    // Note: I cannot destroy the lambda values due to the mesh
    //       initialization
}
