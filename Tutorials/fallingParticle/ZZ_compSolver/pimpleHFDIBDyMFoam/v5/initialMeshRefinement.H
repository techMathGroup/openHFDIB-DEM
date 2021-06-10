{
    // read the objects (stl files) and initialize IB
    openHFDIB  HFDIB2(mesh);
    HFDIB2.initialize();
    
    HFDIB2.update(lambda,f);
    // Note: from this, I have the lambda field necessary for mesh
    //       refinement.
    // Note: I am trying to use different HFDIB object to update lambda
    //       values for mesh refinement
        
    //~ forAll(Ui,cellI)
    //~ {
        //~ Ui[cellI].component(0)=0.0;
        //~ if( lambda[cellI]>0.001 ) surface[cellI] =1.0;
        //~ else surface[cellI] =0.0;
    //~ }
    
    // update mesh
    Info << "Updating the mesh" << endl;
    mesh.update();
    
    //~ mesh.update();
    
    //~ // Calculate absolute flux from the mapped surface velocity
    //~ phi = mesh.Sf() & Uf;
//~ 
    //~ if (mesh.changing() && correctPhi)
    //~ {
        //~ #include "correctPhi.H"
    //~ }
//~ 
    //~ // Make the flux relative to the mesh motion
    //~ fvc::makeRelative(phi, U);
    //~ 
    //~ if (mesh.changing() && checkMeshCourantNo)
    //~ {
        //~ #include "meshCourantNo.H"
    //~ }
    
    lambda *= scalar(0);
    f      *= scalar(0);
    // Note: I cannot destroy the lambda values due to the mesh
    //       initialization
    
    //~ mesh.write();
    
}
