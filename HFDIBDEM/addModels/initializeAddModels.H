    // Note: if you want to add a new addModel, do it here
    
    scalar thrSurf(readScalar(HFDIBDEMDict_.lookup("surfaceThreshold")));
    
    forAll (bodyNames_,nameI)
    {
        word bodyName(bodyNames_[nameI]);
        Info << "-- addModelMessage-- " << "Creating addition model for body: " << bodyName << endl;
        
        // initialize the addition models
        dictionary addModelDict(HFDIBDEMDict_.subDict(bodyName).subDict("bodyAddition"));
        word howToAdd(word(addModelDict.lookup("addModel")));
        word bodyGeom;
        autoPtr<geomModel> bodyGeomModel;
        // check if the immersedDict_ contains bodyGeom
        if (HFDIBDEMDict_.subDict(bodyName).found("bodyGeom"))
        {  
            word input = HFDIBDEMDict_.subDict(bodyName).lookup("bodyGeom");
            bodyGeom = input;
            Info << "Found bodyGeom for " << bodyName << ", the body is: " << bodyGeom << endl;
        }
        else
        {
            bodyGeom = "convex";
            Info << "Did not find bodyGeom for " << bodyName << ", using bodyGeom: " << bodyGeom << endl;
        }
        
        if(bodyGeom == "convex")
        {
            word stlPath("constant/triSurface/"+bodyName+".stl");
            bodyGeomModel.set(new convexBody(mesh_,stlPath,thrSurf,geometricD_));
        }
        else if(bodyGeom == "nonConvex")
        {            
            word stlPath("constant/triSurface/"+bodyName+".stl");
            bodyGeomModel.set(new nonConvexBody(mesh_,stlPath,thrSurf,geometricD_));
        }
        else if(bodyGeom == "sphere")
        {
            vector startPosition = HFDIBDEMDict_.subDict(bodyName).subDict("sphere").lookup("startPosition");
            scalar radius = readScalar(HFDIBDEMDict_.subDict(bodyName).subDict("sphere").lookup("radius"));
            
            bodyGeomModel.set(new sphereBody(mesh_,startPosition,radius,thrSurf,geometricD_));
        }
        else
        {
            word stlPath("constant/triSurface/"+bodyName+".stl");
            Info << "bodyGeom: " << bodyGeom << " not supported, using bodyGeom nonConvex" << endl;
            bodyGeom = "nonConvex";
            bodyGeomModel.set(new nonConvexBody(mesh_,stlPath,thrSurf,geometricD_));
        }
        
        Info << "-- addModelMessage-- " << "trying to set addModel: " << howToAdd << endl;
        
        if (howToAdd == "once")
        {
            addModels_.set
            (
                nameI,
                new addModelOnce
                (
                    addModelDict,
                    mesh_,
                    startTime0,
                    bodyGeomModel.ptr()
                )
            );
        }
        else if (howToAdd == "onceScatter")
        {
            addModels_.set
            (
                nameI,
                new addModelOnceScatter
                (
                    addModelDict,
                    mesh_,
                    startTime0,
                    geometricD_,
                    bodyGeomModel.ptr()
                )
            );
        }
        else if (howToAdd == "repeatSamePosition")
        {
            addModels_.set
            (
                nameI,
                new addModelRepeatSamePosition
                (
                    addModelDict,
                    mesh_,
                    bodyGeomModel.ptr()
                )
            );
        }
        else if (howToAdd == "repeatRandomPosition")
        {
            addModels_.set
            (
                nameI,
                new addModelRepeatRandomPosition
                (
                    addModelDict,
                    mesh_,
                    geometricD_,
                    bodyGeomModel.ptr()
                )
            );
        }
        else if (howToAdd == "distribution")
        {
            addModels_.set
            (
                nameI,
                new addModelDistribution
                (
                    addModelDict,
                    mesh_,
                    geometricD_,
                    bodyGeomModel.ptr()
                )
            );
        }
        else
        {
            Info << "-- addModelMessage-- " << "No addition model specified, defaulting to addModelOnce" << endl;
            
            addModels_.set
            (
                nameI,
                new addModelOnce
                (
                    addModelDict,
                    mesh_,
                    startTime0,
                    bodyGeomModel.ptr()
                )
            );
        }
    }
