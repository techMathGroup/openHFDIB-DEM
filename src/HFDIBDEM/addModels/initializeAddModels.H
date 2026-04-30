    // Note: if you want to add a new addModel, do it here

    scalar thrSurf(readScalar(HFDIBDEMDict_.lookup("surfaceThreshold")));

    forAll (bodyNames_,nameI)
    {
        word bodyName(bodyNames_[nameI]);
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "Creating addition model for body: " << bodyName << endl;

        // initialize the addition models
        dictionary addModelDict(HFDIBDEMDict_.subDict(bodyName).subDict("bodyAddition"));
        word howToAdd(word(addModelDict.lookup("addModel")));
        word bodyGeom;
        std::unique_ptr<geomModel> bodyGeomModel;
        // check if the immersedDict_ contains bodyGeom
        if (HFDIBDEMDict_.subDict(bodyName).found("bodyGeom"))
        {

            word input = word(HFDIBDEMDict_.subDict(bodyName).lookup("bodyGeom"));
            bodyGeom = input;
            InfoH << iB_Info << "Found bodyGeom for " << bodyName
                << ", the body is: " << bodyGeom << endl;
        }
        else
        {
            bodyGeom = "convex";
            InfoH << iB_Info << "Did not find bodyGeom for "
                << bodyName << ", using bodyGeom: " << bodyGeom << endl;
        }

        if(bodyGeom == "convex")
        {
            word stlPath("constant/triSurface/"+bodyName+".stl");
            bodyGeomModel = std::unique_ptr<convexBody>(
                new convexBody(mesh_,stlPath,thrSurf));
        }
        else if(bodyGeom == "nonConvex")
        {
            word stlPath("constant/triSurface/"+bodyName+".stl");
            bodyGeomModel = std::unique_ptr<nonConvexBody>(
                new nonConvexBody(mesh_,stlPath,thrSurf));
        }
        else if(bodyGeom == "sphere")
        {

            vector startPosition = vector(HFDIBDEMDict_.subDict(bodyName).subDict("sphere").lookup("startPosition"));
            scalar radius = readScalar(HFDIBDEMDict_.subDict(bodyName).subDict("sphere").lookup("radius"));

            bodyGeomModel = std::unique_ptr<sphereBody>(
                new sphereBody(mesh_,startPosition,radius,thrSurf));
        }
        else
        {
            word stlPath("constant/triSurface/"+bodyName+".stl");
            InfoH << iB_Info << "bodyGeom: " << bodyGeom
                << " not supported, using bodyGeom nonConvex" << endl;
            bodyGeom = "nonConvex";
            bodyGeomModel = std::unique_ptr<nonConvexBody>(
                new nonConvexBody(mesh_,stlPath,thrSurf));
        }

        InfoH << addModel_Info << "-- addModelMessage-- "
            << "trying to set addModel: " << howToAdd << endl;

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
                    std::move(bodyGeomModel),
                    cellPoints_
                )
            );
        }
        else if (howToAdd == "onceFromFile")
        {
            addModels_.set
            (
                nameI,
                new addModelOnceFromFile
                (
                    addModelDict,
                    mesh_,
                    startTime0,
                    std::move(bodyGeomModel),
                    cellPoints_,
                    bodyGeom,
                    thrSurf
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
                    std::move(bodyGeomModel),
                    cellPoints_
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
                    std::move(bodyGeomModel),
                    cellPoints_
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
                    std::move(bodyGeomModel),
                    cellPoints_
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
                    std::move(bodyGeomModel),
                    cellPoints_
                )
            );
        }
        else
        {
            InfoH << addModel_Info << "-- addModelMessage-- "
                << "No addition model specified, defaulting to addModelOnce"
                << endl;

            addModels_.set
            (
                nameI,
                new addModelOnce
                (
                    addModelDict,
                    mesh_,
                    startTime0,
                    std::move(bodyGeomModel),
                    cellPoints_
                )
            );
        }
    }
