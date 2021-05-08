    // Note: if you want to add a new addModel, do it here
    
    forAll (stlNames_,nameI)
    {
        word stlName(stlNames_[nameI]);
        Info << "-- addModelMessage-- " << "Creating addition model for body: " << stlName << endl;
        
        // initialize the addition models
        dictionary addModelDict(HFDIBDEMDict_.subDict(stlName).subDict("bodyAddition"));
        word howToAdd(word(addModelDict.lookup("addModel")));
        
        Info << "-- addModelMessage-- " << "trying to set addModel: " << howToAdd << endl;
        
        if (howToAdd == "once")
        {
            
            addModels_.set
            (
                nameI,
                new addModelOnce
                (
                    addModelDict,
                    stlName,
                    mesh_,
                    startTime0
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
                    stlName,
                    mesh_,
                    startTime0,
                    geometricD_
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
                    stlName,
                    mesh_
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
                    stlName,
                    mesh_,
                    geometricD_
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
                    stlName,
                    mesh_,
                    geometricD_
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
                    stlName,
                    mesh_,
                    startTime0
                )
            );
        }
    }
