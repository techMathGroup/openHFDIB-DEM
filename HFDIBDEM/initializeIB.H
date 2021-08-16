if (mesh_.nGeometricD() < 3)
    {
        const vector validDirs = (mesh_.geometricD() + Vector<label>::one)/2;
        Axis_ -= validDirs;
    }
    Axis_ /= mag(Axis_);
    
    // read the body operation name from immersedDict_
    if (immersedDict_.found("staticBody"))
    {
        bodyOperation_ = 0;
        Info << bodyName_ << " is a static body." << endl;
    }
    else if (immersedDict_.found("prescribedTransBody"))
    {
        bodyOperation_ = 1;
        
        Vel_   = immersedDict_.subDict("prescribedTransBody").lookup("velocity");
        
        Info << bodyName_ << " is a freely rotating body with prescribed linear velocity." << endl;
    }
    else if (immersedDict_.found("prescribedRotBody"))
    {
        bodyOperation_ = 2;
        
        Axis_ = immersedDict_.subDict("prescribedRotBody").lookup("axis");       
        omega_  = readScalar(
                    immersedDict_.subDict("prescribedRotBody").lookup("omega")
                );
    
        Info << bodyName_ << " is a freely moving body with prescribed rotation." << endl;
    }
    else if (immersedDict_.found("prescribedTransRotBody"))
    {
        bodyOperation_ = 3;
        
        Vel_   = immersedDict_.subDict("prescribedTransRotBody").lookup("velocity");
        Axis_ = immersedDict_.subDict("prescribedTransRotBody").lookup("axis");       
        omega_  = readScalar(
                    immersedDict_.subDict("prescribedTransRotBody").lookup("omega")
                );
    
        Info << bodyName_ << " has prescribed both movement and rotation." << endl;
    }
    else if (immersedDict_.found("prescribedTransFixedAxisRotBody"))
    {
        bodyOperation_ = 4;
        
        Vel_   = immersedDict_.subDict("prescribedTransFixedAxisRotBody").lookup("velocity");
        Axis_ = immersedDict_.subDict("prescribedTransFixedAxisRotBody").lookup("axis");
        
        Info << bodyName_ << " has prescribed movement and axis of rotation." << endl;
    }
    else if (immersedDict_.found("fullyCoupledBody"))
    {
        bodyOperation_ = 5;
        
        Info << bodyName_ << " is fully coupled with fluid phase." << endl;
    }
    else
    {
        Info << "No body operation was found for " << bodyName_ << endl
          << "Assuming static body.";
    }
    
    // check if the immersedDict_ contains switch for updateTorque_
    if (immersedDict_.found("updateTorque"))
    {
        updateTorque_ = readBool(immersedDict_.lookup("updateTorque"));
        Info << "Found updateTorque: " << updateTorque_ << endl;
    }
    else
    {
        Info << "Did not find updateTorque, using updateTorque: " << updateTorque_ << endl;
    }
    
    // do I want the IB to start as in sync with the flow?
    if (immersedDict_.found("startSynced"))
    {
        startSynced_ = readBool(immersedDict_.lookup("startSynced"));
        if (startSynced_)
        {
            Info << "Will try to sync the body with the flow upon creation" << endl;
        }
        else
        {
            Info << "The body will be created as static" << endl;
        }
    }
    else
    {
        Info << "startSynced was not specified, using startSynced: " << startSynced_ << endl;
    }
    
    // body auxiliary settings
    if (immersedDict_.found("sdBasedLambda"))
    {
        sdBasedLambda_ = readBool(immersedDict_.lookup("sdBasedLambda"));
    }
    if (immersedDict_.found("interfaceSpan"))
    {
        intSpan_ = readScalar(immersedDict_.lookup("interfaceSpan"));
    }
    if (immersedDict_.found("velRelaxFac"))
    {
        velRelaxFac_ = readScalar(immersedDict_.lookup("velRelaxFac"));
    }
    if (immersedDict_.found("refineBuffers"))
    {
        refineBuffers_ = readLabel(immersedDict_.lookup("refineBuffers"));
    }
    if (immersedDict_.found("useInterpolation"))
    {
        useInterpolation_ = readBool(immersedDict_.lookup("useInterpolation"));
    }
    if (immersedDict_.found("timesToSetStatic"))
    {
        timesToSetStatic_ = readLabel(immersedDict_.lookup("timesToSetStatic"));
    }
    
    // set up body geometry
    geomModel_->setIntSpan(intSpan_);
    geomModel_->setSdBasedLambda(sdBasedLambda_);
    geomModel_->setRhoS(immersedDict_.lookup("rho"));
    
    // set up contact properties
    scalar kN(readScalar(immersedDict_.lookup("kN")));
    scalar kt(readScalar(immersedDict_.lookup("kt")));
    scalar gammaN(readScalar(immersedDict_.lookup("gammaN")));
    scalar gammat(readScalar(immersedDict_.lookup("gammat")));
    scalar mu(readScalar(immersedDict_.lookup("mu")));
    scalar adhN(0.0);
    if (immersedDict_.found("adhN"))
    {
        adhN = readScalar(immersedDict_.lookup("adhN"));
    }
    scalar adhEqui(0.0);
    if (immersedDict_.found("adhEqui"))
    {
        adhEqui = readScalar(immersedDict_.lookup("adhEqui"));
    }
    
    contactInfo_.set(new contactInfo(
        geomModel_(),
        surfCells_,
        intCells_,
        kN,
        kt,
        gammaN,
        gammat,
        mu,
        adhN,
        adhEqui)
    );

    // set sizes for parallel runs
    surfCells_.setSize(Pstream::nProcs());
    intCells_.setSize(Pstream::nProcs());
    interpolationInfo_.setSize(Pstream::nProcs());
    interpolationVecReqs_.setSize(Pstream::nProcs());
    
    // initialize bounding points
    boundBox bound(geomModel_->getBounds());
    minBoundPoint_ = bound.min();
    maxBoundPoint_ = bound.max();    
    
    // initialize list of labels for bounding points based on bodyID
    // -> this list holds index of both bounding points for each dimension
    // -> thanks to this information the update of bounding points positions is more efficient
    for (label coord = 0; coord < 3; coord = coord + 1)
    {
        List<label> labeli(2);
        labeli[0] = 2 * bodyId_;
        labeli[1] = 2 * bodyId_ + 1;
        boundIndList_[coord] = labeli;
    } 
