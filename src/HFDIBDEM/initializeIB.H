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
    InfoH << iB_Info << bodyName_ << " is a static body." << endl;
}
else if (immersedDict_.found("prescribedTransBody"))
{
    bodyOperation_ = 1;
    Vel_   = vector(immersedDict_.subDict("prescribedTransBody").lookup("velocity"));

    InfoH << iB_Info << bodyName_ << " is a freely rotating body with prescribed linear velocity." << endl;
}
else if (immersedDict_.found("prescribedRotBody"))
{
    bodyOperation_ = 2;
    Axis_ = vector(immersedDict_.subDict("prescribedRotBody").lookup("axis"));

    omega_  = readScalar(
                immersedDict_.subDict("prescribedRotBody").lookup("omega")
            );

    InfoH << iB_Info << bodyName_ << " is a freely moving body with prescribed rotation." << endl;
}
else if (immersedDict_.found("prescribedTransRotBody"))
{
    bodyOperation_ = 3;
    Vel_   = vector(immersedDict_.subDict("prescribedTransRotBody").lookup("velocity"));
    Axis_  = vector(immersedDict_.subDict("prescribedTransRotBody").lookup("axis"));
    omega_ = readScalar(
                immersedDict_.subDict("prescribedTransRotBody").lookup("omega")
            );

    InfoH << iB_Info << bodyName_ << " has prescribed both movement and rotation." << endl;
}
else if (immersedDict_.found("prescribedTransFixedAxisRotBody"))
{
    bodyOperation_ = 4;
    Vel_  = vector(immersedDict_.subDict("prescribedTransFixedAxisRotBody").lookup("velocity"));
    Axis_ = vector(immersedDict_.subDict("prescribedTransFixedAxisRotBody").lookup("axis"));

    InfoH << iB_Info << bodyName_ << " has prescribed movement and axis of rotation." << endl;
}
else if (immersedDict_.found("fullyCoupledBody"))
{
    bodyOperation_ = 5;

    if (immersedDict_.isDict("fullyCoupledBody"))
    {
        if(immersedDict_.subDict("fullyCoupledBody").found("velocity"))
        {
            Vel_ = vector(immersedDict_.subDict("fullyCoupledBody").lookup("velocity"));
        }
        if(immersedDict_.subDict("fullyCoupledBody").found("axis"))
        {
            Axis_ = vector(immersedDict_.subDict("fullyCoupledBody").lookup("axis"));
        }
        if(immersedDict_.subDict("fullyCoupledBody").found("omega"))
        {
            omega_  = readScalar(
                immersedDict_.subDict("fullyCoupledBody").lookup("omega")
            );
        }
    }

    InfoH << iB_Info << bodyName_ << " is fully coupled with fluid phase." << endl;
}
else
{
    InfoH << iB_Info << "No body operation was found for " << bodyName_ << endl
      << "Assuming static body.";
}

// check if the immersedDict_ contains switch for updateTorque_
if (immersedDict_.found("updateTorque"))
{
    updateTorque_ = readBool(immersedDict_.lookup("updateTorque"));
    InfoH << iB_Info << "Found updateTorque: " << updateTorque_ << endl;
}
else
{
    InfoH << iB_Info << "Did not find updateTorque, using updateTorque: " << updateTorque_ << endl;
}

// do I want the IB to start as in sync with the flow?
if (immersedDict_.found("startSynced"))
{
    startSynced_ = readBool(immersedDict_.lookup("startSynced"));
    if (startSynced_)
    {
        InfoH << iB_Info << "Will try to sync the body with the flow upon creation" << endl;
    }
    else
    {
        InfoH << iB_Info << "The body will be created as static" << endl;
    }
}
else
{
    InfoH << iB_Info << "startSynced was not specified, using startSynced: " << startSynced_ << endl;
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
if (immersedDict_.found("refineBuffers"))
{
    refineBuffers_ = readLabel(immersedDict_.lookup("refineBuffers"));
}
if (immersedDict_.found("timesToSetStatic"))
{
    timesToSetStatic_ = readLabel(immersedDict_.lookup("timesToSetStatic"));
}

// set up body geometry
geomModel_->setIntSpan(intSpan_);
geomModel_->setSdBasedLambda(sdBasedLambda_);
geomModel_->setRhoS(dimensionedScalar(immersedDict_.get<dimensionedScalar>("rho")));

string material = word(immersedDict_.lookup("material"));

if(ibIntp.valid())
{
    intpInfo_.reset(ibIntp().createIntpInfo(                            //OF.com: set -> reset
        mesh_,
        geomModel_
    ));
}

ibContactClass_.reset(new ibContactClass(
    geomModel_,
    material
));

ibCntVars_.reset(new ibContactVars
(
    bodyId_,
    Vel_,
    omega_,
    Axis_,
    geomModel_->getM0(),
    geomModel_->getM(),
    geomModel_->getRhoS()
));

wallCntInfo_.reset(new wallContactInfo
(
    ibContactClass_(),
    ibCntVars_()
));
