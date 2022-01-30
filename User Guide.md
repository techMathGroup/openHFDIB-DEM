# Documentation for HFDIBDEMdict

**bodyNames** - *required* > List of body names (***bodyName***)

**surfaceThreshold** - *required* > Cutoff value for determining body presence in a cell (based on lambda fraction)

**stepDEM** - *required* > Fraction of timestep used as a sub timestep for DEM calculation

**geometricD** - *optional* > Define empty direction for pseudo 2D simulation
> *default value is taken from fvMesh i.e. empty patch is recognized.*

**recordSimulation** - *required* > Should be particles recorded into bodiesInfo directory
> Possible values: {*true*, *false*}

**interpolationSchemes** - *optional* > Interpolation settings used for velocity reconstruction near bodies

> **U** - *required* > Scheme for velocity interpolation inside a cell. Values according to Foam::interpolation
>> Possible values: {*cell*, *cellPoint*, *cellPointFace*}

> **method** - *required* > Method for interpolation point searching
>> Possible values: {*line*}

**outputSetup** - *optional* > Output settings. When not provided, everything is outputted

> **basic** - *required* > Basic output
>> Possible values: {*true*, *false*}

> **iB** - *required* > Detail output for bodies
>> Possible values: {*true*, *false*}

> **DEM** - *required* > Output for DEM
>> Possible values: {*true*, *false*}

> **addModel** - *required* > Output for body addition
>> Possible values: {*true*, *false*}

**DEM** - *required* > Input for DEM

> **materials** - *required* > Input for materials  
>> ***materialName*** - *required* > Custom name of material Multiple blocks possible (see example)
>>> **Y** - *required* > Young's modulus  
>>> **nu** - *required* > Poisson ratio  
>>> **gamma** - *required* > Viscoelastic damping constant  
>>> **mu** - *required* > Tangential force truncation  
>>> **adhN** - *required* > Adhesive force coefficient

> **interfaceAdh** - *optional* > Truncation for adhesive for between materials  
>> ***name*** - *required* > Custom id name
>>> **materials** - *required* > Two names of affected materials  
>>> **value** - *required* > Value is used for both materials

> **collisionPatches** - *required* > List of patches with which the bodies collide  
>> ***patchName*** - *required* > Name of the patch
>>> Possible values: ***materialName***

***bodyName*** - *required* > Body name which correspond to value in **bodyNames**

> ***bodyType*** - *required* > Type of body
>> Possible values:
>> - *staticBody* > Body is static
>> - *prescribedTransBody* > Body has prescribed translational movement
>>> **velocity** - *required* > Translational velocity. (*in subdictionary*)
>> - *prescribedRotBody* > Body has prescribed rotational movement
>>> **axis** - *required* > Axis of rotation. (*in subdictionary*)  
>>> **omega** - *required* > Angular velocity . (*in subdictionary*)
>> - *prescribedTransRotBody* > Body has prescribed translational and rotational movement
>>> **velocity** - *required* > Translational velocity. (*in subdictionary*)  
>>> **axis** - *required* > Axis of rotation. (*in subdictionary*)  
>>> **omega** - *required* > Angular velocity . (*in subdictionary*)
>> - *prescribedTransFixedAxisRotBody* > Body has prescribed translational movement. Rotational movement has fixed axis.
>>> **velocity** - *required* > Translational velocity. (*in subdictionary*)  
>>> **axis** - *required* > Axis of rotation. (*in subdictionary*)
>> - *fullyCoupledBody* > Body is fully coupled with fluid phase
>>> **velocity** - *optional* > Translational velocity. (*in subdictionary*)

> **rho** - *required* > Density of the body

> **U** - *required* > Dictionary for boundary condition
>> **BC** - *required* > Name of boundary condition
>>> Possible values: {*noSlip*}

> **material** - *required* > Material of the body
>> Possible values: ***materialName***

> **bodyGeom** - *required* > Body geometry type
>> Possible values:
>> - *convex* > Body is defined according to STL file: ./constant/triSurface/***bodyName***.stl. Solver is optimized for this type of geometry.
>> - *nonConvex* > Body is defined according to STL file: ./constant/triSurface/***bodyName***.stl
>> - *sphere* > Body is defined by center and radius.
>> **sphere** - *required*  
>>> **startPosition** - *required* > Starting position of the center. (*in subdictionary*. May not be used in some **addModel**)  
>>> **radius** - *required* > Radius of the sphere. (*in subdictionary*)

> **updateTorque** - *optional* > Should be rotational movement updated according to acting forces?
>> Possible values: {*true*, *false*}

> **sdBasedLambda** - *optional* > sdBased algorithm for lambda reconstruction. Default false
>> Possible values: {*true*, *false*}  

> **interfaceSpan** - *optional* > For sdBased algorithm. Best results with 1.0

> **startSynced** - *optional* > Should be body movement sync with fluid movement when added? Default false
>> Possible values: {*true*, *false*}

> **refineBuffers** - *optional* > Number of layers around cell to be refined. dynamicRefineFvMesh required when value > 0. Default 0

> **timesToSetStatic** - *optional* > Number of timesteps to set body as static. Value -1 means never. Default -1