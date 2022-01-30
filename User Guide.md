# Documentation for HFDIBDEMdict

**bodyNames** - *required* > List of body names (***bodyName***)

**surfaceThreshold** - *required* > Cutoff value for determining body presence in a cell (based on lambda fraction)

**stepDEM** - *required* > Fraction of timestep used as a sub timestep for DEM calculation

**geometricD** - *optional* > Define empty direction for pseudo 2D simulation.
> *default value is taken from fvMesh i.e. empty patch is recognized.*

**recordSimulation** - *required* > Should be particles recorded into bodiesInfo directory
> Possible values: {*true*, *false*}

**interpolationSchemes** - *optional* > Interpolation settings used for velocity reconstruction near bodies

> **U** - *required* > Scheme for velocity interpolation inside a cell. Values according to Foam::interpolation
>> Possible values: {*cell*, *cellPoint*, *cellPointFace*}

> **method** - *required* > Method for interpolation point searching
>> Possible values: {*line*}

**outputSetup** - *optional* > Output settings. When not provided, everything is outputted.

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

>> ***materialName*** - *required* > Custom name of material. Multiple blocks possible (see example)
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

> ***bodyType*** - *required* > type of body
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