# Documentation for HFDIBDEMdict

**bodyNames** - *required* > List of body names

**surfaceThreshold** - *required* > Cutoff value for determining body presence in a cell (based on lambda fraction)

**stepDEM** - *required* > fraction of timestep used as a sub timestep for DEM calculation

**geometricD** - *optional* > define empty direction for pseudo 2D simulation.
> *default value is taken from fvMesh i.e. empty patch is recognized.*

**recordSimulation** - *required* > should be particles recorded into bodiesInfo directory
> Possible values: {*true*, *false*}

**interpolationSchemes** - *optional* > interpolation settings used for velocity reconstruction near bodies

> **U** - *required* > scheme for velocity interpolation inside a cell. Values according to Foam::interpolation
>> Possible values: {*cell*, *cellPoint*, *cellPointFace*}

> **method** - *required* > method for interpolation point searching
>> Possible values: {*line*}

**outputSetup** - *optional* > output settings. When not provided, everything is outputted.

> **basic** - *required* > basic output
>> Possible values: {*true*, *false*}

> **iB** - *required* > detail output for bodies
>> Possible values: {*true*, *false*}

> **DEM** - *required* > output for DEM
>> Possible values: {*true*, *false*}

> **addModel** - *required* > output for body addition
>> Possible values: {*true*, *false*}

**DEM** - *required* > Input for DEM

> **materials** - *required* > Input for materials