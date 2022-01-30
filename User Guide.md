# Documentation for HFDIBDEMdict

**bodyNames** - *required* > List of body names

**surfaceThreshold** - *required* > Cutoff value for determining body presence in a cell (based on lambda fraction)

**stepDEM** - *required* > fraction of timestep used as a sub timestep for DEM calculation

**geometricD** - *optional* > define empty direction for pseudo 2D simulation.
> *default value is taken from fvMesh i.e. empty patch is recognized.*

**recordSimulation** - *required* > should be particles recorded into bodiesInfo directory
> Possible values: *true*/*false*