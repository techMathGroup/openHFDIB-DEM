# AMR Library

This library provides adaptive mesh refinement (AMR) support for openHFDIB-DEM.
It must be compiled separately before running any case that uses AMR.

## Compilation

From this directory:

wmake libso


Make sure OpenFOAM v2412 (or compatible .com version) is sourced before compiling.

## Usage

Add the following to your case's `system/controlDict`:

libs ("$FOAM_USER_LIBBIN/libAMRfield.so");

functions
{
AMRfield
{
type AMRfield;
libs ("$FOAM_USER_LIBBIN/libAMRfield.so");
executeControl timeStep;
executeInterval 1;
writeControl writeTime;
}
}