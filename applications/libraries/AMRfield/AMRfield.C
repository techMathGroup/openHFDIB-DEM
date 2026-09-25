/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Library
    AMRfield

Description
    AMR field for refineCombined in pimplHFDIBDEM

\*---------------------------------------------------------------------------*/

#include "AMRfield.H"
#include "fvcCurl.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

defineTypeNameAndDebug(AMRfield, 0);
addToRunTimeSelectionTable(functionObject, AMRfield, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

AMRfield::AMRfield
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{}

bool AMRfield::read(const dictionary& dict)
{
    return fvMeshFunctionObject::read(dict);
}

bool AMRfield::execute()
{
    const volVectorField& U =
        mesh_.lookupObject<volVectorField>("U");

    volScalarField magVort = mag(fvc::curl(U));

    dimensionedScalar maxVal
    (
        "maxVal",
        magVort.dimensions(),
        max(gMax(magVort), SMALL)
    );

    // Create or update normMagVort
    if (!mesh_.foundObject<volScalarField>("normMagVort"))
    {
        auto* normMagVortPtr = new volScalarField
        (
            IOobject
            (
                "normMagVort",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT, //NO_READ,
                IOobject::AUTO_WRITE
            ),
            magVort / maxVal
        );
        normMagVortPtr->store();
    }
    else
    {
        volScalarField& normMagVort =
            mesh_.lookupObjectRef<volScalarField>("normMagVort");
        normMagVort = magVort / maxVal;
        normMagVort.correctBoundaryConditions();
    }

    const volScalarField& normMagVort =
        mesh_.lookupObject<volScalarField>("normMagVort");

    const volScalarField& refineF =
        mesh_.lookupObject<volScalarField>("refineF");

    // Create or update refineCombined
    if (!mesh_.foundObject<volScalarField>("refineCombined"))
    {
        auto* refineCombinedPtr = new volScalarField
        (
            IOobject
            (
                "refineCombined",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT, //NO_READ,
                IOobject::AUTO_WRITE
            ),
            max(normMagVort, refineF)
        );
        refineCombinedPtr->store();
    }
    else
    {
        volScalarField& refineCombined =
            mesh_.lookupObjectRef<volScalarField>("refineCombined");
        refineCombined = max(normMagVort, refineF);
        refineCombined.correctBoundaryConditions();
    }

    return true;
}

bool AMRfield::write()
{
    return true;
}

} 
} 

// ************************************************************************* //
