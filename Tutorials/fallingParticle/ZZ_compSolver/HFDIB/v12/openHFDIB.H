/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________
                       | | | ||  ___|  _  \_   _| ___ \     H ybrid
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /     F ictitious
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \     D omain
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /     I mmersed
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/      B oundary
      | |
      |_|
-------------------------------------------------------------------------------
License

    openHFDIB is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

Description
    implementation of the HFDIBM method (Municchi and Radl, 2016) in OpenFOAM
    
Note: (MI)
    there are significant problems with dynamic meshes - the program
    crashes or freezes. I believe it has something to do with unlucky
    treatment of dynamically allocated lists and pointers
    (after an object destruction, something lives on and is shared with
    another object of the same instance)
Note: (MI)
    I removed almost all the pointers and objects with dynamically
    allocated memory from the code. I would say that it resulted in a
    slightly better stability and speed and it now works even on
    dynamic meshes (a strong statement which is probably not completely
    true)

SourceFiles
    openHFDIB.C

Contributors
    Federico Municchi (2016), Martin Isoz (2019)
\*---------------------------------------------------------------------------*/

#ifndef openHFDIB_H
#define openHFDIB_H

#include "dictionary.H"
#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedTensor.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "fvMatricesFwd.H"
//~ #include "fvMesh.H"
#include "dynamicFvMesh.H"
#include "triSurfaceMesh.H"
//~ #include <vector>
#include "immersedBody.H"

namespace Foam
{

 class openHFDIB
 {

    private:
  
    // Private data types

    typedef DynamicList<label>  DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point>  DynamicPointList;

    //Basic members
    const  dynamicFvMesh&            mesh_;
    IOdictionary                HFDIBDict_;
    dictionary            HFDIBinterpDict_;
    IOdictionary      transportProperties_;

    //Immersed bodies objects
    PtrList<immersedBody> immersedBodies_;
    
    //Contact variables (global)
    
    
    public:
    
    openHFDIB(const Foam::dynamicFvMesh& mesh);
    ~openHFDIB();
    
    void initialize();
    void update
    (
        volScalarField& body,
        volVectorField& f
    );
    
    //Interpolate vector V using Vs for the surface values
    void interpolateIB
    (
        volVectorField& V,
        volVectorField& Vs,
        volScalarField& body
    );
    
    //Interpolate scalar S using Ss for the surface values
    void interpolateIB
    (
        volScalarField& S,
        volScalarField& Ss,
        volScalarField& body
    )
    {};
    
    // Reinitialize in the case of a mesh change
    void forceReinit
    (
        volScalarField & body,
        volVectorField & f
    );

 };

}

#endif
