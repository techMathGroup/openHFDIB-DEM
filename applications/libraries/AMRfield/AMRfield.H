#ifndef AMRfield_H
#define AMRfield_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*/
//                         Class AMRfield Declaration
/*---------------------------------------------------------------------------*/

class AMRfield
:
    public fvMeshFunctionObject
{
public:

    TypeName("AMRfield");

    AMRfield
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    virtual ~AMRfield() = default;

    // Member functions
    virtual bool read(const dictionary& dict);
    virtual bool execute();
    virtual bool write();
};

} 
} 

#endif
