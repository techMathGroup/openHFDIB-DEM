/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  | |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___/  |_____|_|  |_|
      | |                     H ybrid F ictitious D omain - I mmersed B oundary
      |_|                                        and D iscrete E lement M ethod
-------------------------------------------------------------------------------
License

    openHFDIB-DEM is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

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
    class for storage of contact variables

SourceFiles
    wallContactInfo.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef wallContactInfo_H
#define wallContactInfo_H

#include "dictionary.H"
#include "fvCFD.H"

#include "geomModel.H"
#include "ibContactVars.H"
#include "ibContactClass.H"
#include "wallContactVars.H"
#include "materialInfo.H"
#include "physicalProperties.H"
#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class wallContactInfo Declaration
\*---------------------------------------------------------------------------*/

class wallContactInfo
{
    private:

    typedef DynamicList<label> DynamicLabelList;
    typedef DynamicList<scalar> DynamicScalarList;
    typedef DynamicList<vector> DynamicVectorList;
    typedef DynamicList<point> DynamicPointList;

    ibContactClass& ibContactClass_;
    ibContactVars& ibContactVars_;

    label bodyId_;

    // contacts with walls
    DynamicLabelList curUsedHashes_;
    DynamicLabelList wallCntHashes_;
    HashTable<wallContactVars,label,Hash<label>> wallcVarsTable_;

    // Mean model parameters
    HashTable<physicalProperties,string,Hash<string>> wallMeanPars_;

    scalar reduceM_;    // Reduced mass

    forces outForce_; // Out variable for force

    public:

    // Constructors ------------------------------------------------
    wallContactInfo
    (
        ibContactClass& cClass,
        ibContactVars& cVars
    );

    // Destructors -------------------------------------------------
    ~wallContactInfo();

    ibContactClass& getcClass()
    {
        return ibContactClass_;
    }

    ibContactVars& getcVars()
    {
        return ibContactVars_;
    }

    label getBodyId()
    {
        return bodyId_;
    }

    wallContactVars& getWallCntVars
    (
        const vector contactCenter,
        const vector contactNormal,
        const scalar deltaT
    );

    void clearWallCntVars();

    DynamicLabelList& getCurUsedHashes()
    {
        return curUsedHashes_;
    }

    DynamicLabelList& getWallCntHashes()
    {
        return wallCntHashes_;
    }

    HashTable<wallContactVars,label,Hash<label>>& getWallcVarsTable()
    {
        return wallcVarsTable_;
    }

    scalar getMu(wallContactVars& wallCntvar)
    {
        return wallCntvar.getMeanCntPar().aMu_;
    }

    vector getLVec(wallContactVars& wallCntvar);

    vector getVeli(wallContactVars& wallCntvar);

    void evalVariables(wallContactVars& wallCntvar);

    vector getFNe(wallContactVars& wallCntvar);

    vector getFA(wallContactVars& wallCntvar);

    vector getFNd(wallContactVars& wallCntvar);

    vector getFt(wallContactVars& wallCntvar, scalar deltaT);

    forces& getOutForce()
    {
        return outForce_;
    }

    bool isContactPatch(string patchName)
    {
        return wallMeanPars_.found(patchName);
    }

    HashTable<physicalProperties,string,Hash<string>>& getWallMeanPars()
    {
        return wallMeanPars_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
