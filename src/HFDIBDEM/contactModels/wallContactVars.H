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
    class for storage of particle-wall contact variables

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef wallContactVars_H
#define wallContactVars_H

#include "physicalProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class wallContactVars Declaration
\*---------------------------------------------------------------------------*/

class wallContactVars
{
    public:
    vector contactCenter_;
    scalar contactVolume_;
    vector contactNormal_;
    scalar contactArea_;

    physicalProperties physicalProperties_;

    vector lVec_;
    vector Veli_;
    scalar Vn_;
    scalar Lc_;
    scalar curAdhN_;    // Current contact adhesive coef
    
    vector lVector_;    // Vector from contact center to the contact point

    vector FtPrev_;     // Tangential force from previous time-step

    // Constructors ------------------------------------------------
    wallContactVars()
    {
        contactCenter_ = vector::zero;
        contactVolume_ = 0;
        contactNormal_ = vector::zero;
        contactArea_ = 0;

        lVec_ = vector::zero;
        Veli_ = vector::zero;
        Vn_ = 0;
        Lc_ = 0;
        curAdhN_ = 0;
        FtPrev_ = vector::zero;
    }

    // Destructors -------------------------------------------------
    ~wallContactVars(){}

    void setMeanCntPars
    (
        const fvMesh&   mesh,
        DynamicList<Tuple2<label,string>>& contactFaces,
        HashTable<physicalProperties,string,Hash<string>>& wallMeanPars
    );

    void setMeanCntPars_Plane
    (
        List<scalar>& contactAreas,
        List<string> contactFaces,
        HashTable<physicalProperties,string,Hash<string>>& wallMeanPars
    );

    physicalProperties& getMeanCntPar()
    {
        return physicalProperties_;
    }

    void clear()
    {
        contactCenter_ = vector::zero;
        contactVolume_ = 0;
        contactNormal_ = vector::zero;
        contactArea_ = 0;

        lVec_ = vector::zero;
        Veli_ = vector::zero;
        Vn_ = 0;
        Lc_ = 0;
        curAdhN_ = 0;
        FtPrev_ = vector::zero;
    }    
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
