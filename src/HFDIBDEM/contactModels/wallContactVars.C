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

    openHFDIB-DEM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (Version 3) as published
    by the Free Software Foundation.

    openHFDIB-DEM is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with openHFDIB-DEM. If not, see <http://www.gnu.org/licenses/>.

InNamespace
    Foam

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "wallContactVars.H"

using namespace Foam;

//---------------------------------------------------------------------------//
void wallContactVars::setMeanCntPars
(
    const fvMesh&   mesh,
    DynamicList<Tuple2<label,string>>& contactFaces,
    HashTable<physicalProperties,string,Hash<string>>& wallMeanPars
)
{
    scalar overallArea = 0;
    physicalProperties_.aY_ = 0;
    physicalProperties_.aG_ = 0;
    physicalProperties_.aMu_ = 0;
    physicalProperties_.maxAdhN_ = 0;
    physicalProperties_.reduceBeta_ = 0;

    forAll(contactFaces,faceI)
    {
        scalar area = mesh.magSf()[contactFaces[faceI].first()];
        overallArea += area;

        physicalProperties& cMeanCntPars(
            wallMeanPars[contactFaces[faceI].second()]
        );

        physicalProperties_.aY_ += (cMeanCntPars.aY_*area);
        physicalProperties_.aG_ += (cMeanCntPars.aG_*area);
        physicalProperties_.aMu_ += (cMeanCntPars.aMu_*area);
        physicalProperties_.maxAdhN_ += (cMeanCntPars.maxAdhN_*area);
        physicalProperties_.reduceBeta_ += (cMeanCntPars.reduceBeta_*area);
    }

    physicalProperties_.aY_ /= overallArea;
    physicalProperties_.aG_ /= overallArea;
    physicalProperties_.aMu_ /= overallArea;
    physicalProperties_.maxAdhN_ /= overallArea;
    physicalProperties_.reduceBeta_ /= overallArea;
}
//---------------------------------------------------------------------------//
void wallContactVars::setMeanCntPars_Plane
(
    List<scalar>& contactAreas,
    List<string> contactFaces,
    HashTable<physicalProperties,string,Hash<string>>& wallMeanPars
)
{
    scalar overallArea = 0;
    physicalProperties_.aY_ = 0;
    physicalProperties_.aG_ = 0;
    physicalProperties_.aMu_ = 0;
    physicalProperties_.maxAdhN_ = 0;
    physicalProperties_.reduceBeta_ = 0;

    forAll(contactFaces,faceI)
    {
        scalar area = contactAreas[faceI];
        overallArea += area;

        physicalProperties& cMeanCntPars(
            wallMeanPars[contactFaces[faceI]]
        );

        physicalProperties_.aY_ += (cMeanCntPars.aY_*area);
        physicalProperties_.aG_ += (cMeanCntPars.aG_*area);
        physicalProperties_.aMu_ += (cMeanCntPars.aMu_*area);
        physicalProperties_.maxAdhN_ += (cMeanCntPars.maxAdhN_*area);
        physicalProperties_.reduceBeta_ += (cMeanCntPars.reduceBeta_*area);
    }

    physicalProperties_.aY_ /= overallArea;
    physicalProperties_.aG_ /= overallArea;
    physicalProperties_.aMu_ /= overallArea;
    physicalProperties_.maxAdhN_ /= overallArea;
    physicalProperties_.reduceBeta_ /= overallArea;
}
//---------------------------------------------------------------------------//
