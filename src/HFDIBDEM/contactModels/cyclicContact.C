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
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-2025),
    Ondřej Studeník (2020-*), Lucie Kubíčková (2026-*)
\*---------------------------------------------------------------------------*/
#include "cyclicContact.H"

#include "cyclicPlaneInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace contactModel
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//---------------------------------------------------------------------------//
bool detectCyclicContact(
    wallContactInfo& wallCntInfo,
    vector& transVec
)
{
    // bool inContact = false;
    // int numOfContact = 0;

    // label nCells = mesh.nCells();
    // List<DynamicLabelList>& surfCells(wallCntInfo.getcClass().getSurfCells());
    // // go through all surfCells and check if there is any surfCell whose face is a boundary face
    // forAll (surfCells[Pstream::myProcNo()],sCellI)
    // {
    //     label cCell(surfCells[Pstream::myProcNo()][sCellI]);
    //     if(cCell < nCells)
    //     {
    //         const labelList& cFaces = mesh.cells()[cCell];
    //         forAll (cFaces,faceI)
    //         {
    //             if (!mesh.isInternalFace(cFaces[faceI]))
    //             {
    //                 // get reference to the patch which is in contact with IB. There is contact only if the patch is marked as a wall
    //                 label facePatchId = mesh.boundaryMesh().whichPatch(cFaces[faceI]);
    //                 const polyPatch& cPatch = mesh.boundaryMesh()[facePatchId];
    //                 if (cPatch.type()=="cyclic")
    //                 {
    //                     forAll (cyclicPatches, patchI)
    //                     {
    //                         if (cyclicPatches[patchI] == cPatch.name())
    //                         {
    //                             const cyclicPolyPatch& cyclicPatch = refCast<const cyclicPolyPatch>(cPatch);
    //                             transVec = cyclicPatch.transform().invTransformPosition(wallCntInfo.getcClass().getGeomModel().getCoM());
    //                             inContact = true;
    //                             numOfContact = 1;
    //                             break;
    //                         }
    //                     }
    //                 }
    //                 else if(cPatch.type()=="processorCyclic")
    //                 {
    //                     forAll (cyclicPatches, patchI)
    //                     {
    //                         if(cPatch.name().find(cyclicPatches[patchI]) != string::npos)
    //                         {
    //                             const coupledPolyPatch& cyclicPatch = refCast<const coupledPolyPatch>(cPatch);
    //                             transVec = cyclicPatch.transform().invTransformPosition(wallCntInfo.getcClass().getGeomModel().getCoM());
    //                             inContact = true;
    //                             numOfContact = 1;
    //                             break;
    //                         }
    //                     }
    //                 }
    //             }
    //             if(inContact)
    //                 break;
    //         }
    //     }
    //     if(inContact)
    //         break;
    // }

    if (wallCntInfo.detectWallContact(cyclicPlaneInfo::getCyclicPlaneInfo()))
    {
        transVec = cyclicPlaneInfo::getCyclicTransVec(wallCntInfo.getContactPatches()[0]);
        return true;
    }

    return false;
}
//---------------------------------------------------------------------------//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace contactModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
