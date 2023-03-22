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

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "prtContact.H"

#include "virtualMeshLevel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace contactModel
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//---------------------------------------------------------------------------//
void findSubContacts(
    const fvMesh&   mesh,
    prtContactInfo& prtcInfo
)
{
    if (
        prtcInfo.getcClass().getGeomModel().getcType() == sphere
        &&
        prtcInfo.gettClass().getGeomModel().getcType() == sphere
        )
    {
        findSubContacts_Sphere(prtcInfo);
    }
    else if (
        prtcInfo.getcClass().getGeomModel().getcType() == cluster
        ||
        prtcInfo.gettClass().getGeomModel().getcType() == cluster
        )
    {
        findSubContacts_Cluster(mesh, prtcInfo);
    }
    else
    {
        findSubContacts_ArbShape(mesh, prtcInfo);
    }

    prtcInfo.swapSubContactLists();
}
//---------------------------------------------------------------------------//
void findSubContacts_Sphere(
    prtContactInfo& prtcInfo
)
{
    prtcInfo.setSubContacts_Sphere();
}
//---------------------------------------------------------------------------//
void findSubContacts_ArbShape(
    const fvMesh&   mesh,
    prtContactInfo& prtcInfo
)
{
    // DynamicLabelList    commonCells;

    // List<DynamicLabelList> cSurfCells(prtcInfo.getcClass().getSurfCells());
    if (Pstream::myProcNo() == 0)
    {
        scalar ranCellVol = mesh.V()[0];
        prtcInfo.setSubContacts_ArbShape(ranCellVol);
    }
    // List<DynamicLabelList> tSurfCells(prtcInfo.gettClass().getSurfCells());

    // // iterate over surfCells to find common cells
    // forAll (cSurfCells[Pstream::myProcNo()],cSCellI)
    // {
    //     forAll (tSurfCells[Pstream::myProcNo()],tSCellI)
    //     {
    //         if (mag(cSurfCells[Pstream::myProcNo()][cSCellI]-tSurfCells[Pstream::myProcNo()][tSCellI]) < SMALL)
    //         {
    //             commonCells.append(cSurfCells[Pstream::myProcNo()][cSCellI]);
    //         }
    //     }
    // }

    // prtcInfo.setSubContacts_ArbShape(
    //     mesh,
    //     detectContactCells(
    //         mesh,
    //         prtcInfo.getcClass().getGeomModel(),
    //         prtcInfo.gettClass().getGeomModel(),
    //         commonCells
    //     )
    // );
}
//---------------------------------------------------------------------------//
void findSubContacts_Cluster(
    const fvMesh&   mesh,
    prtContactInfo& prtcInfo
)
{
    std::vector<std::shared_ptr<geomModel>> cBodies;
    std::vector<std::shared_ptr<geomModel>> tBodies;

    if(prtcInfo.getcClass().getGeomModel().isCluster())
    {
        periodicBody& cCluster = dynamic_cast<periodicBody&>(prtcInfo.getcClass().getGeomModel());
        cBodies = cCluster.getClusterBodies();
    }
    else
    {
        cBodies.push_back(prtcInfo.getcClass().getGeomModelPtr());
    }

    if(prtcInfo.gettClass().getGeomModel().isCluster())
    {
        periodicBody& tCluster = dynamic_cast<periodicBody&>(prtcInfo.gettClass().getGeomModel());
        tBodies = tCluster.getClusterBodies();
    }
    else
    {
        tBodies.push_back(prtcInfo.gettClass().getGeomModelPtr());
    }

    for(std::shared_ptr<geomModel>& cgModel : cBodies)
    {
        for(std::shared_ptr<geomModel>& tgModel : tBodies)
        {
            ibContactClass cIbClassI(
                cgModel,
                prtcInfo.getcClass().getMatInfo().getMaterial()
            );

            ibContactClass tIbClassI(
                tgModel,
                prtcInfo.gettClass().getMatInfo().getMaterial()
            );

            prtContactInfo tmpPrtCntInfo
            (
                cIbClassI,
                prtcInfo.getcVars(),
                tIbClassI,
                prtcInfo.gettVars()
            );

            findSubContacts(mesh, tmpPrtCntInfo);
            tmpPrtCntInfo.swapSubContactLists();
            prtcInfo.getPrtSCList().swap(tmpPrtCntInfo.getPrtSCList());
        }
    }
}
//---------------------------------------------------------------------------//
List<DynamicList<label>> detectContactCells
(
    const fvMesh&   mesh,
    geomModel& cGeomModel,
    geomModel& tGeomModel,
    DynamicLabelList & commonCells
)
{
    labelHashSet checkedOctreeCells;

    autoPtr<DynamicLabelList> nextToCheck(
            new DynamicLabelList);

    autoPtr<DynamicLabelList> auxToCheck(
            new DynamicLabelList);

    DynamicLabelList subContactCells;

    List<DynamicLabelList> baseSubContactList;

    label iterMax(mesh.nCells());
    label iterCount(0);

    while((commonCells.size() > SMALL) && iterCount++ < iterMax)
    {
        label iterCount2(0);
        nextToCheck().clear();
        subContactCells.clear();

        nextToCheck().append(commonCells[0]);

        while ((nextToCheck().size() > 0) && iterCount2++ < iterMax)
        {
            auxToCheck().clear();
            forAll(nextToCheck(), nCell)
            {
                if (!checkedOctreeCells.found(nextToCheck()[nCell]))
                {
                    checkedOctreeCells.insert(nextToCheck()[nCell]);
                    if(isCellContactCell(
                        mesh,
                        cGeomModel,
                        tGeomModel,
                        nextToCheck()[nCell]))
                    {
                        subContactCells.append(nextToCheck()[nCell]);
                        auxToCheck().append(mesh.cellCells()[nextToCheck()[nCell]]);
                    }
                }
            }
            const autoPtr<DynamicLabelList> helpPtr(nextToCheck.ptr());
            nextToCheck.set(auxToCheck.ptr());
            auxToCheck = helpPtr;
        }

        if(subContactCells.size() > 0)
        {
            baseSubContactList.append(subContactCells);
        }

        // Remove checked cells from commonCells
        for (auto it = commonCells.begin(); it != commonCells.end();)
        {
            if (checkedOctreeCells.found(*it))
            {
                it = commonCells.erase(it);
            }
            else
            {
                ++it;
            }
        }
    }

    return baseSubContactList;
}
//---------------------------------------------------------------------------//
bool isCellContactCell
(
    const fvMesh&   mesh,
    geomModel& cGeomModel,
    geomModel& tGeomModel,
    label cellId
)
{
    const pointField& pp = mesh.points();

    const labelList& vertexLabels = mesh.cellPoints()[cellId];
    const pointField vertexPoints(pp,vertexLabels);

    const boolList cvertexesInside = cGeomModel.pointInside(vertexPoints);
    const boolList tvertexesInside = tGeomModel.pointInside(vertexPoints);

    bool tBody(false);
    bool cBody(false);

    forAll(tvertexesInside,vertexId)
    {
        if(tvertexesInside[vertexId])
        {
            tBody = true;
        }
        if(cvertexesInside[vertexId])
        {
            cBody = true;
        }

    }
    return tBody && cBody;
}
//---------------------------------------------------------------------------//
bool detectPrtPrtContact(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass,
    prtSubContactInfo& subCInfo
)
{
    if
    (
        cClass.getGeomModel().getcType() == sphere
        &&
        tClass.getGeomModel().getcType() == sphere
    )
    {
        return detectPrtPrtContact_Sphere
        (
            mesh,
            cClass,
            tClass
        );
    }
    else if
    (
        cClass.getGeomModel().getcType() == cluster
        ||
        tClass.getGeomModel().getcType() == cluster
    )
    {
        return detectPrtPrtContact_Cluster
        (
            mesh,
            cClass,
            tClass,
            subCInfo
        );
    }
    else
    {
        return detectPrtPrtContact_ArbShape
        (
            mesh,
            cClass,
            tClass,
            subCInfo
        );
    }
}
//---------------------------------------------------------------------------//
bool detectPrtPrtContact_ArbShape(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass,
    prtSubContactInfo& subCInfo
)
{
    autoPtr<virtualMeshInfo>& vmInfo = subCInfo.getVMInfo();
    if (!vmInfo.valid())
    {
        return false;
    }

    virtualMesh virtMesh(
        vmInfo(),
        cClass.getGeomModel(),
        tClass.getGeomModel()
    );

    return virtMesh.detectFirstContactPoint();
}
//---------------------------------------------------------------------------//
bool detectPrtPrtContact_Sphere
(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass
)
{
    if
    (
        mag(cClass.getGeomModel().getCoM()-tClass.getGeomModel().getCoM())
        <
        ((cClass.getGeomModel().getDC() / 2) + (tClass.getGeomModel().getDC() / 2))
    )
    {
        return true;
    }
    return false;
}
//---------------------------------------------------------------------------//
bool detectPrtPrtContact_Cluster
(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass,
    prtSubContactInfo& subCInfo
)
{
    std::vector<std::shared_ptr<geomModel>> cBodies;
    std::vector<std::shared_ptr<geomModel>> tBodies;

    if(cClass.getGeomModel().isCluster())
    {
        periodicBody& cCluster = dynamic_cast<periodicBody&>(cClass.getGeomModel());
        cBodies = cCluster.getClusterBodies();
    }
    else
    {
        cBodies.push_back(cClass.getGeomModelPtr());
    }

    if(tClass.getGeomModel().isCluster())
    {
        periodicBody& tCluster = dynamic_cast<periodicBody&>(tClass.getGeomModel());
        tBodies = tCluster.getClusterBodies();
    }
    else
    {
        tBodies.push_back(tClass.getGeomModelPtr());
    }

    for(std::shared_ptr<geomModel>& cgModel : cBodies)
    {
        for(std::shared_ptr<geomModel>& tgModel : tBodies)
        {
            ibContactClass cIbClassI(
                cgModel,
                cClass.getMatInfo().getMaterial()
            );

            ibContactClass tIbClassI(
                tgModel,
                tClass.getMatInfo().getMaterial()
            );

            prtSubContactInfo tmpSubCInfoI(
                subCInfo.getCPair(),
                subCInfo.getPhysicalProperties()
            );

            if (detectPrtPrtContact(
                mesh,
                cIbClassI,
                tIbClassI,
                subCInfo
            ))
            {
                return true;
            }
        }
    }
    return false;
}
//---------------------------------------------------------------------------//
void getPrtContactVars_ArbShape(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass,
    prtSubContactInfo& subCInfo
)
{
    autoPtr<virtualMeshInfo>& vmInfo = subCInfo.getVMInfo();
    if (!vmInfo.valid())
    {
        return;
    }

    virtualMesh virtMesh(
        vmInfo(),
        cClass.getGeomModel(),
        tClass.getGeomModel()
    );

    scalar intersectedVolume(0);
    vector contactCenter(vector::zero);
    vector normalVector = vector::zero;
    scalar contactArea(0);

    intersectedVolume = virtMesh.evaluateContact();

    if(virtMesh.getEdgeSVPoints().size() <= 4)
    {
        intersectedVolume = 0;
    }

    if (intersectedVolume > 0)
    {
        Tuple2<scalar,vector> surfaceAndNormal = virtMesh.get3DcontactNormalAndSurface();
        contactArea = surfaceAndNormal.first();
        normalVector = surfaceAndNormal.second();
        contactCenter = virtMesh.getContactCenter();
    }

    if (intersectedVolume > 0 && contactArea > 0)
    {
        subCInfo.getprtCntVars().contactCenter_ = contactCenter;
        subCInfo.getprtCntVars().contactVolume_ = intersectedVolume;
        subCInfo.getprtCntVars().contactNormal_ = normalVector;
        subCInfo.getprtCntVars().contactArea_ = contactArea;
    }
    else
    {
        subCInfo.getprtCntVars().contactCenter_ = vector::zero;
        subCInfo.getprtCntVars().contactVolume_ = 0;
        subCInfo.getprtCntVars().contactNormal_ = vector::zero;
        subCInfo.getprtCntVars().contactArea_ = 0;
    }
}
//---------------------------------------------------------------------------//
void getPrtContactVars_Sphere
(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass,
    prtSubContactInfo& subCInfo
)
{
    scalar cRadius = cClass.getGeomModel().getDC() / 2;
    scalar tRadius = tClass.getGeomModel().getDC() / 2;
    vector tCenter = tClass.getGeomModel().getCoM();

    vector centerDir = cClass.getGeomModel().getCoM()
                        - tCenter;

    scalar d = mag(centerDir);

    if(mag(centerDir) < SMALL || d > (cRadius + tRadius))
    {
        subCInfo.getprtCntVars().contactCenter_ = vector::zero;
        subCInfo.getprtCntVars().contactVolume_ = 0;
        subCInfo.getprtCntVars().contactNormal_ = vector::zero;
        subCInfo.getprtCntVars().contactArea_ = 0;
        return;
    }

    scalar xLength = (sqr(d) - sqr(cRadius) + sqr(tRadius))/(2*d);

    if(sqr(xLength) > sqr(tRadius))
    {
        subCInfo.getprtCntVars().contactCenter_
            = tCenter + (centerDir/d)*xLength;
        subCInfo.getprtCntVars().contactVolume_
            = (4/3)*Foam::constant::mathematical::pi*pow(tRadius,3);
        subCInfo.getprtCntVars().contactNormal_ = centerDir/d;
        subCInfo.getprtCntVars().contactArea_
            = Foam::constant::mathematical::pi*sqr(tRadius);
        return;
    }

    if(case3D)
    {

        scalar cSphCapV = (Foam::constant::mathematical::pi
                            *sqr(cRadius - d + xLength)
                            *(3*cRadius - (cRadius - d + xLength))) / 3;
        scalar tSphCapV = (Foam::constant::mathematical::pi
                            *sqr(tRadius - xLength)
                            *(3*tRadius - (tRadius - xLength))) / 3;

        subCInfo.getprtCntVars().contactCenter_
            = tCenter + (centerDir/d)*xLength;
        subCInfo.getprtCntVars().contactVolume_
            = cSphCapV + tSphCapV;
        subCInfo.getprtCntVars().contactNormal_
            = centerDir/d;
        subCInfo.getprtCntVars().contactArea_
            = Foam::constant::mathematical::pi
                *sqr(sqrt(sqr(tRadius) - sqr(xLength)));
    }
    else
    {
        boundBox meshBounds = mesh.bounds();
        scalar emptyLength = meshBounds.max()[emptyDim]
                            - meshBounds.min()[emptyDim];

        scalar cCirSeg = sqr(cRadius)*acos((d - xLength)/cRadius)
                         - (d - xLength)*sqrt(sqr(cRadius)
                         - sqr(d - xLength));
        scalar tCirSeg = sqr(tRadius)*acos(xLength/tRadius)
                         - xLength*sqrt(sqr(tRadius) - sqr(xLength));

        subCInfo.getprtCntVars().contactCenter_
            = tCenter + (centerDir/d)*xLength;
        subCInfo.getprtCntVars().contactVolume_
            = (cCirSeg + tCirSeg)*emptyLength;
        subCInfo.getprtCntVars().contactNormal_
            = centerDir/d;
        subCInfo.getprtCntVars().contactArea_
            = 2*sqrt(sqr(tRadius) - sqr(xLength))*emptyLength;
    }
}
//---------------------------------------------------------------------------//
void getPrtContactVars_Cluster
(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass,
    prtSubContactInfo& subCInfo
)
{
    subCInfo.getprtCntVars().contactCenter_ = vector::zero;
    subCInfo.getprtCntVars().contactVolume_ = 0;
    subCInfo.getprtCntVars().contactNormal_ = vector::zero;
    subCInfo.getprtCntVars().contactArea_ = 0;

    std::vector<std::shared_ptr<geomModel>> cBodies;
    std::vector<std::shared_ptr<geomModel>> tBodies;

    if(cClass.getGeomModel().isCluster())
    {
        periodicBody& cCluster = dynamic_cast<periodicBody&>(cClass.getGeomModel());
        cBodies = cCluster.getClusterBodies();
    }
    else
    {
        cBodies.push_back(cClass.getGeomModelPtr());
    }

    if(tClass.getGeomModel().isCluster())
    {
        periodicBody& tCluster = dynamic_cast<periodicBody&>(tClass.getGeomModel());
        tBodies = tCluster.getClusterBodies();
    }
    else
    {
        tBodies.push_back(tClass.getGeomModelPtr());
    }

    for(std::shared_ptr<geomModel>& cgModel : cBodies)
    {
        for(std::shared_ptr<geomModel>& tgModel : tBodies)
        {
            ibContactClass cIbClassI(
                cgModel,
                cClass.getMatInfo().getMaterial()
            );

            ibContactClass tIbClassI(
                tgModel,
                tClass.getMatInfo().getMaterial()
            );

            prtSubContactInfo tmpSubCInfoI(
                subCInfo.getCPair(),
                subCInfo.getPhysicalProperties()
            );

            getPrtContactVars(
                mesh,
                cIbClassI,
                tIbClassI,
                tmpSubCInfoI
            );

            if (tmpSubCInfoI.getprtCntVars().contactVolume_
                > subCInfo.getprtCntVars().contactVolume_)
            {
                subCInfo.getprtCntVars() = tmpSubCInfoI.getprtCntVars();
            }
        }
    }
}
//---------------------------------------------------------------------------//
void getPrtContactVars
(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass,
    prtSubContactInfo& subCInfo
)
{
    if
    (
        cClass.getGeomModel().getcType() == sphere
        &&
        tClass.getGeomModel().getcType() == sphere
    )
    {
        getPrtContactVars_Sphere(
            mesh,
            cClass,
            tClass,
            subCInfo
        );
    }
    else if
    (
        cClass.getGeomModel().getcType() == cluster
        ||
        tClass.getGeomModel().getcType() == cluster
    )
    {
        getPrtContactVars_Cluster(
            mesh,
            cClass,
            tClass,
            subCInfo
        );
    }
    else
    {
        getPrtContactVars_ArbShape(
            mesh,
            cClass,
            tClass,
            subCInfo
        );
    }
}
//---------------------------------------------------------------------------//
bool solvePrtContact(
    const fvMesh&   mesh,
    prtContactInfo& cInfo,
    prtSubContactInfo& subCInfo,
    scalar deltaT
)
{
    getPrtContactVars(
        mesh,
        cInfo.getcClass(),
        cInfo.gettClass(),
        subCInfo
    );

    if (!(subCInfo.getprtCntVars().contactVolume_ > 0))
    {
        return false;
    }

    InfoH << parallelDEM_Info << "-- Detected Particle-particle contact: -- body "
            << subCInfo.getCPair().first() << " & -- body "
            << subCInfo.getCPair().second() << endl;
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact cBody pos: "
            << cInfo.getcClass().getGeomModel().getCoM() << " & tBody pos: "
            << cInfo.gettClass().getGeomModel().getCoM() << endl;
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact center "
            << subCInfo.getprtCntVars().contactCenter_ << endl;
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact normal "
            << subCInfo.getprtCntVars().contactNormal_ << endl;
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact volume "
            << subCInfo.getprtCntVars().contactVolume_ << endl;
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact area "
            << subCInfo.getprtCntVars().contactArea_ << endl;

    subCInfo.evalVariables(
        cInfo.getcClass().getGeomModel().getCoM(),
        cInfo.gettClass().getGeomModel().getCoM(),
        cInfo.getcVars(),
        cInfo.gettVars()
    );

    // compute the normal force
    vector F = subCInfo.getFNe();
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact FNe " << F << endl;

    vector FNd = subCInfo.getFNd();
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact FNd " << FNd << endl;

    // clamp FNd if opposite direction to FNe
    if ((F & FNd) < 0 && mag(FNd) > mag(F))
    {
        FNd *= mag(F) / mag(FNd);
    }

    F += FNd;
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact FN " << F << endl;

    vector Ft = subCInfo.getFt(deltaT);
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact Ft " << Ft << endl;

    if (mag(Ft) > cInfo.getMu() * mag(F))
    {
        Ft *= cInfo.getMu() * mag(F) / mag(Ft);
    }
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact Ft clamped " << Ft << endl;
    F += Ft;

    vector FA = subCInfo.getFA();
    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact FA " << FA << endl;
    F -= FA;

    InfoH << parallelDEM_Info << "-- Particle-particle " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact F " << F << endl;

    InfoH << parallelDEM_Info << "-- Resolved Particle-particle contact: -- body "
            << subCInfo.getCPair().first() << " & -- body "
            << subCInfo.getCPair().second() << endl;

    // add the computed force to the affected bodies
    subCInfo.getOutForce().first().F = F;
    subCInfo.getOutForce().first().T = subCInfo.getcLVec() ^  F;
    subCInfo.getOutForce().second().F = -F;
    subCInfo.getOutForce().second().T = subCInfo.gettLVec() ^ -F;
    return true;
}
//---------------------------------------------------------------------------//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace contactModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
