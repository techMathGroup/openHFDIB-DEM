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
#include "prtContact.H"

#include "virtualMeshLevel.H"
#include "pspPairControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace contactModel
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//---------------------------------------------------------------------------//
void getContacts(
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
        prtcInfo.getContacts_Sphere();
    }
    else if (
        prtcInfo.getcClass().getGeomModel().getcType() == cluster
        ||
        prtcInfo.gettClass().getGeomModel().getcType() == cluster
        )
    {
        getContacts_Cluster(mesh, prtcInfo);
    }
    else
    {
        // if (Pstream::myProcNo() == 0)
        // {
            // scalar ranCellVol = mesh.V()[0];
            // prtcInfo.getContacts_ArbShape(ranCellVol);
        // }
        scalar ranCellVol = pow(virtualMeshLevel::getCharCellSize(),3);
        prtcInfo.getContacts_ArbShape(ranCellVol);        
    }

    prtcInfo.swapContactLists();
}
//---------------------------------------------------------------------------//
void getContacts_Cluster(
    const fvMesh&   mesh,
    prtContactInfo& prtcInfo
)
{
    prtcInfo.swapContactLists();
    std::vector<std::shared_ptr<geomModel>> cBodies;
    std::vector<std::shared_ptr<geomModel>> tBodies;

    bool cIsCluster = prtcInfo.getcClass().getGeomModel().isCluster();

    if(cIsCluster)
    {
        clusterBody& cCluster = dynamic_cast<clusterBody&>(prtcInfo.getcClass().getGeomModel());
        cBodies = cCluster.getClusterBodies();
    }
    else
    {
        cBodies.push_back(prtcInfo.getcClass().getGeomModelPtr());
    }

    if(prtcInfo.gettClass().getGeomModel().isCluster())
    {
        clusterBody& tCluster = dynamic_cast<clusterBody&>(prtcInfo.gettClass().getGeomModel());

        if (cIsCluster)
        {
            tBodies.push_back(tCluster.getClusterBodies()[0]);
        }
        else
        {
            tBodies = tCluster.getClusterBodies();
        }
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

            getContacts(mesh, tmpPrtCntInfo);

            if (tmpPrtCntInfo.getPrtSCList().size() > 0)
            {
                prtcInfo.getPrtSCList().insert(
                    prtcInfo.getPrtSCList().end(),
                    tmpPrtCntInfo.getPrtSCList().begin(),
                    tmpPrtCntInfo.getPrtSCList().end()
                );
            }
        }
    }
    prtcInfo.swapContactLists();
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
    std::shared_ptr<virtualMeshInfo>& vmInfo = subCInfo.getVMInfo();
    if (!vmInfo)
    {
        return false;
    }

    virtualMesh virtMesh(
        *vmInfo,
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

    bool isCCluster = cClass.getGeomModel().isCluster();

    if(isCCluster)
    {
        clusterBody& cCluster = dynamic_cast<clusterBody&>(cClass.getGeomModel());
        cBodies = cCluster.getClusterBodies();
    }
    else
    {
        cBodies.push_back(cClass.getGeomModelPtr());
    }

    if(tClass.getGeomModel().isCluster())
    {
        clusterBody& tCluster = dynamic_cast<clusterBody&>(tClass.getGeomModel());
        if (isCCluster)
        {
            tBodies.push_back(tCluster.getClusterBodies()[0]);
        }
        else
        {
            tBodies = tCluster.getClusterBodies();
        }
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
// synthesise the footprint of a contact patch as a polygon fan: a disc of
// radius r around the given centre (or in the pseudo-2D case a rectangle
// chord x emptyLength), in the plane perpendicular to the contact normal
// carried by vars
static void setContactFootprint
(
    prtContactVars& vars,
    const point& centre,
    const scalar r,
    const scalar emptyLength
)
{
    vars.contactPolygons_.clear();

    const vector n(vars.contactNormal_);
    if (mag(n) < SMALL || r <= 0) return;

    vector e1(1, 0, 0);
    if (mag(e1 & n) > 0.9) e1 = vector(0, 1, 0);
    vector u(n ^ e1);
    u /= mag(u);
    vector v(n ^ u);

    if (case3D)
    {
        const label nPts(24);
        pointField fan(nPts + 1);
        fan[0] = centre;
        for (label i = 0; i < nPts; i++)
        {
            const scalar a(2*Foam::constant::mathematical::pi*i/nPts);
            fan[i + 1] = centre + r*(cos(a)*u + sin(a)*v);
        }
        vars.contactPolygons_.append(fan);
    }
    else
    {
        vector eDir(vector::zero);
        eDir[emptyDim] = 1;
        vector t(n ^ eDir);
        t /= mag(t) + VSMALL;

        pointField fan(5);
        fan[0] = centre;
        fan[1] = centre + r*t + 0.5*emptyLength*eDir;
        fan[2] = centre - r*t + 0.5*emptyLength*eDir;
        fan[3] = centre - r*t - 0.5*emptyLength*eDir;
        fan[4] = centre + r*t - 0.5*emptyLength*eDir;
        vars.contactPolygons_.append(fan);
    }
}
//---------------------------------------------------------------------------//
void getPrtContactVars_ArbShape(
    const fvMesh&   mesh,
    ibContactClass& cClass,
    ibContactClass& tClass,
    prtSubContactInfo& subCInfo
)
{
    std::shared_ptr<virtualMeshInfo>& vmInfo = subCInfo.getVMInfo();
    if (!vmInfo)
    {
        return;
    }

    virtualMesh virtMesh(
        *vmInfo,
        cClass.getGeomModel(),
        tClass.getGeomModel()
    );

    scalar intersectedVolume(0);
    vector contactCenter(vector::zero);
    vector normalVector = vector::zero;
    scalar contactArea(0);

    // emptyScale restores the full extruded extent of the contact patch
    // along the (pseudo-2D) empty direction if the virtual mesh was
    // clipped there; it is 1 otherwise. The edge points span only the
    // clipped layer, so the swept-hull area is the 2D footprint and must
    // be scaled the same way as the volume (cf. the plane area formula in
    // wallContact.C).
    intersectedVolume = virtMesh.evaluateContact()*vmInfo->getEmptyScale();

    const label cId(subCInfo.getCPair().first());
    const label tId(subCInfo.getCPair().second());

    // psp automatic pass-through state of the pair (episode state is
    // tracked for every body pair; STL-wall pair contacts are excluded)
    const bool pspEligible
    (
        !cClass.getIsSTLWall() && !tClass.getIsSTLWall()
    );
    pspPairControl::forcedNormal pspRec;
    const bool hasPspRec
    (
        pspEligible && pspPairControl::getRecord(cId, tId, pspRec)
    );

    // a pair whose penetrator crossed the backside of the other body is
    // permanently force-free (the standard exit-face detection would
    // otherwise repel the penetrator out through the back face)
    if (hasPspRec && pspRec.passedThrough_)
    {
        subCInfo.getprtCntVars().contactPolygons_.clear();
        subCInfo.getprtCntVars().contactCenter_ = vector::zero;
        subCInfo.getprtCntVars().contactVolume_ = 0;
        subCInfo.getprtCntVars().contactNormal_ = vector::zero;
        subCInfo.getprtCntVars().contactArea_ = 0;
        subCInfo.getprtCntVars().synthesized_ = false;
        return;
    }

    const boundBox& cBBox(cClass.getGeomModel().getBounds());
    const boundBox& tBBox(tClass.getGeomModel().getBounds());
    const vector cCoM(cClass.getGeomModel().getCoM());
    const vector tCoM(tClass.getGeomModel().getCoM());

    // containment: one bounding box fully inside the other
    const bool cInT(tBBox.contains(cBBox));
    const bool tInC(cBBox.contains(tBBox));

    // psp: embedded pass-through state. When one body is fully inside
    // the other the contact hull degenerates - there is no rim between
    // the two surfaces (edge points <= 4) although the octree still
    // measures the intersection volume from the INSIDE&&INSIDE leaves.
    // The standard evaluation would drop the contact here. Keep it
    // alive instead and synthesise the contact variables: the normal is
    // the frozen entry normal, the volume is the measured intersection
    // volume grown linearly by the (signed) penetration displacement
    // since the embed onset times the maximum registered contact area
    // of the episode, and the area follows from the contained body's
    // extent along the penetration direction. The signed displacement
    // makes the embedded contact spring-like: the force grows while the
    // body buries itself and relaxes again on rebound.
    const bool pspEmbedded
    (
        pspEligible
        && intersectedVolume > 0
        && virtMesh.getEdgeSVPoints().size() <= 4
        && (cInT || tInC)
    );

    if (pspEmbedded)
    {
        // the contained body is the penetrator
        const label pen(cInT ? cId : tId);
        const vector& penCoM(cInT ? cCoM : tCoM);
        const vector& tarCoM(cInT ? tCoM : cCoM);
        const boundBox& contained(cInT ? cBBox : tBBox);

        pspPairControl::forcedNormal& pspR
        (
            pspPairControl::record(cId, tId)
        );

        if (!pspR.embedded_)
        {
            // embed onset (or re-embed after a rebound): freeze the
            // entry direction and the relative-CoM datum. An entry
            // direction registered by an application is kept.
            if (pspR.penetrator_ == -1 || !pspR.appControlled_)
            {
                if (mag(pspR.lastNormal_) > SMALL)
                {
                    // lastNormal_ points from the pair's second member
                    // towards its first member; e_ points from the
                    // penetrator towards the target
                    pspR.e_ = (pen == cId)
                      ? -pspR.lastNormal_
                      : pspR.lastNormal_;
                }
                else
                {
                    // contact first seen already embedded (the entry was
                    // swallowed by the detection gap): fall back to the
                    // CoM axis of the pair
                    pspR.e_ = tarCoM - penCoM;
                }
                pspR.e_ /= mag(pspR.e_) + VSMALL;
            }
            // psp: freeze the static footprint of the embedded contact -
            // the area of the last real contact patch (the max registered
            // geometric area, the same factor the overlap-volume
            // extrapolation uses), anchored at the patch centre of the
            // last geometric contact, i.e. the entry point on the
            // target's surface. The pressure carries the time dependence
            // (p = F/A_frozen). Without a preceding geometric contact
            // the anchor falls back to the CoM midpoint and the area to
            // the first synthesized area (assigned below)
            pspR.footprintArea_ = pspR.maxArea_;
            pspR.footprintCentre_ =
                pspR.maxArea_ > SMALL
              ? pspR.lastCentre_
              : 0.5*(cCoM + tCoM);

            // freeze the shape of the last real contact patch as well:
            // keep its polygon set (the logger projects these onto the
            // mesh) and translate it so that it sits at the frozen
            // anchor - for a footprint taken from the same sub-step the
            // shift is zero, it only moves a cached patch whose centre
            // differs from the anchor
            pspR.footprintPolygons_.clear();
            if (pspR.lastPolygons_.size() && pspR.maxArea_ > SMALL)
            {
                const vector fpShift
                (
                    pspR.footprintCentre_ - pspR.lastCentre_
                );

                pspR.footprintPolygons_.setSize(pspR.lastPolygons_.size());
                forAll(pspR.lastPolygons_, i)
                {
                    pointField& poly(pspR.footprintPolygons_[i]);
                    poly = pspR.lastPolygons_[i];
                    forAll(poly, j)
                    {
                        poly[j] += fpShift;
                    }
                }
            }

            pspR.penetrator_ = pen;
            pspR.relCoM0_ = penCoM - tarCoM;
            pspR.sMax_ = 0;
            pspR.embedded_ = true;
            pspR.episodeActive_ = true;

            InfoH << parallelDEM_Info << "pspHFDIBDEM: body " << pen
                    << " fully immersed in body "
                    << (pen == cId ? tId : cId)
                    << " at time " << mesh.time().value()
                    << " (geometric contact ended: full immersion, "
                       "embedded contact continues)" << endl;
        }

        // signed penetration since the embed onset: positive while the
        // penetrator keeps burying itself, negative on rebound
        const vector pspPenRelCoM
        (
            pspR.penetrator_ == cId ? cCoM - tCoM : tCoM - cCoM
        );
        const scalar pspS((pspPenRelCoM - pspR.relCoM0_) & pspR.e_);
        pspR.sMax_ = max(pspR.sMax_, pspS);

        // extrapolated overlap volume: measured contained volume
        // increased linearly by (displacement x max registered contact
        // area during the contact)
        const scalar pspVC
        (
            max(intersectedVolume + pspS*pspR.maxArea_, 0)
        );

        // cross-section estimate: overlap volume over the extent
        // of the contained body along the penetration direction
        const scalar pspProj
        (
            mag(pspR.e_ & contained.span())
        );
        const scalar pspArea
        (
            pspProj > SMALL
          ? pspVC/(pspProj*vmInfo->getEmptyScale())
          : pspVC/pow(vmInfo->subVolumeV, 1.0/3.0)
        );

        // psp: fallback of the frozen static footprint for a contact
        // first seen already embedded (no geometric patch was registered
        // in this episode): use the first synthesised area
        if (pspR.footprintArea_ <= SMALL)
        {
            pspR.footprintArea_ = pspArea;
        }

        subCInfo.getprtCntVars().contactCenter_ = 0.5*(cCoM + tCoM);
        subCInfo.getprtCntVars().contactVolume_ = pspVC;
        // the normal force on the penetrator acts along -e_ (the entry
        // direction points into the target) and always decelerates the
        // relative approach, whichever pair member is c
        subCInfo.getprtCntVars().contactNormal_ =
            (pspR.penetrator_ == cId) ? -pspR.e_ : pspR.e_;
        // psp: static footprint - the reported contact area is the frozen
        // footprint area, so the derived pressure p = F/A_frozen carries
        // the whole time dependence of the embedded load
        subCInfo.getprtCntVars().contactArea_ = pspR.footprintArea_;

        // psp: apply the static footprint of the embedded contact. Its
        // polygons were frozen at the embed onset from the last real
        // contact patch (same shape, orientation and anchor as the
        // geometric contact, so the projected cell set does not change).
        // Only without a cached patch (contact first seen already
        // embedded) synthesize an equivalent-area disc (3D) / chord
        // rectangle (2D) at the frozen anchor. The application
        // contact-pressure logger projects these polygons onto the mesh,
        // so without them the pCont field would stay empty while a body
        // is submerged
        if (pspR.footprintPolygons_.size())
        {
            subCInfo.getprtCntVars().contactPolygons_
                = pspR.footprintPolygons_;
        }
        else
        {
            scalar r(0);
            scalar emptyLength(0);
            if (case3D)
            {
                r = sqrt
                (
                    pspR.footprintArea_
                   /Foam::constant::mathematical::pi
                );
            }
            else
            {
                emptyLength = contained.span()[emptyDim];
                r = 0.5*pspR.footprintArea_/max(emptyLength, VSMALL);
            }

            setContactFootprint
            (
                subCInfo.getprtCntVars(),
                pspR.footprintCentre_,
                r,
                emptyLength
            );
        }

        subCInfo.getprtCntVars().synthesized_ = true;

        InfoH << parallelDEM_Info << "-- Detected Particle-particle contact [embedded]: -- body "
                << cId << " & -- body "
                << tId << endl;
        InfoH << parallelDEM_Info << "-- Particle-particle " << cId
                << "-" << tId << " contact volume "
                << pspVC << endl;
        InfoH << parallelDEM_Info << "-- Particle-particle " << cId
                << "-" << tId << " contact area "
                << pspArea << endl;
        InfoH << parallelDEM_Info << "-- Particle-particle " << cId
                << "-" << tId << " static footprint area "
                << pspR.footprintArea_ << endl;
        return;
    }

    // psp: the embedded phase just ended (the penetrator is no longer
    // fully contained). Distinguish how: still at (or beyond) the
    // deepest registered point means the leading edge crossed the
    // backside - switch the contact force off permanently, unless an
    // application controls the pair and runs its own pass-through
    // trigger; a clearly shallower depth means a rebound - the standard
    // contact takes over again
    if (hasPspRec && pspRec.embedded_)
    {
        pspPairControl::forcedNormal& pspR
        (
            pspPairControl::record(cId, tId)
        );

        const vector pspPenRelCoM
        (
            pspR.penetrator_ == cId ? cCoM - tCoM : tCoM - cCoM
        );
        const scalar pspS((pspPenRelCoM - pspR.relCoM0_) & pspR.e_);

        if (pspS >= pspR.sMax_ - SMALL && !pspR.appControlled_)
        {
            // backside crossed with the movement direction unchanged
            pspR.passedThrough_ = true;
            pspR.embedded_ = false;

            InfoH << parallelDEM_Info << "pspHFDIBDEM: body "
                    << pspR.penetrator_ << " penetrated the backside of body "
                    << (pspR.penetrator_ == cId ? tId : cId)
                    << " at time " << mesh.time().value()
                    << " (contact force off)" << endl;

            subCInfo.getprtCntVars().contactPolygons_.clear();
            subCInfo.getprtCntVars().contactCenter_ = vector::zero;
            subCInfo.getprtCntVars().contactVolume_ = 0;
            subCInfo.getprtCntVars().contactNormal_ = vector::zero;
            subCInfo.getprtCntVars().contactArea_ = 0;
            subCInfo.getprtCntVars().synthesized_ = false;
            return;
        }

        if (pspS < pspR.sMax_ - SMALL)
        {
            InfoH << parallelDEM_Info << "pspHFDIBDEM: body "
                    << pspR.penetrator_ << " rebounded inside body "
                    << (pspR.penetrator_ == cId ? tId : cId)
                    << " at time " << mesh.time().value()
                    << " (embedded contact ended by rebound: standard "
                       "contact resumed)" << endl;
        }
        pspR.embedded_ = false;
    }

    if(virtMesh.getEdgeSVPoints().size() <= 4)
    {
        intersectedVolume = 0;
    }

    if (intersectedVolume > 0)
    {
        Tuple2<scalar,vector> surfaceAndNormal =
            virtMesh.get3DcontactNormalAndSurface
            (
                cClass.getGeomModel().getcType() == nonConvex
                ||
                tClass.getGeomModel().getcType() == nonConvex
            );
        contactArea = surfaceAndNormal.first()*vmInfo->getEmptyScale();
        normalVector = surfaceAndNormal.second();
        contactCenter = virtMesh.getContactCenter();

        // hull area degenerated (collinear edge points, flat-face
        // touch): estimate the area from the contact shape the way
        // the wall side does - n leaves of volume svVol carry
        // roughly n*svVol^(2/3), with V = n*svVol this is V^(2/3)
        // scaled by the leaf volume
        const scalar relASMALL
        (
            SMALL*pow(vmInfo->subVolumeV, 2.0/3.0)
        );

        if (contactArea < relASMALL)
        {
            contactArea = intersectedVolume
                /pow(vmInfo->subVolumeV, 1.0/3.0)
               *vmInfo->getEmptyScale();
        }
    }


    if (intersectedVolume > 0 && contactArea > 0)
    {
        subCInfo.getprtCntVars().contactCenter_ = contactCenter;
        subCInfo.getprtCntVars().contactVolume_ = intersectedVolume;
        subCInfo.getprtCntVars().contactNormal_ = normalVector;
        subCInfo.getprtCntVars().contactArea_ = contactArea;
        subCInfo.getprtCntVars().contactPolygons_ = virtMesh.getContactPolygons();
        subCInfo.getprtCntVars().synthesized_ = false;

        // psp: track the contact episode of the pair for the automatic
        // pass-through support - the maximum registered contact area (the
        // overlap-volume extrapolation factor of the embedded contact)
        // and the last standard contact normal (the entry-direction
        // fallback at the embed onset)
        if (pspEligible)
        {
            pspPairControl::forcedNormal& pspR
            (
                pspPairControl::record(cId, tId)
            );

            if (!pspR.episodeActive_)
            {
                pspR.episodeActive_ = true;
                pspR.maxArea_ = 0;
            }

            if (!pspR.embedded_)
            {
                pspR.maxArea_ = max(pspR.maxArea_, contactArea);
                pspR.lastNormal_ = normalVector;
                pspR.lastCentre_ = contactCenter;
                // psp: remember the shape of the last real contact patch
                // (the static footprint of a following embed)
                pspR.lastPolygons_
                    = subCInfo.getprtCntVars().contactPolygons_;
            }
        }
    }
    else
    {
        subCInfo.getprtCntVars().contactPolygons_.clear();
        subCInfo.getprtCntVars().contactCenter_ = vector::zero;
        subCInfo.getprtCntVars().contactVolume_ = 0;
        subCInfo.getprtCntVars().contactNormal_ = vector::zero;
        subCInfo.getprtCntVars().contactArea_ = 0;
        subCInfo.getprtCntVars().synthesized_ = false;
    }
}
//---------------------------------------------------------------------------//
// sphere contacts carry no hull: synthesise the footprint of the contact
// patch around the contact centre
static void setSphereFootprint
(
    prtContactVars& vars,
    const scalar r,
    const scalar emptyLength
)
{
    setContactFootprint(vars, vars.contactCenter_, r, emptyLength);
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

    subCInfo.getprtCntVars().contactPolygons_.clear();

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
        setSphereFootprint(subCInfo.getprtCntVars(), tRadius, 0);
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
        setSphereFootprint
        (
            subCInfo.getprtCntVars(),
            sqrt(sqr(tRadius) - sqr(xLength)),
            0
        );
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
        setSphereFootprint
        (
            subCInfo.getprtCntVars(),
            sqrt(sqr(tRadius) - sqr(xLength)),
            emptyLength
        );
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

    bool isCCluster = cClass.getGeomModel().isCluster();

    if(isCCluster)
    {
        clusterBody& cCluster = dynamic_cast<clusterBody&>(cClass.getGeomModel());
        cBodies = cCluster.getClusterBodies();
    }
    else
    {
        cBodies.push_back(cClass.getGeomModelPtr());
    }

    if(tClass.getGeomModel().isCluster())
    {
        clusterBody& tCluster = dynamic_cast<clusterBody&>(tClass.getGeomModel());
        if (isCCluster)
        {
            tBodies.push_back(tCluster.getClusterBodies()[0]);
        }
        else
        {
            tBodies = tCluster.getClusterBodies();
        }
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

            if (subCInfo.getVMInfo())
            {
                tmpSubCInfoI.setVMInfo(*(subCInfo.getVMInfo()));
            }

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

    // pair contacts of stlBased collision patch bodies are reported
    // as particle-wall, not particle-particle
    const word cntType
    (
        cInfo.getcClass().getIsSTLWall() || cInfo.gettClass().getIsSTLWall()
        ? "Particle-wall"
        : "Particle-particle"
    );

    InfoH << parallelDEM_Info << "-- Detected " << cntType << " contact: -- body "
            << subCInfo.getCPair().first() << " & -- body "
            << subCInfo.getCPair().second() << endl;
    InfoH << parallelDEM_Info << "-- " << cntType << " " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact cBody pos: "
            << cInfo.getcClass().getGeomModel().getCoM() << " & tBody pos: "
            << cInfo.gettClass().getGeomModel().getCoM() << endl;
    InfoH << parallelDEM_Info << "-- body "<< subCInfo.getCPair().first() <<"  linear velocity:"
        << cInfo.getcVars().Vel_ << " magnitude: " << mag(cInfo.getcVars().Vel_) <<endl;
    InfoH << parallelDEM_Info << "-- body "<< subCInfo.getCPair().first() <<"  angular velocity:"
        << cInfo.getcVars().omega_ << " magnitude: " << mag(cInfo.getcVars().omega_) <<endl;

    InfoH << parallelDEM_Info << "-- body "<< subCInfo.getCPair().second() <<"  linear velocity:"
        << cInfo.gettVars().Vel_ << " magnitude: " << mag(cInfo.gettVars().Vel_) <<endl;
    InfoH << parallelDEM_Info << "-- body "<< subCInfo.getCPair().second() <<"  angular velocity:"
        << cInfo.gettVars().omega_ << " magnitude: " << mag(cInfo.gettVars().omega_) <<endl;
    InfoH << parallelDEM_Info << "-- " << cntType << " " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact center "
            << subCInfo.getprtCntVars().contactCenter_ << endl;
    InfoH << parallelDEM_Info << "-- " << cntType << " " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact normal "
            << subCInfo.getprtCntVars().contactNormal_ << endl;
    InfoH << parallelDEM_Info << "-- " << cntType << " " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact volume "
            << subCInfo.getprtCntVars().contactVolume_ << endl;
    InfoH << parallelDEM_Info << "-- " << cntType << " " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact area "
            << subCInfo.getprtCntVars().contactArea_ << endl;

    subCInfo.evalVariables(
        cInfo.getcClass(),
        cInfo.gettClass(),
        cInfo.getcVars(),
        cInfo.gettVars()
    );

    // compute the normal force
    vector F = subCInfo.getFNe();
    InfoH << parallelDEM_Info << "-- " << cntType << " " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact FNe " << F << endl;

    vector FNd = subCInfo.getFNd();
    InfoH << parallelDEM_Info << "-- " << cntType << " " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact FNd " << FNd << endl;

    // clamp FNd if opposite direction to FNe
    if ((F & FNd) < 0 && mag(FNd) > mag(F))
    {
        FNd *= mag(F) / mag(FNd);
    }

    F += FNd;
    InfoH << parallelDEM_Info << "-- " << cntType << " " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact FN " << F << endl;

    // the coulomb ceiling caps the applied force and the stored
    // tangential state inside getFt
    vector Ft = subCInfo.getFt(deltaT, cInfo.getMu()*mag(F));
    InfoH << parallelDEM_Info << "-- " << cntType << " " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact Ft " << Ft << endl;

    F += Ft;

    vector FA = subCInfo.getFA();
    InfoH << parallelDEM_Info << "-- " << cntType << " " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact FA " << FA << endl;
    F -= FA;

    // add the computed force to the affected bodies
    subCInfo.getOutForce().first().F = F;
    subCInfo.getOutForce().first().T = subCInfo.getcLVec() ^  F;
    subCInfo.getOutForce().second().F = -F;
    subCInfo.getOutForce().second().T = subCInfo.gettLVec() ^ -F;

    InfoH << parallelDEM_Info << "-- " << cntType << " " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact F " << F << endl;

    InfoH << parallelDEM_Info << "-- " << cntType << " " <<subCInfo.getCPair().first() <<"-"<<subCInfo.getCPair().second() << " contact T " << subCInfo.getOutForce().first().T << endl;

    InfoH << parallelDEM_Info << "-- Resolved " << cntType << " contact: -- body "
            << subCInfo.getCPair().first() << " & -- body "
            << subCInfo.getCPair().second() << endl;

    return true;
}
//---------------------------------------------------------------------------//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace contactModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
