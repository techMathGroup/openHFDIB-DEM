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
    by the Free Software Foundation, either version 3 of the License, or (at
    your option) any later version.

    openHFDIB-DEM is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with openHFDIB-DEM.  If not, see <http://www.gnu.org/licenses/>.

Application
    vmWallContactBenchmark

Description
    Benchmark of the wall-contact geometry accuracy: an icosphere
    (STL) pressed against a plane wall, evaluated through the
    production wall-contact chain (wallContactInfo ->
    findContactAreas -> getWallContactVars -> ArbShape VM path).

    References: the smooth sphere cap

        V = pi * xH^2 * (3*R - xH) / 3,  xH = R - d
        A = pi * (R^2 - d^2)

    AND the exact faceted cap of the STL polyhedron (clip every
    triangle to y <= 0, divergence-theorem volume; convex section
    polygon area). The faceted reference is the primary one - the
    virtual mesh resolves the STL body, not the smooth sphere,
    and at icosphere faceting the two differ by percent-level
    offsets that would otherwise mask the svEdge convergence.

    For each (d, level) the benchmark reports the measured
    contact volume and wetted area of both arms (exactSubVolume
    true and false - the identical-binary A/B pair) with errors
    against both references. The gate is the convergence ORDER
    only (evaluator script): the exact arm must halve-plus its
    error per level refinement (second order), the legacy arm is
    expected first order by design. The tilted-wall leg is
    reported but not gated (known open issue).

    Each evaluation additionally prints the split of the
    measured contact volume into the vm-flood term and the
    internal-box term (both re-derived per the production
    loops of getWallContactVars_ArbShape) on a line of its own:
    the tilted-wall bias is level-independent, pointing at the
    charCellSize-scale scaffolding, and this split separates
    the two candidate sources. The evaluator script ignores
    lines it does not parse, so the plotting output is
    unchanged.

    The ccs-sweep leg repeats the tilted wall at
    charCellSize/2 with the levels shifted one down: both
    legs then share the svEdge range (matched svEdge), and the
    error comparison at matched svEdge isolates the pure
    charCellSize dependence of the tilted-wall error from the
    svEdge-order virtual-mesh machinery (a first-order-in-ccs
    floor halves when ccs halves; an svEdge-order term stays).

    The exact arm additionally prints a missed-element scan
    (first level of each leg only, the SM keep-set is
    level-independent): the SM lattice is re-derived locally
    (the floor/ceil bbox snap and ceil matrix sizing of
    constructSM) and every element with behind-wall body
    material is matched against the sub-contact boxes. This
    attributes the residual tilted-wall volume error to either
    the detection gate (elements the possibleSMContact flood
    never visited, missedBody > 0) or the SM-scale box
    quantization (missedBody == 0 with a flat idealBodyVol
    error) - the two candidate mechanisms of the phase-2
    decision.

    Run from the case directory of this benchmark (needs
    constant/triSurface/sphere.stl and a blockMesh-generated mesh).

Contributors
    Martin Isoz (2019-*)
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallContactInfo.H"
#include "wallContact.H"
#include "nonConvexBody.H"
#include "virtualMeshLevel.H"
#include "materialProperties.H"
#include "wallMatInfo.H"
#include "wallPlaneInfo.H"
#include "interAdhesion.H"
#include "outputHFDIBDEM.H"
#include "wallSubContactInfo.H"
#include "virtualMeshWall.H"
#include "triSurface.H"
#include <memory>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

//- re-derive the SM lattice exactly as wallContactInfo::
// constructSM does (floor/ceil bbox snap at ccs, ceil matrix
// sizing, vertex centers at bBox.min + i*ccs) and scan every
// element with behind-wall body material against the sub-contact
// boxes. the raw element boxes of the exact arm (post a1b563e)
// are min-corner-identified, so covered elements are matched
// exactly. prints the missed-element attribution the phase-2
// decision needs: idealBodyVol (the behind-wall body material
// of the SM slab, wall-clipped), coveredVol (the part the
// contact carries), missedVol, and the element-count split
void missedElementScan
(
    geomModel& geom,
    const List<planePolyClip::halfSpace>& wallPlanes,
    wallSubContactInfo& sC,
    const scalar charCellSize,
    const vector& wallN,
    const point& wallP
)
{
    // the SM lattice of constructSM: body bbox snapped to the
    // ccs grid (floor at min, ceil at max), ceil matrix sizing
    const boundBox bodyBB(geom.getBounds());
    point minP(bodyBB.min());
    point maxP(bodyBB.max());
    for (int i = 0; i < 3; i++)
    {
        minP[i] = std::floor(minP[i]/charCellSize)*charCellSize;
        maxP[i] = std::ceil(maxP[i]/charCellSize)*charCellSize;
    }
    const vector nCells
    (
        std::ceil((maxP.x() - minP.x())/charCellSize),
        std::ceil((maxP.y() - minP.y())/charCellSize),
        std::ceil((maxP.z() - minP.z())/charCellSize)
    );

    // the element boxes the contact actually carries: the raw
    // SM element boxes of the export and internal lists
    List<boundBox> carried;
    forAll(sC.getInternalElements(), iI)
    {
        carried.append(sC.getInternalElements()[iI].second());
    }
    for (label vmI = 0; vmI < sC.getVMContactSize(); vmI++)
    {
        autoPtr<virtualMeshWallInfo>& vmWInfo
        (
            sC.getVMContactInfo(vmI)
        );
        if (vmWInfo.valid())
        {
            carried.append(boundBox(vmWInfo->getbBox()));
        }
    }

    scalar idealBodyVol(0);
    scalar coveredVol(0);
    label nMissed(0);
    label nCovered(0);

    // the scan: every element with a corner strictly behind the
    // wall and body material (a vertex or the centroid inside,
    // or a facet crossing the box - the count-all rule of the
    // VM leaf classification)
    for (label i = 0; i < nCells.x(); i++)
    {
        for (label j = 0; j < nCells.y(); j++)
        {
            for (label k = 0; k < nCells.z(); k++)
            {
                const point eMin
                (
                    minP
                  + vector(i, j, k)*charCellSize
                );
                const boundBox element
                (
                    eMin,
                    eMin + vector::one*charCellSize
                );

                // some corner strictly behind the wall (the
                // keep-set of boxPenetratesWall)
                bool penetrates(false);
                {
                    const pointField cPts(element.points());
                    forAll(cPts, cP)
                    {
                        if (((cPts[cP] - wallP) & wallN) > 0)
                        {
                            penetrates = true;
                            break;
                        }
                    }
                }
                if (!penetrates)
                {
                    continue;
                }

                // body material in the element: centroid or
                // vertex inside, or a facet crossing the box
                bool hasBody(geom.pointInside(element.midpoint()));
                if (!hasBody)
                {
                    const pointField ePts(element.points());
                    forAll(ePts, eP)
                    {
                        if (geom.pointInside(ePts[eP]))
                        {
                            hasBody = true;
                            break;
                        }
                    }
                }
                if (!hasBody)
                {
                    hasBody = geom.boxIntersectsSurface(element);
                }
                if (!hasBody)
                {
                    continue;
                }

                // the ideal behind-wall body volume of this
                // element: the box clipped by the wall planes.
                // the body part inside is bounded by this clip
                // from above, so idealBodyVol is an upper bound
                // of the true behind-wall body material
                scalar eV(0);
                vector eC(vector::zero);
                planePolyClip::volumeAndCentroid
                (
                    element,
                    wallPlanes,
                    eV,
                    eC
                );
                idealBodyVol += eV;

                // covered: the min corner falls in some carried
                // box (the raw element boxes of the exact arm
                // are the element boxes themselves)
                bool isCovered(false);
                forAll(carried, cI)
                {
                    if
                    (
                        mag(carried[cI].min() - eMin)
                      < 0.5*SMALL*charCellSize
                    )
                    {
                        isCovered = true;
                        break;
                    }
                }
                if (isCovered)
                {
                    nCovered++;
                    coveredVol += eV;
                }
                else
                {
                    nMissed++;
                }
            }
        }
    }

    Info << "    missed scan: idealBodyVol " << idealBodyVol
         << " coveredVol " << coveredVol
         << " missedVol " << idealBodyVol - coveredVol
         << " covered/missed elements " << nCovered
         << "/" << nMissed << endl;
}

//- exact faceted references of the STL cap below the wall plane:
// the body is the STL polyhedron shifted to center (0, d, 0); the
// wall is y = 0. volume by clipping every triangle to y <= 0 and
// integrating with the divergence theorem (clipped faces keep
// their winding), area of the convex section polygon by sorting
// its vertices around the centroid
void facetedReference
(
    const triSurface& surf,
    const scalar& d,
    scalar& capV,
    scalar& capA
)
{
    capV = 0;
    capA = 0;

    DynamicList<point> sectionPts;

    forAll(surf, fI)
    {
        const triFace& f(surf[fI]);
        pointField pts(3);
        for (label i = 0; i < 3; i++)
        {
            const point& p(surf.localPoints()[f[i]]);
            pts[i] = point(p.x(), p.y() + d, p.z());
        }

        // sutherland-hodgman: keep y <= 0
        DynamicList<point> poly;
        for (label i = 0; i < 3; i++)
        {
            const point& p1(pts[i]);
            const point& p2(pts[(i + 1) % 3]);
            const bool in1(p1.y() <= 0);
            const bool in2(p2.y() <= 0);
            if (in1)
            {
                poly.append(p1);
            }
            if (in1 != in2)
            {
                const scalar t((0 - p1.y())/(p2.y() - p1.y()));
                poly.append(p1 + t*(p2 - p1));
            }
        }

        // fan triangulation from vertex 0 keeps the winding
        for (label k = 1; k < poly.size() - 1; k++)
        {
            const point& a(poly[0]);
            const point& b(poly[k]);
            const point& c(poly[k + 1]);
            capV += ((a ^ b) & c)/6.0;
        }

        // section segments: edges crossing y = 0
        for (label i = 0; i < 3; i++)
        {
            const point& p1(pts[i]);
            const point& p2(pts[(i + 1) % 3]);
            if ((p1.y() > 0) != (p2.y() > 0))
            {
                const scalar t((0 - p1.y())/(p2.y() - p1.y()));
                const point x(p1 + t*(p2 - p1));
                bool dup(false);
                forAll(sectionPts, sI)
                {
                    if (mag(sectionPts[sI] - x) < 1e-12)
                    {
                        dup = true;
                        break;
                    }
                }
                if (!dup)
                {
                    sectionPts.append(x);
                }
            }
        }
    }

    // convex section polygon: sort around the centroid by angle
    if (sectionPts.size() >= 3)
    {
        point c(vector::zero);
        forAll(sectionPts, sI)
        {
            c += sectionPts[sI];
        }
        c /= sectionPts.size();

        // selection sort by atan2 in the (x, z) plane
        DynamicList<point> sorted;
        sorted.transfer(sectionPts);
        for (label i = 0; i < sorted.size() - 1; i++)
        {
            label best(i);
            scalar bestA(std::atan2(sorted[i].z() - c.z(), sorted[i].x() - c.x()));
            for (label j = i + 1; j < sorted.size(); j++)
            {
                const scalar a(std::atan2(sorted[j].z() - c.z(), sorted[j].x() - c.x()));
                if (a < bestA)
                {
                    best = j;
                    bestA = a;
                }
            }
            Swap(sorted[i], sorted[best]);
        }

        // shoelace on (x, z)
        for (label i = 0; i < sorted.size(); i++)
        {
            const point& p1(sorted[i]);
            const point& p2(sorted[(i + 1) % sorted.size()]);
            capA += p1.x()*p2.z() - p2.x()*p1.z();
        }
        capA = 0.5*mag(capA);
    }
}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "convergence study of the exact wall-contact geometry "
        "(sphere cap) against the faceted and analytic references"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // silence the production-chain diagnostics
    outputHFDIBDEM::InfoH.setOutput(false, false, false, false, false);

    // --------------------------------------------------------------
    // case parameters of the study
    // --------------------------------------------------------------

    // analytic radius of the test sphere (the STL was generated
    // with exactly this radius)
    const scalar R(1.0);

    // center-to-wall distance d of the sweep: near touching to
    // half immersed in the wall - the deep end keeps the cap
    // well-resolved at coarse levels
    const FixedList<scalar,3> dList({0.95*R, 0.75*R, 0.5*R});

    // virtualMesh levels of the sweep: svEdge halves each step
    const FixedList<label,3> levelListBase({3, 4, 5});

    const scalar charCellSize(1.0/15.0);
    const scalar maxSubVolumes(1000000000);

    // virtualMesh levels of the ccs-sweep leg: shifted one down
    // from levelList so the svEdge range matches the 1/15 legs
    // (svEdge = ccs/2^level; halving ccs and lowering the level
    // by one keeps svEdge) - the matched-svEdge comparison of
    // the sweep leg and the tilted leg isolates the pure
    // charCellSize dependence of the tilted-wall error
    const FixedList<label,3> sweepLevelList({2, 3, 4});

    // wall: y = 0, fluid above. the stored normal points OUT of
    // the fluid, into the wall (the convention of every plane
    // entry of the fallingSphereOnWall case); "inside the plane"
    // ((x - p) & n < 0) is then the fluid side
    const word wallName("bot");
    const vector wallN(0, -1, 0);
    const point wallP(0, 0, 0);

    // --------------------------------------------------------------
    // static setup normally done by openHFDIBDEM from HFDIBDEMDict
    // --------------------------------------------------------------

    wallPlaneInfo::wallPlaneInfo_insert(wallName, wallN, wallP);

    materialProperties::matProps_insert
    (
        "particle",
        materialInfo("particle", 1e8, 0.3, 0.0, 0.0, 1.0)
    );
    wallMatInfo::wallMatInfo_insert
    (
        wallName,
        materialProperties::getMatProps()["particle"]
    );

    virtualMeshLevel::setMaxSubVolumes(maxSubVolumes);

    // body: the icosphere STL, centered at the origin of the STL
    // coordinates; positioned by shifting the surface points
    autoPtr<nonConvexBody> sphereBody
    (
        new nonConvexBody
        (
            mesh,
            word("constant/triSurface/sphere.stl")
        )
    );

    // raw surface for the faceted references
    const triSurface refSurf("constant/triSurface/sphere.stl");

    // --------------------------------------------------------------
    // per (d, level) evaluation
    // --------------------------------------------------------------

    label nFailed(0);

    Info << "vmWallContactBenchmark: sphere-cap convergence study"
         << endl;

    // wall-orientation sweep: leg 0 is the lattice-aligned wall,
    // leg 1 the 30-degree tilt, leg 2 the off-origin aligned wall,
    // leg 3 the ccs-sweep (30-degree tilt at ccs/2, levels shifted
    // one down - matched svEdge) (same cap references - the
    // analytic and faceted cap of a sphere cut by a plane are
    // invariant to a shared rotation of sphere and plane)
    const FixedList<vector,4> wallOrientations
    ({
        vector(0, -1, 0),
        vector(0.5, -Foam::sqrt(3.0)/2.0, 0),
        vector(0, -1, 0),
        vector(0.5, -Foam::sqrt(3.0)/2.0, 0)
    });

    // plane points of the legs: the tilted legs rotate about the
    // origin, so they share wallP; the off-origin leg sits at
    // y = 1/30 (inside the mesh, below every sphere placement,
    // noisy for every SM grid line of this ccs)
    const FixedList<point,4> wallPoints
    ({
        wallP,
        wallP,
        point(0, 1.0/30.0, 0),
        wallP
    });

    // charCellSize of the legs: every leg but the sweep runs at
    // the base ccs; the sweep halves it (levels compensate, so
    // the svEdge range matches)
    const FixedList<scalar,4> legCcs
    ({
        charCellSize,
        charCellSize,
        charCellSize,
        0.5*charCellSize
    });

    for (label wI = 0; wI < wallOrientations.size(); wI++)
    {
    // the level list of this leg (matched svEdge for the sweep)
    const FixedList<label,3>& levelList
    (
        wI == 3 ? sweepLevelList : levelListBase
    );

    // the ccs of this leg
    const scalar ccs(legCcs[wI]);

    if (wI > 0)
    {
        // re-set the wall registry entry for this leg (leg 0
        // entered it at the origin through wallPlaneInfo_insert;
        // the HashTable insert would not overwrite).
        // leg 1 rotates the wall and the sphere placement by the
        // same 30-degree tilt about z: the tilted leg sees
        // diagonal wall half-spaces in every band leaf.
        // leg 2 shifts the aligned wall along its normal to
        // y = 1/30: an off-origin flat wall, the projected plane
        // box carries 1-ULP noise in its normal span.
        // leg 3 repeats the tilted wall at ccs/2 with the levels
        // shifted one down (matched svEdge - the ccs-order probe)
        wallPlaneInfo::wallPlaneInfo_.set
        (
            wallName,
            List<vector>{{wallOrientations[wI], wallPoints[wI]}}
        );

        Info << nl
             << "=== "
             << (wI == 1
                 ? "tilted-wall leg (30 degrees about z)"
                 : wI == 2
                 ? "off-origin aligned-wall leg (y = 1/30)"
                 : "ccs-sweep tilted-wall leg (30 degrees, "
                   "ccs/2, matched svEdge)")
             << " ===" << endl;
    }

    forAll(dList, dI)
    {
        const scalar d(dList[dI]);
        // the sweep leg keeps the sphere placement of the tilted
        // leg: the wall and sphere rotate together, the ccs is
        // the only change (the level shift compensates svEdge)
        const scalar xH(R - d);
        const scalar VrefAna
        (
            constant::mathematical::pi*xH*xH*(3*R - xH)/3
        );
        const scalar ArefAna
        (
            constant::mathematical::pi*(R*R - d*d)
        );

        scalar VrefF(0);
        scalar ArefF(0);
        facetedReference(refSurf, d, VrefF, ArefF);

        Info << nl << "--- center distance d = " << d
             << " (d/R = " << d/R << ")" << endl;
        Info << "    analytic cap volume " << VrefAna
             << ", faceted cap volume " << VrefF
             << " (faceting " << VrefF/VrefAna - 1 << ")" << endl;
        Info << "    analytic disk area  " << ArefAna
             << ", faceted section area " << ArefF
             << " (faceting " << ArefF/ArefAna - 1 << ")" << endl;

        // place the sphere: STL center to (0, d, 0) on the aligned
        // leg; the tilted legs rotate sphere and wall together
        // (rotation about z through the origin), so the cap
        // geometry relative to the plane is identical; the
        // off-origin leg shifts sphere and wall together, so the
        // cap geometry relative to the plane is identical too
        {
            pointField pos(sphereBody->getBodyPoints());
            const vector shift
            (
                vector(0, d, 0)
              + vector(0, 1.0/30.0, 0)*(wI == 2 ? 1 : 0)
              - sphereBody->getCoM()
            );
            forAll(pos, pI)
            {
                pos[pI] += shift;
            }
            if (wI == 1 || wI == 3)
            {
                const scalar a(30.0*constant::mathematical::pi/180.0);
                const tensor rot
                (
                    Foam::cos(a), -Foam::sin(a), 0,
                    Foam::sin(a),  Foam::cos(a), 0,
                    0, 0, 1
                );
                forAll(pos, pI)
                {
                    pos[pI] = rot & pos[pI];
                }
            }
            sphereBody->setBodyPosition(pos);
        }

        forAll(levelList, lI)
        {
            const label level(levelList[lI]);
            virtualMeshLevel::setVirtualMeshLevel(level, ccs);
            const scalar svEdge
            (
                ccs/virtualMeshLevel::getLevelOfDivision()
            );

            for (label exactArm = 0; exactArm <= 1; exactArm++)
            {
                const bool exact(exactArm == 1);
                virtualMeshLevel::setExactSubVolume(exact);

                // contact variables of the body (velocities are
                // irrelevant to the geometric evaluation)
                vector Vel(vector::zero);
                scalar omega(0);
                vector Axis(0, 0, 1);
                scalar M0(1);
                scalar M(1);
                dimensionedScalar rhoS("rhoS", dimDensity, 4000);

                ibContactVars cVars
                (
                    0, Vel, omega, Axis, M0, M, rhoS
                );

                std::shared_ptr<geomModel> gm
                (
                    sphereBody.get(),
                    [](geomModel*){}                           // no ownership
                );
                ibContactClass ibClass
                (
                    gm,
                    "particle"
                );

                wallContactInfo wallCntInfo(ibClass, cVars);

                wallCntInfo.detectWallContact();

                if (wallCntInfo.getContactPatches().size() == 0)
                {
                    Info << "  level " << level << " "
                         << (exact ? "exact" : "legacy")
                         << ": [FAIL] no wall contact detected"
                         << endl;
                    nFailed++;
                    continue;
                }

                wallCntInfo.findContactAreas();

                DynamicList<wallSubContactInfo*> subContacts;
                wallCntInfo.registerSubContactList(subContacts);

                if (subContacts.size() == 0)
                {
                    Info << "  level " << level << " "
                         << (exact ? "exact" : "legacy")
                         << ": [FAIL] no sub-contact found"
                         << endl;
                    nFailed++;
                    continue;
                }

                // volume split BEFORE the production call: the
                // production loop of getWallContactVars_ArbShape
                // copies the info autoPtrs by value
                // (autoPtr<virtualMeshWallInfo> vmWInfo = ...)
                // and the v2412 autoPtr copy ctor is a move in
                // disguise, so after the production call the
                // list entries are empty. the split uses a
                // reference and never steals; it re-derives the
                // two volume terms (vm flood vs internal boxes)
                // to attribute the tilted-wall bias
                scalar splitVmFlood(0);
                scalar splitInternalBox(0);
                {
                    wallSubContactInfo& sC(*subContacts[0]);

                    // wall half-spaces of the sub-contact, as the
                    // production call builds them (wallContact.C)
                    List<planePolyClip::halfSpace> wallPlanes;
                    forAll(sC.getContactPatches(), cP)
                    {
                        List<vector> planeInfo
                        (
                            wallPlaneInfo::getWallPlaneInfo()
                            [sC.getContactPatches()[cP]]
                        );
                        wallPlanes.append
                        (
                            planePolyClip::halfSpace
                            (
                                planeInfo[1],
                                planeInfo[0]
                            )
                        );
                    }

                    // the vm-flood term: every contact VM of the
                    // sub-contact re-flooded (same lattices and
                    // wall planes as the production call)
                    for
                    (
                        label vmI = 0;
                        vmI < sC.getVMContactSize();
                        vmI++
                    )
                    {
                        autoPtr<virtualMeshWallInfo>& vmWInfo
                        (
                            sC.getVMContactInfo(vmI)
                        );
                        if (!vmWInfo.valid())
                        {
                            continue;
                        }

                        virtualMeshWall virtMeshWall
                        (
                            vmWInfo(),
                            wallCntInfo.getcClass().getGeomModel()
                        );

                        virtMeshWall.setWallPlanes(wallPlanes);

                        if (virtMeshWall.detectFirstContactPoint())
                        {
                            splitVmFlood +=
                                virtMeshWall.evaluateContact()
                               *vmWInfo->getEmptyScale();
                        }
                    }

                    // the internal-box term: full reboxed
                    // volumes, exactly as the production sum
                    forAll(sC.getInternalElements(), sCII)
                    {
                        splitInternalBox +=
                            sC.getInternalElements()[sCII]
                            .second().volume();
                    }

                    // the missed-element scan: first level of the
                    // leg only (the SM keep-set is level-
                    // independent) and the exact arm only (the
                    // carried boxes are the raw element boxes
                    // post a1b563e). run on the tilted and the
                    // sweep legs (the decision inputs) and on the
                    // aligned leg (the scan self-check: detection
                    // is clean there, so missedVol must be ~0)
                    if (exact && lI == 0 && wI != 2)
                    {
                        missedElementScan
                        (
                            wallCntInfo.getcClass().getGeomModel(),
                            wallPlanes,
                            sC,
                            ccs,
                            wallOrientations[wI],
                            wallPoints[wI]
                        );
                    }
                }

                // geometric evaluation only: getWallContactVars
                // fills wallCntVars (no force integration needed)
                contactModel::getWallContactVars
                (
                    mesh,
                    wallCntInfo,
                    1e-4,
                    *subContacts[0]
                );

                const wallContactVars& vars
                (
                    subContacts[0]->getWallCntVars()
                );

                Info << "    split level " << level << " "
                     << (exact ? "exact" : "legacy")
                     << ": vmFlood " << splitVmFlood
                     << " internalBox " << splitInternalBox
                     << " sum " << splitVmFlood + splitInternalBox
                     << " (production V "
                     << vars.contactVolume_ << ")"
                     << endl;

                Info << "  level " << level << " "
                     << (exact ? "exact" : "legacy")
                     << ": svEdge " << svEdge
                     << "  V " << vars.contactVolume_
                     << " (err vs faceted " << vars.contactVolume_/VrefF - 1
                     << ")"
                     << "  A " << vars.contactArea_
                     << " (err vs faceted " << vars.contactArea_/ArefF - 1
                     << ")";
                Info << endl;
            }
        }
    }
    } // end wall-orientation sweep

    Info << nl << "vmWallContactBenchmark: " << nFailed
         << " failure(s)" << endl;

    return nFailed == 0 ? 0 : 1;
}

// ************************************************************************* //
