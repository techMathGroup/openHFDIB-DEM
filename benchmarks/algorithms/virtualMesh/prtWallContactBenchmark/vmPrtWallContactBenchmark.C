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
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    openHFDIB-DEM is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with openHFDIB-DEM. If not, see <http://www.gnu.org/licenses/>.

Application

    vmPrtWallContactBenchmark

Description

    Benchmark of the particle-particle contact geometry accuracy on
    the wallContactBenchmark geometric setup: the icosphere (STL)
    pressed against a planar STL body (a closed box whose top face
    lies at y = 0), evaluated through virtualMesh directly - the
    same detect/evaluate pair used by prtContact.C, but on two
    STL bodies, so the leaf rules run through the shapesIn path of
    stlBased (single-facet planes, rim leaves where both bodies are
    MIXED). this complements wallContactBenchmark (the production
    wall chain, plane-clipped band) and prtPrtContactBenchmark
    (sphere on sphere, both bodies faceted everywhere).

    References: the smooth sphere cap

        V = pi * xH^2 * (3*R - xH) / 3,  xH = R - d
        A = pi * (R^2 - d^2)

    AND the exact faceted cap of the STL polyhedron below the
    plane (clip every triangle to y <= 0, divergence-theorem volume;
    convex section polygon area) - identical to the reference of
    wallContactBenchmark, since the plate is thicker than the
    deepest penetration, the overlap of sphere and plate is exactly
    the sphere cap.

    For each (d, level) the benchmark reports the measured contact
    volume and area of both arms (exactSubVolume true and false -
    the identical-binary A/B pair) with errors against both
    references. the area comes from the production prt-prt sector
    sweep (get3DcontactNormalAndSurface, nonConvex path); its error
    is reported but the gate is the convergence ORDER of the volume
    only (evaluator script). The tilted leg (30 degrees about z,
    sphere and plate rotated together) is reported but not gated.

    Run from the case directory of this benchmark (needs
    constant/triSurface/{sphere,plate}.stl and a blockMesh mesh).

Contributors
    Martin Isoz (2019-*)
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "nonConvexBody.H"
#include "triSurface.H"
#include "virtualMesh.H"
#include "virtualMeshLevel.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

namespace
{

//- exact faceted references of the STL cap below y = 0: the
// sphere STL shifted to center (0, d, 0), clipped to y <= 0
// (identical to the facetedReference of vmWallContactBenchmark -
// the plate removes nothing the plane would not)
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

}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "convergence study of the prt-prt contact geometry "
        "(sphere cap against a planar STL body) against the "
        "faceted and analytic references"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // --------------------------------------------------------------
    // case parameters of the study
    // --------------------------------------------------------------

    // analytic radius of the test sphere (the STL was generated
    // with exactly this radius)
    const scalar R(0.0075);

    // center-to-plane distance d of the sweep (the plate top face
    // plays the wall of wallContactBenchmark)
    const FixedList<scalar,3> dList({0.95*R, 0.75*R, 0.5*R});

    // virtualMesh levels of the sweep: svEdge halves each step
    const FixedList<label,3> levelList({3, 4, 5});

    const scalar charCellSize(0.001);
    const scalar maxSubVolumes(1000000000);

    // leg orientations: leg 0 is the lattice-aligned plate; leg 1
    // the 30-degree tilt about z (sphere and plate rotated
    // together, so the cap geometry relative to the plane is
    // identical, but every leaf sees a diagonal plate plane)
    const FixedList<scalar,2> legAngles({0.0, 30.0});

    virtualMeshLevel::setMaxSubVolumes(maxSubVolumes);

    // the sphere and the plate share the STL coordinates; both are
    // placed by shifting their surface points
    autoPtr<nonConvexBody> sphereBody
    (
        new nonConvexBody
        (
            mesh,
            word("constant/triSurface/sphere.stl")
        )
    );
    autoPtr<nonConvexBody> plateBody
    (
        new nonConvexBody
        (
            mesh,
            word("constant/triSurface/plate.stl")
        )
    );

    // raw surface for the faceted references
    const triSurface refSurf("constant/triSurface/sphere.stl");

    // --------------------------------------------------------------
    // per (leg, d, level) evaluation
    // --------------------------------------------------------------

    label nFailed(0);

    Info << "vmPrtWallContactBenchmark: sphere-on-plate "
         << "convergence study" << endl;

    forAll(legAngles, wI)
    {
        if (wI > 0)
        {
            Info << nl
                 << "=== tilted leg (30 degrees about z) ===" << endl;
        }

        // shared rotation of sphere and plate about z through the
        // origin: the cap relative to the plate top face is the
        // same on both legs, only the leaf planes differ
        const scalar a
        (
            legAngles[wI]*constant::mathematical::pi/180.0
        );
        const tensor rot
        (
            Foam::cos(a), -Foam::sin(a), 0,
            Foam::sin(a),  Foam::cos(a), 0,
            0, 0, 1
        );

        // place the plate once per leg: rotate its STL coordinates
        // (top face at y = 0, containing the origin) about z through
        // the origin; on the aligned leg rot is the identity
        {
            pointField platePos(plateBody->getBodyPoints());
            forAll(platePos, pI)
            {
                platePos[pI] = rot & platePos[pI];
            }
            plateBody->setBodyPosition(platePos);
        }

        forAll(dList, dI)
        {
            const scalar d(dList[dI]);
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

            // position the sphere for this d: the un-rotated STL
            // center goes to (0, d, 0), then the shared leg
            // rotation maps it to rot & (0, d, 0) - moving the
            // body points along the rotated wall normal
            {
                const vector target
                (
                    rot & vector(0, d, 0)
                );
                const vector delta
                (
                    target - sphereBody->getCoM()
                );
                pointField pos(sphereBody->getBodyPoints());
                forAll(pos, pI)
                {
                    pos[pI] += delta;
                }
                sphereBody->setBodyPosition(pos);
            }

            forAll(levelList, lI)
            {
                const label level(levelList[lI]);
                virtualMeshLevel::setVirtualMeshLevel
                (
                    level,
                    charCellSize
                );
                const scalar svEdge
                (
                    charCellSize
                    /virtualMeshLevel::getLevelOfDivision()
                );

                for (label exactArm = 0; exactArm <= 1; exactArm++)
                {
                    const bool exact(exactArm == 1);
                    virtualMeshLevel::setExactSubVolume(exact);

                    // pair bbox: sphere bbox clipped to the plate
                    // bbox (the limitBBox rule of prtContactInfo)
                    boundBox pairBBox
                    (
                        max
                        (
                            sphereBody->getBounds().min(),
                            plateBody->getBounds().min()
                        ),
                        min
                        (
                            sphereBody->getBounds().max(),
                            plateBody->getBounds().max()
                        )
                    );

                    if (!pairBBox.valid())
                    {
                        Info << "  level " << level << " "
                             << (exact ? "exact" : "legacy")
                             << ": [FAIL] no pair-bbox overlap"
                             << endl;
                        nFailed++;
                        continue;
                    }

                    const scalar subVolumeV(sqr(svEdge)*svEdge);

                    virtualMeshInfo vmInfo(pairBBox, subVolumeV);

                    virtualMesh virtMesh
                    (
                        vmInfo,
                        *(sphereBody.get()),
                        *(plateBody.get())
                    );

                    if (!virtMesh.detectFirstContactPoint())
                    {
                        Info << "  level " << level << " "
                             << (exact ? "exact" : "legacy")
                             << ": [FAIL] no contact detected"
                             << endl;
                        nFailed++;
                        continue;
                    }

                    const scalar V(virtMesh.evaluateContact());

                    if (V < VSMALL)
                    {
                        Info << "  level " << level << " "
                             << (exact ? "exact" : "legacy")
                             << ": [FAIL] contact detected "
                             << "but V = 0" << endl;
                        nFailed++;
                        continue;
                    }

                    // area: the production prt-prt sector sweep
                    // (both bodies are nonConvex, so the
                    // sub-contact clustering path runs)
                    const Tuple2<scalar, vector> areaAndNormal
                    (
                        virtMesh.get3DcontactNormalAndSurface(true)
                    );
                    const scalar A(areaAndNormal.first());

                    Info << "  level " << level << " "
                         << (exact ? "exact" : "legacy")
                         << ": svEdge " << svEdge
                         << "  V " << V
                         << " (err vs faceted " << V/VrefF - 1 << ")"
                         << "  A " << A
                         << " (err vs faceted " << A/ArefF - 1 << ")"
                         << endl;
                }
            }
        }
    }

    Info << nl << "vmPrtWallContactBenchmark: " << nFailed
         << " failure(s)" << endl;

    return nFailed == 0 ? 0 : 1;
}

// ************************************************************************* //
