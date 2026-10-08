/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| |  \/ |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | |_/ /---| | | || |_  | |\/| |
 / _ \| "_ \ / _ \| "_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  /\| |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___/ |_____|_|  |_|
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

Application
    vmPrtPrtContactBenchmark

Description
    Benchmark of the particle-particle contact geometry accuracy:
    two icospheres (same STL) overlapping along x, evaluated through
    virtualMesh directly (the detect/evaluate pair used by
    prtContact.C). No contact model is involved - this isolates the
    evaluateLeafExact leaf rule, in particular the two-half-space clip
    of rim leaves (both bodies MIXED).

    References: the smooth sphere lens

        V = pi * (2R - d)^2 * (d^2 + 4 d R) / (12 d)

    AND the exact faceted lens of the STL polyhedron pair (clip every
    triangle to x <= 0, divergence-theorem volume; the overlap is the
    reflected union of the two clipped halves). The faceted reference
    is the primary one - the virtual mesh resolves the STL bodies,
    not the smooth spheres, and at icosphere faceting the two differ
    by percent-level offsets that would otherwise mask the svEdge
    convergence.

    For each (d, level) the benchmark reports the measured contact
    volume of both arms (exactSubVolume true and false - the
    identical-binary A/B pair) with errors against both references.
    The gate is the convergence ORDER only (evaluator script): the
    exact arm must halve-plus its error per level refinement (second
    order), the legacy arm is expected first order by design.

    Cross-check: tests/virtualMesh/vmReplica.py run on the identical
    STL pair and settings must reproduce the exact-arm numbers.

Contributors
    Martin Isoz (2019-*)
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "nonConvexBody.H"
#include "triSurface.H"
#include "virtualMesh.H"
#include "virtualMeshLevel.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

//- exact faceted lens reference of the STL pair: the overlap of the
// polyhedron centered at (-d/2, 0, 0) with its mirror image. the
// overlap's right half is the part of the left polyhedron with
// x >= 0 (every point of the mirrored body with x <= 0 already
// lies inside the left body), so the lens volume is twice the
// volume of the polyhedron clipped to x >= 0; volume by clipping
// every triangle to x >= 0 and integrating with the divergence
// theorem (clipped faces keep their winding)
void facetedLensReference
(
    const triSurface& surf,
    const scalar& d,
    scalar& lensV
)
{
    lensV = 0;

    forAll(surf, fI)
    {
        const triFace& f(surf[fI]);
        pointField pts(3);
        for (label i = 0; i < 3; i++)
        {
            const point& p(surf.localPoints()[f[i]]);
            pts[i] = point(p.x() - 0.5*d, p.y(), p.z());
        }

        // sutherland-hodgman: keep x >= 0
        DynamicList<point> poly;
        for (label i = 0; i < 3; i++)
        {
            const point& p1(pts[i]);
            const point& p2(pts[(i + 1) % 3]);
            const bool in1(p1.x() >= 0);
            const bool in2(p2.x() >= 0);
            if (in1)
            {
                poly.append(p1);
            }
            if (in1 != in2)
            {
                const scalar t((0 - p1.x())/(p2.x() - p1.x()));
                poly.append(p1 + t*(p2 - p1));
            }
        }

        // fan triangulation from vertex 0 keeps the winding
        for (label k = 1; k < poly.size() - 1; k++)
        {
            const point& a(poly[0]);
            const point& b(poly[k]);
            const point& c(poly[k + 1]);
            lensV += ((a ^ b) & c)/6.0;
        }
    }

    lensV *= 2;
}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "benchmark of the prt-prt contact geometry accuracy "
        "(sphere lens) against the faceted and analytic references"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // analytic radius of the test sphere (the STL was generated with
    // exactly this radius, centered at the origin)
    const scalar R(1.0);

    // center distance sweep: near touching to half overlap
    const FixedList<scalar,4> dList({1.9*R, 1.6*R, 1.3*R, 1.0*R});

    // virtualMesh levels of the sweep
    const FixedList<label,4> levelList({2, 3, 4, 5});

    const scalar charCellSize(2.0/15.0);
    const scalar maxSubVolumes(1000000000);

    virtualMeshLevel::setMaxSubVolumes(maxSubVolumes);

    // the two bodies share one STL; each is positioned by shifting
    // the surface points to its center
    autoPtr<nonConvexBody> cBody
    (
        new nonConvexBody
        (
            mesh,
            word("constant/triSurface/sphere.stl")
        )
    );
    autoPtr<nonConvexBody> tBody
    (
        new nonConvexBody
        (
            mesh,
            word("constant/triSurface/sphere.stl")
        )
    );

    // raw surface for the faceted lens reference
    const triSurface refSurf("constant/triSurface/sphere.stl");

    label nFailed(0);

    Info << "vmPrtPrtContactBenchmark: sphere-lens convergence study"
         << endl;

    forAll(dList, dI)
    {
        const scalar d(dList[dI]);
        const scalar VrefAna
        (
            constant::mathematical::pi*sqr(2*R - d)
           *(d*d + 4*d*R)/(12*d)
        );
        scalar VrefF(0);
        facetedLensReference(refSurf, d, VrefF);

        Info << nl << "--- center distance d = " << d
             << " (d/R = " << d/R << ")" << endl;
        Info << "    analytic lens volume " << VrefAna
             << ", faceted lens volume " << VrefF
             << " (faceting " << VrefF/VrefAna - 1 << ")" << endl;

        forAll(levelList, lI)
        {
            const label level(levelList[lI]);
            virtualMeshLevel::setVirtualMeshLevel(level, charCellSize);
            const scalar svEdge
            (
                charCellSize/virtualMeshLevel::getLevelOfDivision()
            );

            for (label exactArm = 0; exactArm <= 1; exactArm++)
            {
                const bool exact(exactArm == 1);
                virtualMeshLevel::setExactSubVolume(exact);

                // position: STL center to (+/- d/2, 0, 0)
                {
                    const vector cShift
                    (
                        vector(-0.5*d, 0, 0) - cBody->getCoM()
                    );
                    const vector tShift
                    (
                        vector(0.5*d, 0, 0) - tBody->getCoM()
                    );

                    pointField cPos(cBody->getBodyPoints());
                    forAll(cPos, pI)
                    {
                        cPos[pI] += cShift;
                    }
                    cBody->setBodyPosition(cPos);

                    pointField tPos(tBody->getBodyPoints());
                    forAll(tPos, pI)
                    {
                        tPos[pI] += tShift;
                    }
                    tBody->setBodyPosition(tPos);
                }

                // pair bbox: c bbox clipped to the t bbox (the
                // limitBBox rule of prtContactInfo)
                boundBox pairBBox
                (
                    max
                    (
                        cBody->getBounds().min(),
                        tBody->getBounds().min()
                    ),
                    min
                    (
                        cBody->getBounds().max(),
                        tBody->getBounds().max()
                    )
                );

                if (!pairBBox.valid())
                {
                    Info << "  level " << level << " "
                         << (exact ? "exact" : "legacy")
                         << ": [FAIL] no pair-bbox overlap" << endl;
                    nFailed++;
                    continue;
                }

                const scalar subVolumeV(sqr(svEdge)*svEdge);

                virtualMeshInfo vmInfo(pairBBox, subVolumeV);

                virtualMesh virtMesh
                (
                    vmInfo,
                    *(cBody.get()),
                    *(tBody.get())
                );

                if (!virtMesh.detectFirstContactPoint())
                {
                    Info << "  level " << level << " "
                         << (exact ? "exact" : "legacy")
                         << ": [FAIL] no contact detected" << endl;
                    nFailed++;
                    continue;
                }

                const scalar V(virtMesh.evaluateContact());

                if (V < VSMALL)
                {
                    Info << "  level " << level << " "
                         << (exact ? "exact" : "legacy")
                         << ": [FAIL] contact detected but V = 0"
                         << endl;
                    nFailed++;
                    continue;
                }

                Info << "  level " << level << " "
                     << (exact ? "exact" : "legacy")
                     << ": svEdge " << svEdge
                     << "  V " << V
                     << " (err vs faceted " << V/VrefF - 1 << ")"
                     << " (err vs analytic " << V/VrefAna - 1 << ")"
                     << endl;
            }
        }
    }

    Info << nl << "vmPrtPrtContactBenchmark: " << nFailed
         << " failure(s)" << endl;

    return nFailed == 0 ? 0 : 1;
}

// ************************************************************************* //
