// Verification check 1 of the surface-loading model: the law table.
//
// Builds the kapton-flyer / aluminium-sample pair by hand (no OpenFOAM, no
// registry, no mesh) and prints F, p_mean, p_H, f, Z at a set of velocities
// that bracket every regime of the law. Compare with reference.py, which
// implements the same equations independently (Newton instead of bisection).
//
//   g++ -O2 -I<repo>/src/HFDIBDEM/contactModels \
//       -o surfaceLoadingLawTable surfaceLoadingLawTable.C
//   ./surfaceLoadingLawTable > table.txt
//   python3 -I reference.py table.txt

#include "surfaceLoadingLaw.H"

#include <cstdio>

using namespace Foam;

int main()
{
    surfaceLoadingLaw::pairParams p;

    // kapton flyer, Katagiri PRB 105 054103 (2022), two-segment Hugoniot
    p.flyer.rho = 1420.0;
    p.flyer.Yc = 0.12e9;
    p.flyer.capFactor = 1.25;
    p.flyer.hasHugoniot = true;
    p.flyer.c0 = 2327.0;
    p.flyer.s = 1.55;
    p.flyer.c0High = 1790.0;
    p.flyer.sHigh = 1.43;
    p.flyer.uMelt = 2780.0;

    // aluminium sample
    p.sample.rho = 2700.0;
    p.sample.Yc = 0.276e9;
    p.sample.capFactor = 2.7;
    p.sample.hasHugoniot = true;
    p.sample.c0 = 5350.0;
    p.sample.s = 1.34;

    // series modulus of the pair exactly as prtContactInfo computes it
    const double Yf(2.0e9), nuf(0.3), Ys(75.0e9), nus(0.33);
    p.aY = 1.0/((1.0 - nuf*nuf)/Yf + (1.0 - nus*nus)/Ys);
    p.reduceM = 2.2187494e-7;   // thin-case flyer (geometric mass)
    p.reduceBeta = 0.5;         // exercise the damping term

    const double A(6.1527e-6);   // thin-case footprint [m2]
    const double w(1.0e-6);      // 1 um mean overlap
    const double V(A*w);
    const double Lc(1.0e-4);     // Lc of the thin case

    std::printf("# aY %.17g A %.17g w %.17g Lc %.17g\n", p.aY, A, w, Lc);
    std::printf("# v F pMean pHug pStr pEdge f Z kEff\n");

    const double vv[] =
    {
        1.0, 20.0, 56.0, 77.0, 100.0, 300.0,
        1000.0, 1884.0, 2700.0, 4100.0, 4170.0, 4200.0, 4300.0
    };

    for (const double v : vv)
    {
        const surfaceLoadingLaw::result r
        (
            surfaceLoadingLaw::evaluate(p, v, A, V, Lc)
        );

        std::printf
        (
            "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
            v, r.F, r.pMean, r.pHug, r.pStr, r.pEdge, r.blendF, r.zEff,
            r.kEff
        );
    }

    // elastic limit of the same pair with the cap and blend off: the
    // force must be exactly the layer spring aY*V/Lc
    surfaceLoadingLaw::pairParams legacy(p);
    legacy.blendOn = false;
    legacy.capOn = false;
    legacy.halfSpace = false;

    const surfaceLoadingLaw::result rl
    (
        surfaceLoadingLaw::evaluate(legacy, 50.0, A, V, Lc)
    );

    const double damp
    (
        p.reduceBeta*std::sqrt(p.aY*p.reduceM*A/(Lc + 1.0e-15))
    );

    std::printf
    (
        "# legacy %.17g %.17g\n",
        rl.F, p.aY*V/(Lc + 1.0e-15) + damp*50.0
    );

    return 0;
}
