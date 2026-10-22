#!/usr/bin/env python3
"""Independent reference for the surface-loading law table.

Written from the equations, NOT translated from surfaceLoadingLaw.H:
the impedance match is solved in PRESSURE (u = u(P) by inverting the
Hugoniot analytically, bisection in P) whereas the C++ bisects in particle
velocity and evaluates P(u). Different variable, different arithmetic.

usage: python3 -I reference.py table.txt
exit 0 = pass
"""

import math
import sys

# --- pair definition, identical numbers to surfaceLoadingLawTable.C --------
FLYER = dict(rho=1420.0, Yc=0.12e9, cap=1.25,
             c0=2327.0, s=1.55, c0H=1790.0, sH=1.43, uMelt=2780.0)
SAMPLE = dict(rho=2700.0, Yc=0.276e9, cap=2.7,
              c0=5350.0, s=1.34, c0H=0.0, sH=0.0, uMelt=0.0)

A = 6.1527e-6
W = 1.0e-6
LC = 1.0e-4
KAPPA = 1.05
LO, HI = 0.3, 3.0

YF, NUF, YS, NUS = 2.0e9, 0.3, 75.0e9, 0.33
AY = 1.0 / ((1 - NUF**2) / YF + (1 - NUS**2) / YS)
REDUCE_M = 2.2187494e-7
REDUCE_BETA = 0.5


def damping_coeff():
    return REDUCE_BETA * math.sqrt(AY * REDUCE_M * A / (LC + 1e-15))


def p_of_u(m, u):
    if u <= 0.0:
        return 0.0
    high = m["c0H"] > 0 and u >= m["uMelt"]
    c0 = m["c0H"] if high else m["c0"]
    s = m["sH"] if high else m["s"]
    return m["rho"] * (c0 + s * u) * u


def impedance_match(a, b, v):
    """first branch crossing of P_a(u) = P_b(v-u), u in [0, v]

    The Hugoniot jumps DOWN at a melt transition, so P_a(u) is not
    monotone and the equations can have two roots; the physical one is
    the first reached from u = 0. Locate it with a dense scan (the
    C++ uses a coarse scan), then refine by regula falsi.
    """
    if v <= 0.0:
        return 0.0
    n = 20000
    up, fp = 0.0, p_of_u(a, 0.0) - p_of_u(b, v)
    for i in range(1, n + 1):
        u = v * i / n
        f = p_of_u(a, u) - p_of_u(b, v - u)
        if f >= 0.0:
            lo, hi = up, u
            flo = fp
            for _ in range(200):
                # regula falsi keeps the bracket; no derivative needed
                if hi - lo < 1e-18 * max(v, 1.0):
                    break
                mid = lo + (hi - lo) * (-flo) / (f - flo)
                if not (lo < mid < hi):
                    mid = 0.5 * (lo + hi)
                fm = p_of_u(a, mid) - p_of_u(b, v - mid)
                if fm < 0.0:
                    lo, flo = mid, fm
                else:
                    hi, f = mid, fm
            return p_of_u(a, 0.5 * (lo + hi))
        up, fp = u, f
    return p_of_u(a, v)


def elastic_z(a, b):
    za, zb = a["rho"] * a["c0"], b["rho"] * b["c0"]
    return za * zb / (za + zb)


def smoothstep(t):
    t = min(max(t, 0.0), 1.0)
    return t * t * (3.0 - 2.0 * t)


def evaluate(v, blend=True, cap=True, half_space=True):
    pstr = min(FLYER["cap"] * FLYER["Yc"], SAMPLE["cap"] * SAMPLE["Yc"])
    pedge = SAMPLE["cap"] * SAMPLE["Yc"]

    phug = impedance_match(FLYER, SAMPLE, v)
    z = phug / v if v > 0 else elastic_z(FLYER, SAMPLE)
    # no shock on a separating contact (phug <= 0), and log10 would
    # otherwise trap on the non-positive ratio
    f = 0.0
    if blend and phug > 0.0:
        f = smoothstep((math.log10(phug / pstr) - math.log10(LO))
                       / (math.log10(HI) - math.log10(LO)))

    pceil = (1.0 - f) * pstr + f * phug
    k = KAPPA * AY * math.sqrt(A)
    fspring = k * W if half_space else AY * A * W / (LC + 1e-15)
    F = fspring + f * phug * A

    # standard DEM damping, then the ceiling
    fd = damping_coeff() * v
    if fd < 0.0 and -fd > F:
        fd = -F
    F += fd

    if cap and F > pceil * A:
        F = pceil * A
    F = max(F, 0.0)

    return dict(F=F, pmean=F / A, phug=phug, pstr=pstr, pedge=pedge,
                f=f, z=z, k=k)


def main(path):
    rows = []
    legacy = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("# legacy"):
                legacy = [float(x) for x in line.split()[2:]]
                continue
            if line.startswith("#"):
                continue
            rows.append([float(x) for x in line.split()])

    worst = 0.0
    bad = 0
    for r in rows:
        v = r[0]
        ref = evaluate(v)
        names = ["F", "pmean", "phug", "pstr", "pedge", "f", "z", "keff"]
        for i, name in enumerate(names, start=1):
            got = r[i]
            exp = {"F": ref["F"], "pmean": ref["pmean"],
                   "phug": ref["phug"], "pstr": ref["pstr"],
                   "pedge": ref["pedge"], "f": ref["f"],
                   "z": ref["z"], "keff": ref["k"]}[name]
            denom = max(abs(exp), 1e-30)
            rel = abs(got - exp) / denom
            worst = max(worst, rel)
            if rel > 1e-9:
                bad += 1
                print(f"FAIL v={v:g} {name}: got {got:.15g} expected {exp:.15g} rel {rel:.3g}")

    if legacy is not None:
        exp = AY * A * W / (LC + 1e-15) + damping_coeff() * 50.0
        rel = abs(legacy[0] - exp) / exp
        worst = max(worst, rel)
        print(f"legacy spring: got {legacy[0]:.15g} expected {exp:.15g} rel {rel:.3g}")
        if rel > 1e-12:
            bad += 1
            print("FAIL legacy spring is not bit-exact")

    print(f"rows {len(rows)}, worst relative difference {worst:.3g}")
    print("PASS" if bad == 0 else f"FAIL ({bad} values)")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
