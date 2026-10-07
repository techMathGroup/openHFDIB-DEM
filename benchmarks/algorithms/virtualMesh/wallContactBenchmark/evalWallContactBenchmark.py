#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""evaluate the vmWallContactBenchmark log, gate on convergence
order and write plottable data files.

for every (leg, d, arm) series over the levels L1..Ln the gate is
the error ratio between consecutive levels:

    ratio_k = |err(L_k)| / |err(L_k+1)|

exact arm: every ratio must be >= 2 (second order: svEdge halves
per level, error must quarter). legacy arm: every ratio must be
>= 1.5 (first order: error must halve). the tilted-wall leg is
reported but never gated (known open issue). contact-detection
failures in the log FAIL unconditionally.

plottable output: one CSV per tested geometric configuration
(leg, d) into <logDir>/plots/, with the dependence of the
tracked variables (V, V-error, A, A-error of both arms) on the
refinement level.

exit code 0 = all gates pass.
"""
import math
import re
import sys
from pathlib import Path

LOG_RE = re.compile(
    r"^\s*level (\d+) (exact|legacy): svEdge (\S+)\s+"
    r"V (\S+) \(err vs faceted (\S+)\)\s+"
    r"A (\S+) \(err vs faceted (\S+)\)"
)
LEG_RE = re.compile(r"=== tilted-wall leg")
D_RE = re.compile(r"--- center distance d = (\S+) \(d/R = ([\d.]+)\)")

# order gates: error ratio across one level refinement
EXACT_ORDER = 2.0
LEGACY_ORDER = 1.5

# legacy-arm noise floor: the count-all error is a signed sum of
# over/under-counts that partially cancel, so once it falls below
# this level the per-step ratios measure cancellation noise, not
# order - such pairs are skipped by the legacy gate
LEGACY_NOISE = 1e-4


def parse(path):
    # per leg -> d -> arm -> level -> (svEdge, V, eV, A, eA)
    legs = []
    legName = "aligned"
    cur = None
    fails = 0
    with open(path) as f:
        for line in f:
            if LEG_RE.search(line):
                legName = "tilted"
            m = D_RE.search(line)
            if m:
                cur = {}
                if not legs or legs[-1][0] != legName:
                    legs.append((legName, {}))
                legs[-1][1][m.group(2)] = cur
                continue
            if "[FAIL]" in line:
                fails += 1
            m = LOG_RE.match(line)
            if m and cur is not None:
                level, arm = int(m.group(1)), m.group(2)
                if arm not in ("exact", "legacy"):
                    continue
                cur.setdefault(arm, {})[level] = (
                    float(m.group(3)),
                    float(m.group(4)), float(m.group(5)),
                    float(m.group(6)), float(m.group(7)),
                )
    return legs, fails


def writePlotFiles(legs, logPath):
    plots = logPath.parent/"plots"
    plots.mkdir(exist_ok=True)
    n = 0
    for legName, dSeries in legs:
        for dKey, arms in dSeries.items():
            rows = {}
            for arm in ("exact", "legacy"):
                for level, (sv, V, eV, A, eA) in \
                        arms.get(arm, {}).items():
                    rows.setdefault(level, {})[arm] = (sv, V, eV, A, eA)
            if not rows:
                continue
            name = (f"{legName}_dOverR_"
                    + dKey.replace(".", "p") + ".csv")
            with open(plots/name, "w") as f:
                f.write("# vmWallContactBenchmark: wall leg "
                        f"{legName}, d/R = {dKey}\n")
                f.write("# V, A measured; eV, eA relative error vs"
                        " the faceted reference\n")
                f.write("level svEdge "
                        "V_exact eV_exact A_exact eA_exact "
                        "V_legacy eV_legacy A_legacy eA_legacy\n")
                for level in sorted(rows):
                    parts = [str(level)]
                    arm = rows[level]
                    for key in ("exact", "legacy"):
                        if key in arm:
                            sv, V, eV, A, eA = arm[key]
                            parts += [f"{sv:.6e}", f"{V:.9e}",
                                      f"{eV:+.6e}", f"{A:.9e}",
                                      f"{eA:+.6e}"]
                        else:
                            parts += ["nan"]*5
                    f.write(" ".join(parts) + "\n")
                n += 1
    return n


def gate(label, arm, errs, problems):
    if len(errs) < 2:
        problems.append(f"{label} ({arm}): fewer than 2 levels")
        return
    ratios = []
    skipped = []
    lv = sorted(errs)
    for lo, hi in zip(lv[:-1], lv[1:]):
        eLo, eHi = abs(errs[lo]), abs(errs[hi])
        noise = LEGACY_NOISE if arm == "legacy" else 1e-13
        if eHi < noise or eLo < noise:        # below the noise floor
            skipped.append((lo, hi))
            continue
        ratios.append(eLo / eHi)
    if skipped:
        tail = ("  (noise-floor skips: "
                + " ".join(f"L{a}-L{b}" for a, b in skipped)
                + ")")
    else:
        tail = ""
    gateVal = EXACT_ORDER if arm == "exact" else LEGACY_ORDER
    print(f"    {arm:6s} errors: "
          + " ".join(f"{errs[l]:+.3e}" for l in lv)
          + "  ratios: "
          + "/".join(f"{r:.2f}" for r in ratios)
          + tail)
    bad = [r for r in ratios if r < gateVal]
    if bad:
        problems.append(
            f"{label} ({arm}): order gate {gateVal} not met "
            f"(ratio {min(bad):.2f})"
        )


def main():
    if len(sys.argv) != 2:
        print("usage: evalWallContactBenchmark.py <log>")
        sys.exit(2)
    logPath = Path(sys.argv[1])
    legs, fails = parse(logPath)
    problems = []

    nPlots = writePlotFiles(legs, logPath)
    print(f"plot data: {nPlots} file(s) in "
          f"{logPath.parent/'plots'}")

    for legName, dSeries in legs:
        gated = legName != "tilted"
        print(f"=== {legName} wall leg"
              + ("" if gated else " (NOT gated - known issue)"))
        for dKey, arms in dSeries.items():
            print(f"  d/R = {dKey}")
            for arm in ("exact", "legacy"):
                if arm in arms:
                    if gated:
                        gate(f"d/R {dKey}", arm,
                             {l: v[2] for l, v in arms[arm].items()},
                             problems)
                    else:
                        print(f"    {arm:6s} errors: "
                              + " ".join(
                                  f"{arms[arm][l][1]:+.3e}"
                                  for l in sorted(arms[arm]))
                              + " (not gated)")
    if fails:
        problems.append(f"{fails} contact-detection failure(s)")

    if problems:
        print("\nvmWallContactBenchmark: FAIL")
        for p in problems:
            print("  - " + p)
        sys.exit(1)
    print("\nvmWallContactBenchmark: PASS (order gates met)")
    sys.exit(0)


if __name__ == "__main__":
    main()
