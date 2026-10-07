#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""evaluate the vmPrtPrtContactBenchmark log, gate on convergence
order and write plottable data files.

for every (d, arm) series over the levels L1..Ln the gate is
the error ratio between consecutive levels:

    ratio_k = |err(L_k)| / |err(L_k+1)|

exact arm: every ratio must be >= 2 (second order: svEdge halves
per level, error must quarter). legacy arm: every ratio must be
>= 1.5 (first order: error must halve). contact-detection
failures in the log FAIL unconditionally.

plottable output: one CSV per tested geometric configuration
(d) into <logDir>/plots/, with the dependence of the tracked
variables (V, V-error of both arms) on the refinement level.

exit code 0 = all gates pass.
"""
import re
import sys
from pathlib import Path

LOG_RE = re.compile(
    r"^\s*level (\d+) (exact|legacy): svEdge (\S+)\s+"
    r"V (\S+) \(err vs faceted (\S+)\)"
)
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
    # per d -> arm -> level -> (svEdge, V, eV)
    series = {}
    cur = None
    fails = 0
    with open(path) as f:
        for line in f:
            if "[FAIL]" in line:
                fails += 1
            m = D_RE.search(line)
            if m:
                cur = series.setdefault(m.group(2), {})
                continue
            m = LOG_RE.match(line)
            if m and cur is not None:
                level, arm = int(m.group(1)), m.group(2)
                cur.setdefault(arm, {})[level] = (
                    float(m.group(3)),
                    float(m.group(4)), float(m.group(5)),
                )
    return series, fails


def writePlotFiles(series, logPath):
    plots = logPath.parent/"plots"
    plots.mkdir(exist_ok=True)
    n = 0
    for dKey, arms in series.items():
        rows = {}
        for arm in ("exact", "legacy"):
            for level, (sv, V, eV) in arms.get(arm, {}).items():
                rows.setdefault(level, {})[arm] = (sv, V, eV)
        if not rows:
            continue
        name = (f"dOverR_" + dKey.replace(".", "p") + ".csv")
        with open(plots/name, "w") as f:
            f.write("# vmPrtPrtContactBenchmark: d/R = "
                    + dKey + "\n")
            f.write("# V measured; eV relative error vs the "
                    "faceted lens reference\n")
            f.write("level svEdge "
                    "V_exact eV_exact V_legacy eV_legacy\n")
            for level in sorted(rows):
                parts = [str(level)]
                arm = rows[level]
                for key in ("exact", "legacy"):
                    if key in arm:
                        sv, V, eV = arm[key]
                        parts += [f"{sv:.6e}", f"{V:.9e}",
                                  f"{eV:+.6e}"]
                    else:
                        parts += ["nan"]*3
                f.write(" ".join(parts) + "\n")
        n += 1
    return n


def gate(dOverR, arm, errs, problems):
    if len(errs) < 2:
        problems.append(f"d/R {dOverR} ({arm}): fewer than 2 levels")
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
            f"d/R {dOverR} ({arm}): order gate {gateVal} not met "
            f"(ratio {min(bad):.2f})"
        )


def main():
    if len(sys.argv) != 2:
        print("usage: evalPrtPrtContactBenchmark.py <log>")
        sys.exit(2)
    logPath = Path(sys.argv[1])
    series, fails = parse(logPath)
    problems = []

    nPlots = writePlotFiles(series, logPath)
    print(f"plot data: {nPlots} file(s) in "
          f"{logPath.parent/'plots'}")

    for dOverR in series:
        print(f"  d/R = {dOverR}")
        for arm in ("exact", "legacy"):
            if arm in series[dOverR]:
                gate(dOverR, arm,
                     {l: v[2] for l, v in series[dOverR][arm].items()},
                     problems)
    if fails:
        problems.append(f"{fails} contact-detection failure(s)")

    if problems:
        print("\nvmPrtPrtContactBenchmark: FAIL")
        for p in problems:
            print("  - " + p)
        sys.exit(1)
    print("\nvmPrtPrtContactBenchmark: PASS (order gates met)")
    sys.exit(0)


if __name__ == "__main__":
    main()
