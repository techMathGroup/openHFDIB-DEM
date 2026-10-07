#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""plot the virtual-mesh contact benchmark results.

discovers the CSV plot-data files written by the benchmark
evaluators (wallContactBenchmark and prtPrtContactBenchmark,
case/plots/*.csv) relative to this script's location and renders
one PNG per configuration plus one error-summary figure per
benchmark:

  - per configuration: measured V (and A, wall) of both arms and
    their errors vs the refinement level
  - summary: all configurations' errors on a common log-log
    axis against svEdge, with first/second-order slope guides

run after Allrun.sh of the benchmarks; requires matplotlib.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt     # noqa: E402

BENCHMARKS = ["wallContactBenchmark", "prtPrtContactBenchmark",
              "prtWallContactBenchmark"]

# order-guide slopes on the log-log error plot: err ~ svEdge^p
GUIDE_SLOPES = [1, 2]


def readCsv(path):
    """read one benchmark plot CSV: returns (configName, colNames,
    rows) where rows is a list of float lists ('nan' allowed)"""
    lines = [l.strip() for l in path.read_text().splitlines()
             if l.strip()]
    header = next(l for l in lines if not l.startswith("#"))
    cols = header.split()
    rows = []
    for l in lines[lines.index(header) + 1:]:
        rows.append([float(x) for x in l.split()])
    config = lines[0].lstrip("# ").strip()
    return config, cols, rows


def col(rows, name, cols):
    i = cols.index(name)
    return [r[i] for r in rows]


def plotConfig(csvPath, outDir):
    """one figure per configuration: measured values and errors of
    both arms vs the refinement level (errors on a symlog axis -
    the legacy arm changes sign)"""
    config, cols, rows = readCsv(csvPath)
    hasArea = "A_exact" in cols
    level = col(rows, "level", cols)

    panels = [("V", "contact volume", False),
              ("eV", "volume error vs faceted", True)]
    if hasArea:
        panels += [("A", "wetted area", False),
                   ("eA", "area error vs faceted", True)]

    fig, axes = plt.subplots(1, len(panels),
                             figsize=(4*len(panels), 3.5),
                             sharex=True)
    if len(panels) == 1:
        axes = [axes]
    fig.suptitle(config)

    for ax, (key, label, symlog) in zip(axes, panels):
        for arm, style in (("exact", "o-"), ("legacy", "s--")):
            name = f"{key}_{arm}"
            if name in cols:
                ax.plot(level, col(rows, name, cols), style,
                        label=arm)
        ax.set_xlabel("refinement level")
        ax.set_ylabel(label)
        if symlog:
            ax.set_yscale("symlog", linthresh=1e-8)
        else:
            ax.grid(True, alpha=0.3)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)

    out = outDir/(csvPath.stem + ".png")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def plotSummary(csvPaths, outDir, benchName):
    """all configurations' volume errors vs svEdge, log-log, with
    order slope guides"""
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    fig.suptitle(benchName + ": contact-volume error vs refinement")

    for arm, ax in zip(("exact", "legacy"), axes):
        for p in csvPaths:
            config, cols, rows = readCsv(p)
            if config.find("tilted") > -1: continue
            sv = col(rows, "svEdge", cols)
            e = [abs(x) for x in
                 col(rows, f"eV_{arm}", cols)]
            ax.loglog(sv, e, "o-",
                      label=config.split(" d/R = ")[-1])
        ax.set_xlabel("svEdge")
        ax.set_ylabel(f"|eV| ({arm})")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=7)
        # slope guides anchored at the top-left data corner
        if len(sv) >= 2:
            x0, y0 = sv[0], max(e[0], 1e-6)
            for p in GUIDE_SLOPES:
                ax.loglog([x0, sv[-1]],
                          [y0, y0*(sv[-1]/x0)**p],
                          "k:", alpha=0.5)
                ax.annotate(f"slope -{p}",
                            (sv[-1], y0*(sv[-1]/x0)**p),
                            fontsize=7, alpha=0.6)

    out = outDir/(benchName + "_summary.png")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main():
    here = Path(__file__).resolve().parent
    missing = []
    nFig = 0
    for bench in BENCHMARKS:
        plotsDir = here/bench/"case"/"plots"
        csvs = sorted(plotsDir.glob("*.csv")) if plotsDir.is_dir() \
            else []
        if not csvs:
            missing.append(bench)
            continue
        print(f"{bench}: {len(csvs)} configuration(s)")
        for c in csvs:
            print(f"  {plotConfig(c, plotsDir)}")
            nFig += 1
        print(f"  {plotSummary(csvs, plotsDir, bench)}")
        nFig += 1
    if missing:
        print("no plot data found for: "
              + ", ".join(missing)
              + " (run Allrun.sh first)")
        if nFig == 0:
            sys.exit(1)
    print(f"{nFig} figure(s) written")


if __name__ == "__main__":
    main()
