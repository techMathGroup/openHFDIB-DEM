#!/usr/bin/env python3
"""Exercise characteristic cell size lifecycle states in pure DEM."""

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


SIZE_RE = re.compile(
    r"characteristic cell size: "
    r"([0-9]+(?:\.[0-9]*)?(?:[eE][+-]?[0-9]+)?)"
)
UNAVAILABLE = "characteristic cell size: unavailable (no represented cells)"
FALLBACK = "has occupied cells but no surface cells;"
NONFINITE_RE = re.compile(r"\b(?:nan|inf)\b", re.IGNORECASE)


def replace_once(path, old, new):
    text = path.read_text()
    if text.count(old) != 1:
        raise AssertionError(f"unexpected fixture setting in {path}: {old}")
    path.write_text(text.replace(old, new))


def run(command, case, log_name, timeout=60):
    log_path = case/log_name
    env = os.environ.copy()
    env["PWD"] = str(case)
    with log_path.open("w") as log:
        try:
            result = subprocess.run(
                command, cwd=case, env=env, stdout=log,
                stderr=subprocess.STDOUT, check=False, timeout=timeout
            )
        except subprocess.TimeoutExpired as error:
            raise RuntimeError(
                f"command timed out: {' '.join(command)}; see {log_path}"
            ) from error
    if log_path.stat().st_size > 20*1024*1024:
        raise RuntimeError(f"log exceeded 20 MiB: {log_path}")
    if result.returncode:
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"see {log_path}"
        )


def execute(template, root, name, end_time, configure=None):
    case = root/name
    shutil.copytree(template, case)
    shutil.copytree(case/"0.org", case/"0")
    replace_once(
        case/"system/controlDict", "endTime         0.1;",
        f"endTime         {end_time};"
    )
    if configure is not None:
        configure(case)
    run(["blockMesh"], case, "log.blockMesh")
    run(["HFDIBDEMFoam"], case, "log.HFDIBDEMFoam")
    text = (case/"log.HFDIBDEMFoam").read_text(errors="replace")
    require_clean_end(case, text)
    return text


def require_clean_end(case, text):
    if not text.rstrip().endswith("End"):
        raise AssertionError(f"solver did not reach End in {case}")
    if "FOAM FATAL" in text:
        raise AssertionError(f"unexpected fatal error in {case}")
    if NONFINITE_RE.search(text):
        raise AssertionError(f"non-finite numeric diagnostic in {case}")
    if "Start DEM pos: 0 DEM step: 0" in text:
        raise AssertionError(f"zero-step DEM loop in {case}")


def numeric_sizes(text):
    return [float(value) for value in SIZE_RE.findall(text)]


def check_exit(text):
    sizes = numeric_sizes(text)
    if not sizes:
        raise AssertionError("exit case never produced a valid cell size")
    if any(abs(size - 0.01) > 1e-12 for size in sizes):
        raise AssertionError(f"unexpected represented-cell size: {sizes}")
    if text.count(UNAVAILABLE) < 1:
        raise AssertionError("exit case never entered the unrepresented state")
    first_unavailable = text.index(UNAVAILABLE)
    if not SIZE_RE.search(text, 0, first_unavailable):
        raise AssertionError("size did not transition from valid to unavailable")
    switched_off = len(re.findall(
        r"-- body [0-9]+ switched off: purged [0-9]+ contact pair", text
    ))
    if switched_off != 1:
        raise AssertionError(
            f"expected one clean body switch-off, found {switched_off}"
        )
    if "zero DEM stability estimate" in text:
        raise AssertionError("unrepresented body left a stale contact pair")
    return {
        "validSizeValues": len(sizes),
        "unavailableStates": text.count(UNAVAILABLE),
        "switchOffEvents": switched_off,
    }


def configure_fallback(case):
    path = case/"constant/HFDIBDEMDict"
    replace_once(
        case/"system/blockMeshDict", "hex (0 1 2 3 4 5 6 7) (10 10 2)",
        "hex (0 1 2 3 4 5 6 7) (1 1 1)"
    )
    replace_once(
        path, "        radius      0.025;", "        radius      0.08;"
    )
    replace_once(
        path, "startPosition (0.05 0.07 0);",
        "startPosition (0.05 0.05 0);"
    )
    replace_once(
        path, "    fullyCoupledBody;",
        "    prescribedTransBody\n    {\n        velocity (0 0 0);\n    }"
    )
    replace_once(
        path, "    bodyGeom    sphere;", "    bodyGeom    sphere;\n    bodyCreation legacy;"
    )
    replace_once(
        case/"constant/g", "value           (0 9.81 0);",
        "value           (0 0 0);"
    )


def check_fallback(text):
    sizes = numeric_sizes(text)
    if len(sizes) < 3:
        raise AssertionError(
            "fallback was not exercised over repeated body creations"
        )
    expected_size = (0.1*0.1*0.02)**(1.0/3.0)
    if any(abs(size/expected_size - 1.0) > 1e-8 for size in sizes):
        raise AssertionError(f"unexpected occupied-cell fallback size: {sizes}")
    warning_count = text.count(FALLBACK)
    if warning_count != 1:
        raise AssertionError(
            f"expected one fallback warning, found {warning_count}"
        )
    if UNAVAILABLE in text:
        raise AssertionError("occupied body was marked unrepresented")
    return {
        "fallbackCreations": len(sizes),
        "fallbackWarnings": warning_count,
        "characteristicSize": sizes[0],
    }


def configure_stl(case):
    path = case/"constant/HFDIBDEMDict"
    replace_once(
        path,
        "adaptiveStepDEM false;\nstepDEM         1;",
        "adaptiveStepDEM true;\n"
        "adaptiveStepDEMDict\n"
        "{\n"
        "    useDEMdtEstimator true;\n"
        "    DEMdtEstimator conservative;\n"
        "}\n"
        "stepDEM         1;"
    )
    replace_once(
        path, "    collisionPatches {}",
        "    collisionPatches\n"
        "    {\n"
        "        bot\n"
        "        {\n"
        "            type stlBased;\n"
        "            material wall;\n"
        "            stlName botWall;\n"
        "        }\n"
        "    }"
    )
    replace_once(
        path, "        radius      0.025;", "        radius      0.01;"
    )
    replace_once(
        path, "startPosition (0.05 0.07 0);",
        "startPosition (0.05 0.05 0);"
    )
    replace_once(
        case/"constant/g", "value           (0 9.81 0);",
        "value           (0 -9.81 0);"
    )


def check_stl(text):
    sizes = numeric_sizes(text)
    if not sizes or any(abs(size - 0.01) > 1e-12 for size in sizes):
        raise AssertionError(
            f"ordinary STL partner did not supply the cell size: {sizes}"
        )
    if "Created stlBased wall body" not in text:
        raise AssertionError("STL wall body was not created")
    if "characteristic cell size requested for an unrepresented body" in text:
        raise AssertionError("unsized STL wall reached the validated getter")
    if not re.search(r"DEM step: [1-9][0-9]*", text):
        raise AssertionError("STL case did not execute a positive DEM step")
    return {
        "ordinaryPartnerSizes": len(sizes),
        "stlWallPresent": True,
        "positiveDEMSteps": len(re.findall(r"DEM step: [1-9][0-9]*", text)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--keep", type=Path,
        help="retain generated cases under this empty/new directory"
    )
    args = parser.parse_args()
    for executable in ("blockMesh", "HFDIBDEMFoam"):
        if shutil.which(executable) is None:
            raise RuntimeError(f"missing executable: {executable}")
    if args.keep:
        args.keep.mkdir(parents=True, exist_ok=False)
        root = args.keep
        cleanup = None
    else:
        cleanup = tempfile.TemporaryDirectory(prefix="charCellSize-regression-")
        root = Path(cleanup.name)

    template = Path(__file__).resolve().parent/"case"
    try:
        summaries = {}
        summaries["validToUnrepresented"] = check_exit(
            execute(template, root, "exit", 0.14)
        )
        summaries["occupiedWithoutSurface"] = check_fallback(
            execute(template, root, "fallback", 0.002, configure_fallback)
        )
        summaries["unsizedStlWall"] = check_stl(
            execute(template, root, "stlWall", 0.002, configure_stl)
        )
        print(json.dumps(summaries, indent=2, sort_keys=True))
        print("PASS: characteristic cell size lifecycle states are guarded")
    finally:
        if cleanup is not None:
            cleanup.cleanup()


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, RuntimeError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)
