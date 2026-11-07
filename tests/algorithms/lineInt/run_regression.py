"""Exercise processor handoffs with the real lineInt implementation."""

import argparse
from collections import Counter
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


PRODUCTION_FIELDS = ("U", "p", "f", "lambda", "Ui")


def run(command, cwd, log_name, env):
    with (cwd/log_name).open("w") as log:
        result = subprocess.run(
            command, cwd=cwd, env=env, stdout=log,
            stderr=subprocess.STDOUT, check=False
        )
    if result.returncode:
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"see {cwd/log_name}"
        )


def read_probes(case):
    owner = []
    receive = []
    final = []
    for path in sorted(case.glob("marchProbe.*.dat")):
        with path.open() as stream:
            for line in stream:
                fields = line.split()
                if not fields:
                    continue
                if fields[0] == "MARCH_OWNER":
                    owner.append({
                        "time": float(fields[1]),
                        "order": int(fields[2]),
                        "origin_rank": int(fields[3]),
                        "origin_label": int(fields[4]),
                        "owner_rank": int(fields[5]),
                        "owner_cell": int(fields[6]),
                        "inside": int(fields[7])
                    })
                elif fields[0] == "MARCH_RECV":
                    receive.append({
                        "running": float(fields[1]) != 0,
                        "order": int(fields[2]),
                        "origin_rank": int(fields[3]),
                        "origin_label": int(fields[4])
                    })
                elif fields[0] == "MARCH_FINAL":
                    final.append({
                        "time": float(fields[1]),
                        "order": int(fields[2]),
                        "center": tuple(
                            round(float(x), 10) for x in fields[3:6]
                        ),
                        "valid": int(fields[6]) != -1,
                        "point": tuple(float(x) for x in fields[7:10])
                    })
    return {"owner": owner, "receive": receive, "final": final}


def histogram(final):
    result = {}
    for row in final:
        valid = result.setdefault(row["center"], [set(), set()])
        valid[int(row["time"] != 0)].add(row["order"])
    # Use sets so duplicate records cannot hide a missing completed point.
    return {
        key: {"initial": sorted(value[0]), "running": sorted(value[1])}
        for key, value in result.items()
    }


def pair_geometry(serial, parallel):
    """Match physical cells without depending on rounded center keys."""
    parallel_buckets = {}
    for row in parallel:
        parallel_buckets.setdefault(
            (row["time"] != 0, row["order"]), []
        ).append(row)
    pairs = []
    for serial_row in serial:
        key = (serial_row["time"] != 0, serial_row["order"])
        candidates = parallel_buckets.get(key, [])
        if not candidates:
            raise AssertionError(f"no parallel record for {key}")
        distance, index = min(
            (math.dist(serial_row["center"], row["center"]), i)
            for i, row in enumerate(candidates)
        )
        if distance > 1e-8:
            raise AssertionError(f"physical cell-center mismatch: {distance}")
        parallel_row = dict(candidates.pop(index))
        parallel_row["center"] = serial_row["center"]
        pairs.append((serial_row, parallel_row))
    if any(candidates for candidates in parallel_buckets.values()):
        raise AssertionError("unmatched parallel final records")
    return pairs


def latest_time(case):
    def is_time(path):
        try:
            float(path.name)
            return path.is_dir()
        except ValueError:
            return False

    times = [path for path in case.iterdir() if is_time(path)]
    if not times:
        raise AssertionError(f"no output time in {case}")
    return max(times, key=lambda path: float(path.name))


def compare_production(serial_case, parallel_case):
    serial_time = latest_time(serial_case)
    parallel_time = latest_time(parallel_case)
    if serial_time.name != parallel_time.name:
        raise AssertionError("production output times differ")
    for field in PRODUCTION_FIELDS:
        def contents(time_dir):
            path = time_dir/field
            if not path.is_file():
                raise AssertionError(f"missing production field: {path}")
            # OpenFOAM records the generating time directory in this header.
            return [
                line for line in path.read_text().splitlines()
                if "location" not in line
            ]

        if contents(serial_time) != contents(parallel_time):
            raise AssertionError(f"production field differs: {field}")


def compare_records(
    serial_case, parallel_case, serial, parallel, hops, geometry_tolerance
):
    if not serial["owner"]:
        raise AssertionError("serial completed no interpolation points")
    outside = [row for row in serial["owner"] if not row["inside"]]
    if outside:
        raise AssertionError(f"serial has {len(outside)} out-of-cell points")
    if not parallel["owner"]:
        raise AssertionError("parallel completed no interpolation points")
    outside = [row for row in parallel["owner"] if not row["inside"]]
    if outside:
        raise AssertionError(
            f"parallel has {len(outside)} out-of-cell completed points"
        )
    if len(serial["owner"]) != len(parallel["owner"]):
        raise AssertionError("serial/parallel owner completion counts differ")

    receive_counts = Counter(
        (
            row["running"], row["order"],
            row["origin_rank"], row["origin_label"]
        )
        for row in parallel["receive"]
    )
    max_hops = max(receive_counts.values(), default=0)
    if max_hops < hops:
        raise AssertionError(
            f"fixture does not exercise {hops}-hop handoffs "
            f"(maximum observed: {max_hops})"
        )

    pairs = pair_geometry(serial["final"], parallel["final"])
    serial_histogram = histogram([row[0] for row in pairs])
    parallel_histogram = histogram([row[1] for row in pairs])
    if serial_histogram != parallel_histogram:
        raise AssertionError(
            "completed-order histogram differs: "
            + json.dumps(
                {
                    str(key): {
                        "serial": serial_histogram.get(key),
                        "parallel": parallel_histogram.get(key)
                    }
                    for key in sorted(
                        serial_histogram.keys() | parallel_histogram.keys()
                    )
                    if serial_histogram.get(key)
                    != parallel_histogram.get(key)
                },
                sort_keys=True
            )
        )

    error = 0.0
    for serial_row, parallel_row in pairs:
        if serial_row["valid"] != parallel_row["valid"]:
            raise AssertionError(
                f"final validity differs at {serial_row['center']}"
            )
        error = max(
            error,
            max(
                abs(x-y)
                for x, y in zip(serial_row["point"], parallel_row["point"])
            )
        )
    if error > geometry_tolerance:
        raise AssertionError(f"final geometry differs by {error}")

    compare_production(serial_case, parallel_case)
    return {
        "ownerRecords": len(parallel["owner"]),
        "receiveEvents": len(parallel["receive"]),
        "multiHopRecords": sum(n > 1 for n in receive_counts.values()),
        "maxHops": max_hops,
        "surfaceCells": len(serial_histogram),
        "maximumGeometryDifference": error,
        "productionFieldsCompared": len(PRODUCTION_FIELDS)
    }


def replace_once(path, old, new):
    text = path.read_text()
    if text.count(old) != 1:
        raise AssertionError(f"unexpected fixture setting in {path}: {old}")
    path.write_text(text.replace(old, new))


def execute(template, root, name, np, decomp, center, env, deformed=False):
    case = root/f"{name}-np{np}"
    shutil.copytree(template, case)
    decomp_path = case/"system/decomposeParDict"
    replace_once(decomp_path, "numberOfSubdomains 2", f"numberOfSubdomains {np}")
    replace_once(decomp_path, "n (2 1 1)", f"n ({decomp})")
    replace_once(
        case/"constant/HFDIBDEMDict",
        "startPosition (1.3 0.8 1)",
        f"startPosition ({center})"
    )

    run(["blockMesh"], case, "log.blockMesh", env)
    nonorthogonality = None
    if deformed:
        script = Path(__file__).with_name("deform_mesh.py")
        run(
            [
                sys.executable, "-I", str(script), str(case),
                "--seed", "17", "--amplitude", "0.25"
            ],
            case, "log.deformMesh", env
        )
        run(
            ["checkMesh", "-allGeometry", "-allTopology"],
            case, "log.checkMesh", env
        )
        mesh_log = (case/"log.checkMesh").read_text()
        if "Mesh OK." not in mesh_log:
            raise AssertionError(f"deformed mesh is not valid in {case}")
        match = re.search(
            r"non-orthogonality Max:\s*([0-9.eE+-]+)", mesh_log
        )
        if match is None:
            raise AssertionError("no measured nonorthogonality angle")
        nonorthogonality = float(match.group(1))
        if nonorthogonality < 5:
            raise AssertionError("mesh deformation is too small for this test")
    if np == 1:
        run(["pimpleHFDIBFoam"], case, "log.pimpleHFDIBFoam", env)
    else:
        run(["decomposePar", "-force"], case, "log.decomposePar", env)
        run(
            ["mpirun", "-np", str(np), "pimpleHFDIBFoam", "-parallel"],
            case, "log.pimpleHFDIBFoam", env
        )
        run(["reconstructPar", "-latestTime"], case, "log.reconstructPar", env)
    log = (case/"log.pimpleHFDIBFoam").read_text()
    if "End\n" not in log:
        raise RuntimeError(f"solver did not reach End in {case}")
    if "Unresolved marched interpolation points" in log:
        raise AssertionError(f"unexpected unresolved-point warning in {case}")
    return case, read_probes(case), nonorthogonality


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--strict", action="store_true",
        help="abort collectively on any unresolved marched point"
    )
    parser.add_argument(
        "--keep", type=Path,
        help="retain generated cases under this empty/new directory"
    )
    args = parser.parse_args()

    for executable in (
        "blockMesh", "decomposePar", "reconstructPar", "checkMesh", "mpirun",
        "pimpleHFDIBFoam"
    ):
        if shutil.which(executable) is None:
            raise RuntimeError(f"missing executable: {executable}")

    template = Path(__file__).resolve().parent/"case"
    if args.keep:
        args.keep.mkdir(parents=True, exist_ok=False)
        root = args.keep
        cleanup = None
    else:
        cleanup = tempfile.TemporaryDirectory(prefix="lineInt-handoff-")
        root = Path(cleanup.name)

    env = os.environ.copy()
    env["HFDIB_MARCH_DUMP"] = "1"
    if args.strict:
        env["HFDIB_LINEINT_STRICT"] = "1"
    try:
        summaries = {}
        serial = {}
        parallel = {}
        for name, np, decomp, deformed in (
            ("one-hop", 2, "2 1 1", False),
            ("multi-hop", 4, "2 2 1", False),
            ("nonorthogonal-multi-hop", 4, "2 2 1", True),
        ):
            center = "1.3 0.8 1" if np == 2 else "1.55 0.55 1"
            serial[name] = execute(
                template, root, name, 1, "1 1 1", center, env, deformed
            )
            parallel[name] = execute(
                template, root, name, np, decomp, center, env, deformed
            )
            summaries[name] = compare_records(
                serial[name][0], parallel[name][0],
                serial[name][1], parallel[name][1],
                hops=1 if np == 2 else 2,
                geometry_tolerance=1e-8 if deformed else 1e-10
            )
            if deformed:
                for records in (serial[name], parallel[name]):
                    if records[2] != serial[name][2]:
                        raise AssertionError("mesh deformation differs by case")
                summaries[name]["nonOrthogonalityMax"] = serial[name][2]
        print(json.dumps(summaries, indent=2, sort_keys=True))
        print("PASS: handoff coverage, geometry, and production fields agree")
    finally:
        if cleanup is not None:
            cleanup.cleanup()


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, RuntimeError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)
