#!/usr/bin/env python3
"""Regression test: legacy createImmersedBody parallel AMR deadlock.

nonConvexBody::createImmersedBodyLegacy walks an octree whose termination
flag (isInsideBB) is rank-local, while the loop body posts collectives
(PstreamBuffers processor-face exchange, reduce(nextSize, maxOp<label>())).
On decompositions where the body bounding box spans only some ranks, ranks
exit the loop at different iterations and the collectives cross — a
permanent hang at body creation (introduced by 99996e8, fixed by making
the termination global via returnReduceOr).

This test runs the pseudo2DAMR_with_repeatSamePosition case (a hollow
nonConvex body under AMR, the deadlock-triggering decomposition) in
parallel with bodyCreation forced to legacy and a hard timeout. Before
the fix the run deadlocks at the first body creation (0 time steps); the
test fails if no time step completes within the timeout.

A solver binary with the fix completes the (shortened) run normally.
"""
import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
CASE = REPO_ROOT / "tests" / "HFDIBDEMFoam" / "pseudo2DAMR_with_repeatSamePosition"
N_PROCS = 4
# generous for a healthy run (fixed binary completes the shortened case in
# ~15 s wall); short enough to catch the deadlock quickly
TIMEOUT_S = 120
# enough steps to pass several addModel re-creations and AMR refinement
# events; the pre-fix binary hangs at step 0 so even 1 step discriminates
END_TIME = "0.05"


def _foam_env():
    env = os.environ.copy()
    bashrc = "/usr/lib/openfoam/openfoam2412/etc/bashrc"
    if Path(bashrc).exists():
        env["FOAM_BASHRC"] = bashrc
    return env


def _run(cmd, cwd, log_name, timeout, env):
    """Run a command in the case dir, writing output to log.<name>."""
    log = Path(cwd) / f"log.{log_name}"
    with log.open("w") as f:
        return subprocess.run(
            cmd,
            cwd=cwd,
            stdout=f,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            env=env,
        )


def _ignore_case_artifacts(directory, contents):
    """Skip run outputs and bulky artifacts when copying the case."""
    keep_out = re.compile(
        r"^(processor\d+|log\..*|[0-9]+(\..*)?|postProcessing.*|ZZ_.*"
        r"|.*\.foam|.*\.OpenFOAM|.*\.pvsm|gdb_out\.log)$"
    )
    # 0.org / 0.org_AMR / 0.org_2D must be kept (start-time templates)
    keep_out_0 = re.compile(r"^0(?!\.org).*$")
    return [
        c for c in contents
        if (keep_out.match(c) or keep_out_0.match(c)) and c != "0.org_AMR"
    ]


class TestPseudo2DAMRLegacyParallel(unittest.TestCase):
    def test_legacy_creation_parallel_amr_no_deadlock(self):
        for tool in ("blockMesh", "decomposePar", "HFDIBDEMFoam", "mpirun"):
            if shutil.which(tool) is None:
                self.skipTest(f"{tool} not in PATH (source OpenFOAM bashrc)")

        case_src = CASE
        self.assertTrue(
            (case_src / "constant" / "HFDIBDEMDict").exists(),
            f"case missing: {case_src}",
        )

        with tempfile.TemporaryDirectory(prefix="hfdib_legacy_amr_") as tmp:
            case = Path(tmp) / "case"
            shutil.copytree(
                case_src,
                case,
                ignore=_ignore_case_artifacts,
            )

            # fresh start time
            shutil.rmtree(case / "0", ignore_errors=True)
            shutil.copytree(case / "0.org_AMR", case / "0")

            # force legacy creation for every body
            dict_path = case / "constant" / "HFDIBDEMDict"
            dict_text = dict_path.read_text()
            self.assertIn("bodyGeom nonConvex", dict_text)
            dict_text = dict_text.replace(
                "bodyGeom nonConvex",
                "bodyCreation legacy;\n    bodyGeom nonConvex",
            )
            dict_path.write_text(dict_text)

            # shorten the run
            control = case / "system" / "controlDict"
            control_text = control.read_text()
            self.assertIn("endTime", control_text)
            lines = []
            for line in control_text.splitlines():
                if line.strip().startswith("endTime"):
                    lines.append(f"endTime         {END_TIME};")
                else:
                    lines.append(line)
            control.write_text("\n".join(lines) + "\n")

            env = _foam_env()
            try:
                _run(
                    ["blockMesh", "-dict", "system/blockMeshDict_AMR"],
                    case, "blockMesh", 300, env,
                )
                _run(["decomposePar", "-force"], case, "decomposePar", 300, env)
                proc = _run(
                    ["mpirun", "-np", str(N_PROCS), "HFDIBDEMFoam", "-parallel"],
                    case, "run", TIMEOUT_S, env,
                )
            except subprocess.TimeoutExpired:
                self.fail(
                    "HFDIBDEMFoam -parallel timed out: legacy createImmersedBody "
                    "deadlocked (rank-local isInsideBB octree termination in "
                    "nonConvexBody.C)"
                )

            log = (case / "log.run").read_text()
            self.assertEqual(
                proc.returncode, 0,
                f"solver exited with {proc.returncode}; see log tail:\n"
                f"{log[-2000:]}",
            )
            self.assertNotIn("FOAM FATAL", log)
            n_steps = log.count("Time = ")
            self.assertGreater(
                n_steps, 0,
                f"solver completed 0 time steps (deadlock?):\n{log[-2000:]}",
            )
            self.assertIn("End", log, f"run did not reach End:\n{log[-2000:]}")


if __name__ == "__main__":
    unittest.main()
