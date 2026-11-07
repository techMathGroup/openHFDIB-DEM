#!/usr/bin/python3
"""Deform shared interior points of an OpenFOAM blockMesh-generated mesh.

The deformation is adapted from Martin Isoz's nonOrthogonalizeMeshV2.py.
It preserves topology and the outer domain while providing deterministic
nonorthogonal processor interfaces for interpolation tests.

Copyright 2022 Martin Isoz, original nonOrthogonalizeMeshV2.py
Adapted in 2026 for deterministic bounded test-mesh deformation.

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 2 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.
"""

import argparse
from pathlib import Path
import random
import re


def read_points(path):
    text = path.read_text()
    # Accept the plain counted vector list emitted by OpenFOAM blockMesh.
    match = re.search(
        r"(?m)^\s*(\d+)\s*\n\(\s*\n(.*?)^\)", text, re.DOTALL
    )
    if match is None:
        raise ValueError(f"no counted point list in {path}")
    count = int(match.group(1))
    values = re.findall(r"\(([^()]*)\)", match.group(2))
    points = []
    for value in values:
        coordinates = value.split()
        if len(coordinates) != 3:
            raise ValueError(f"invalid point entry: {value}")
        points.append(tuple(float(item) for item in coordinates))
    if len(points) != count:
        raise ValueError(f"expected {count} points, found {len(points)}")
    return text, match.span(2), points


def deform(points, seed, amplitude):
    lower = [min(point[i] for point in points) for i in range(3)]
    upper = [max(point[i] for point in points) for i in range(3)]
    spacing = []
    for i in range(3):
        coordinates = sorted({point[i] for point in points})
        steps = [b-a for a, b in zip(coordinates, coordinates[1:]) if b > a]
        if not steps:
            raise ValueError("the mesh has no extent in every coordinate")
        spacing.append(min(steps))

    generator = random.Random(seed)
    result = []
    deformed = 0
    for point in points:
        changed = list(point)
        for i in (0, 1):
            if lower[i] < changed[i] < upper[i]:
                changed[i] += generator.uniform(-amplitude, amplitude)*spacing[i]
                deformed += 1
        result.append(tuple(changed))
    return result, deformed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("case", type=Path)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument(
        "--amplitude", type=float, default=0.20,
        help="maximum coordinate displacement per grid spacing and axis"
    )
    args = parser.parse_args()
    if not 0 <= args.amplitude <= 0.35:
        raise ValueError("amplitude must be in [0, 0.35]")

    path = args.case/"constant/polyMesh/points"
    text, value_span, points = read_points(path)
    deformed, count = deform(points, args.seed, args.amplitude)
    values = "\n".join(
        "(" + " ".join(format(x, ".17g") for x in point) + ")"
        for point in deformed
    )
    path.write_text(text[:value_span[0]] + values + "\n" + text[value_span[1]:])
    print(f"deformed {count} point coordinates with seed {args.seed}")


if __name__ == "__main__":
    main()
