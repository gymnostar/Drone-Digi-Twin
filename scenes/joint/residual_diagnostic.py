#!/usr/bin/env python3
"""Summarize OpenFOAM v2306 solverInfo fields from an instrumented checkpoint replay."""
import argparse
from collections import Counter
from contextlib import nullcontext
from heapq import nlargest
import json
from math import dist, floor, fsum, isclose
from pathlib import Path
import re
from shutil import copy2, copytree
from statistics import median
from subprocess import TimeoutExpired, run
from tempfile import TemporaryDirectory

from study import city, runner


def replay(source, case, checkpoint, image):
    """Replay three outer iterations in an isolated copy of the archived case."""
    for name in ("0", "constant", "system"):
        copytree(source / name, case / name)
    copy2(source / "parameters.json", case)
    for rank in range(3):
        part = f"processor{rank}"
        for name in ("constant", str(checkpoint)):
            copytree(source / part / name, case / part / name)

    path = case / "system/controlDict"
    control = path.read_text().split("functions {", 1)[0]
    for key, value in {
        "startTime": checkpoint, "endTime": checkpoint + 3,
        "writeInterval": 1, "purgeWrite": 0,
    }.items():
        control, count = re.subn(rf"\b{key}\s+[^;]+;", f"{key} {value};", control)
        if count != 1:
            raise ValueError(f"Expected one archived controlDict entry for {key}")
    path.write_text(control + """functions {
    localResiduals {
        type solverInfo;
        libs (utilityFunctionObjects);
        fields (U p);
        writeResidualFields true;
        writeControl timeStep;
        writeInterval 1;
    }
}
""")
    first, last = checkpoint + 2, checkpoint + 3
    command = f"""source /usr/lib/openfoam/openfoam2306/etc/bashrc &&
        mpirun --allow-run-as-root -np 3 simpleFoam -parallel \
            -debug-switch lduMatrix=2 -debug-switch GAMG=2 &&
        reconstructPar -time {first}:{last} -fields "(initialResidual:Ux initialResidual:p)" &&
        postProcess -func writeCellCentres -time {first} &&
        checkMesh -time {first} -writeAllFields"""
    log = case / "log.residual-norms"
    container = f"niigata-nne-residual-{case.name}"
    with log.open("w") as output:
        try:
            process = run([runner.DOCKER, "run", "--rm", "--name", container,
                           "--cpus", "3", "--memory", "5g", "--network", "none",
                           "-v", f"{case}:/case", "-w", "/case", "--entrypoint", "bash",
                           image, "-lc", command], stdout=output, stderr=output, timeout=120)
        except TimeoutExpired:
            run([runner.DOCKER, "rm", "-f", container], capture_output=True, timeout=20)
            raise RuntimeError(f"OpenFOAM replay exceeded 120 seconds. Last log output:\n{log.read_text()[-3000:]}") from None
    if process.returncode:
        raise RuntimeError(f"OpenFOAM replay exited with code {process.returncode}. Last log output:\n{log.read_text()[-3000:]}")


def read_list(text, element):
    count, body = re.search(r"(?m)^\s*(\d+)\s*\n\(\n(.*?)^\)", text, re.S).groups()
    result = [element(row) for row in body.splitlines()]
    if len(result) != int(count):
        raise ValueError("OpenFOAM list length mismatch")
    return result


def vector(row):
    return tuple(map(float, row.strip("()").split()))


def patch_info(mesh):
    text = (mesh / "boundary").read_text()
    owners = read_list((mesh / "owner").read_text(), int)
    return {
        name: set(owners[int(start):int(start) + int(size)])
        for name, size, start in re.findall(
            r"(?m)^\s*(\w+)\s*\{[^}]*nFaces\s+(\d+);[^}]*startFace\s+(\d+);", text
        )
    }


def iteration_log(log, iteration):
    match = re.search(rf"(?ms)^Time = {iteration}\s*\n(.*?)(?=^Time = |\Z)", log)
    if match is None:
        raise ValueError(f"Missing iteration {iteration} in solver log. Choose a checkpoint with three subsequent solver iterations.")
    return match[1]


def first_solve_norms(block, field):
    return tuple(map(float, re.search(
        rf"(?m)^(?:GAMG|smoothSolver):\s+Solving for {field}, "
        r"Initial residual = ([^,]+), Final residual = ([^,]+)", block
    ).groups()))


def normalization(block, field):
    if field == "p":
        factor = re.search(r"(?m)^\[0\]\s+Normalisation factor = (\S+)", block)[1]
    else:
        factor = re.search(
            rf"Normalisation factor = (\S+)\n[^\n]*Solving for {field},", block
        )[1]
    return float(factor)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Archived case directory")
    parser.add_argument("--checkpoint", required=True, type=int)
    parser.add_argument("--reuse", type=Path, help="Summarize a previously instrumented replay")
    parser.add_argument("--geometry-time", type=int, help="Geometry fields time for an existing replay")
    args = parser.parse_args()
    if args.checkpoint <= 0 or (args.geometry_time is not None and args.reuse is None):
        parser.error("Supply a positive checkpoint. --geometry-time requires --reuse.")

    params = json.loads((args.source / "parameters.json").read_text())
    context = nullcontext(args.reuse) if args.reuse else TemporaryDirectory(prefix="nne-residual-")
    with context as directory:
        case = Path(directory).resolve()
        if args.reuse is None:
            replay(args.source, case, args.checkpoint, params["image"])
        if json.loads((case / "parameters.json").read_text()) != params:
            raise ValueError("Replay parameters differ from the archived case")
        print(json.dumps(summarize(args, case, params), indent=2, allow_nan=False))


def summarize(args, case, params):
    geometry_time = args.geometry_time or args.checkpoint + 2
    geometry = case / str(geometry_time)
    angle = params["wind_from_deg"]
    hub = params["rotor_center_solver_m"]
    log = (case / "log.residual-norms").read_text()
    archived_log = (args.source / "log.simpleFoam").read_text()
    mesh = case / "constant/polyMesh"
    patches = patch_info(mesh)
    centers_text = (geometry / "C").read_text()
    centers = read_list(centers_text, vector)
    building_centers = read_list(centers_text.split("    buildings\n")[1], vector)
    rotor = set(read_list((mesh / "sets/rotor").read_text(), int))
    quality = {
        name: read_list((geometry / name).read_text(), float)
        for name in ("nonOrthoAngle", "skewness", "cellDeterminant", "aspectRatio")
    }
    regions = patches | {"rotor_source_cells": rotor}
    all_ids = range(len(centers))
    unassigned = set(all_ids).difference(*regions.values())
    regions["no_boundary_face_and_outside_rotor"] = unassigned
    regions["building_and_ground_faces"] = patches["buildings"] & patches["ground"]
    world = [city.transform(point, angle, inverse=True) for point in centers]
    bins = [tuple(floor(coordinate / 25) * 25 for coordinate in point) for point in world]
    report = {
        "source_case": args.source.name, "checkpoint": args.checkpoint,
        "coupling": "SIMPLE" if params.get("standard_simple_recovery", False) else "SIMPLEC",
        "wind_from_deg": params["wind_from_deg"],
        "image": params["image"], "openfoam": "v2306", "cells": len(centers),
        "residual_convention": "Signed unnormalized b-A*psi on the first linear solve of each outer "
        "iteration, ranked by absolute value. Pressure uses the first non-orthogonal solve. "
        "The logged initial norm is the global L1 divided by the solver's normalization factor.",
        "coordinate_frames": "World x east, y north, z up in metres. Solver x crosswind, y upwind, z up. "
        "Building distance is to the nearest boundary face centre. Patch owners can overlap.",
        "rotor_hub_world_m": params["rotor_center_m"], "rotor_diameter_m": params["rotor_D_m"],
        "geometry_time": geometry_time,
        "quality_mesh_medians": {name: median(values) for name, values in quality.items()},
        "iterations": {},
    }
    # v2306 solverInfo creates residual storage after the first replay solve.
    for iteration in (args.checkpoint + 2, args.checkpoint + 3):
        block = iteration_log(log, iteration)
        archived_block = iteration_log(archived_log, iteration)
        fields = {}
        for field in ("Ux", "p"):
            values = read_list((case / str(iteration) / f"initialResidual:{field}").read_text(), float)
            if len(values) != len(centers):
                raise ValueError(f"Cell count mismatch for {field} at {iteration}")
            absolute = list(map(abs, values))
            total = fsum(absolute)
            initial, final = first_solve_norms(block, field)
            archived_initial, _ = first_solve_norms(archived_block, field)
            factor = normalization(block, field)
            reconstructed = total / factor
            if not isclose(reconstructed, initial, rel_tol=1e-7):
                raise ValueError(f"Field/log initial norm mismatch for {field} at {iteration}")
            top = nlargest(100, all_ids, key=absolute.__getitem__)
            binned = Counter()
            for cell, lower_corner in enumerate(bins):
                binned[lower_corner] += absolute[cell]
            fields[field] = {
                "raw_L1": total, "normalization_factor": factor, "logged_initial": initial,
                "field_reconstructed_initial": reconstructed, "logged_final": final,
                "archived_initial": archived_initial,
                "top_100_L1_fraction": fsum(absolute[cell] for cell in top)/total,
                "regions": {
                    name: {"mesh_cells": len(cells),
                           "L1_fraction": fsum(absolute[cell] for cell in cells)/total,
                           "top_100_cells": sum(cell in cells for cell in top)}
                    for name, cells in regions.items()
                },
                "world_25m_bins_by_L1": [
                    {"lower_corner_m": corner, "L1_fraction": value/total}
                    for corner, value in binned.most_common(5)
                ],
                "top_100_quality_medians": {
                    name: median(values[cell] for cell in top) for name, values in quality.items()
                },
                "top_5_cells": [
                    {
                        "cell": cell, "signed_residual": values[cell], "solver_m": centers[cell],
                        "world_m": world[cell], "building_face_center_distance_m":
                            min(dist(centers[cell], face) for face in building_centers),
                        "rotor_hub_distance_m": dist(centers[cell], hub),
                        "regions": [name for name, cells in regions.items() if cell in cells],
                        "mesh_quality": {name: values[cell] for name, values in quality.items()},
                    }
                    for cell in top[:5]
                ],
            }
        report["iterations"][iteration] = fields
    return report


if __name__ == "__main__":
    main()
