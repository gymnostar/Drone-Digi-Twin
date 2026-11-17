#!/usr/bin/env python3
"""Build the 3 mesh-study cases (coarse / medium / fine) from ../windCase.

WHAT CHANGES BETWEEN THE 3 CASES: only the cell size inside one box around the
DLR weather-station building ("cluster box"). Everything else (domain, inlet,
physics, solver, iterations, probes) is identical, so any difference in the
probe values comes from the mesh alone.

    coarse : box level 2 = 4 m cells (no extra refinement; same as the original run)
    medium : box level 3 = 2 m cells
    fine   : box level 4 = 1 m cells  (drone scale, Arthur's request)

Each level halves the cell size -> refinement ratio r = 2 (needed for the GCI).
Coordinates are CASE coordinates (m): scene_local - (-335, 279). See site_map.py.
Run from the repo root:  .venv/bin/python DLR/meshStudy/setup.py
"""
import re, shutil
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "windCase"
LEVELS = {"coarse": 2, "medium": 3, "fine": 4}

# Cluster box: station building + its eastern neighbour + the gap + the wake (x, y, z in m)
BOX_MIN = (50.0, 145.0, 12.0)
BOX_MAX = (235.0, 262.0, 60.0)

# Probe locations: Arthur's 4 regions. Each one is a vertical line (heights above local ground).
# Wind blows from the west (270 deg), i.e. towards +x.
LOCATIONS = {
    "upstream": ((-150.0, 206.0), [2, 5, 10, 20, 40]),  # open field west of campus: CONTROL, should not change
    "gap":      ((172.0, 182.0), [2, 5, 10, 20, 40]),   # between station building and the tall block east of it
    "wake":     ((140.0, 210.0), [2, 5, 10, 20, 40]),   # ~15-20 m behind (east of) the station building
    "station":  ((122.0, 206.0), [15, 20, 25, 40]),     # DLR weather station xy (roof ~12 m): "above campus"
}
NITER = 1000
NPROC = 8          # i9-9880H: 8 physical cores (hyperthreads don't speed up CFD)

site = np.load(HERE / "site_map.npz")
def ground(px, py):
    return float(site["hg"][np.abs(site["y"] - py).argmin(), np.abs(site["x"] - px).argmin()])

points, index = [], ["i,location,x,y,height_above_ground,z"]
for name, ((px, py), heights) in LOCATIONS.items():
    g = ground(px, py)
    for h in heights:
        index.append(f"{len(points)},{name},{px},{py},{h},{g + h:.2f}")
        points.append(f"({px} {py} {g + h:.2f})")
(HERE / "probes_index.csv").write_text("\n".join(index) + "\n")

FUNCTIONS = f"""
functions
{{
    probes   // virtual anemometers: value at fixed points, every iteration
    {{
        type            probes;
        libs            (sampling);
        writeControl    timeStep;
        writeInterval   1;
        fields          (U p k);
        fixedLocations  true;
        probeLocations
        (
            {chr(10).join("            " + p for p in points).lstrip()}
        );
    }}
    residuals   // residual history, for the convergence column
    {{
        type            solverInfo;
        libs            (utilityFunctionObjects);
        fields          (U p k epsilon);
    }}
}}
"""

for case, lev in LEVELS.items():
    dst = HERE / case
    if dst.exists():
        shutil.rmtree(dst)
    (dst / "constant").mkdir(parents=True)
    shutil.copytree(SRC / "0.orig", dst / "0.orig")
    shutil.copytree(SRC / "system", dst / "system")
    for f in ("transportProperties", "turbulenceProperties"):
        shutil.copy(SRC / "constant" / f, dst / "constant" / f)
    shutil.copytree(SRC / "constant" / "triSurface", dst / "constant" / "triSurface")
    s = dst / "system"

    # caseParams: add the cluster box and its level
    with open(s / "caseParams", "a") as fh:
        fh.write(f"\n// mesh study ({case})\nclusterMin ({BOX_MIN[0]} {BOX_MIN[1]} {BOX_MIN[2]});\n"
                 f"clusterMax ({BOX_MAX[0]} {BOX_MAX[1]} {BOX_MAX[2]});\nclusterLevel {lev};\n")

    # snappyHexMeshDict: declare the box and refine inside it
    t = (s / "snappyHexMeshDict").read_text()
    t = t.replace("    turbineRefine { type searchableSphere;",
                  "    clusterRefine { type searchableBox; min $clusterMin; max $clusterMax; }\n"
                  "    turbineRefine { type searchableSphere;", 1)
    t = t.replace("        turbineRefine { mode inside;",
                  "        clusterRefine { mode inside; levels ((1e15 $clusterLevel)); }\n"
                  "        turbineRefine { mode inside;", 1)
    assert t.count("clusterRefine") == 2, "snappyHexMeshDict layout changed"
    (s / "snappyHexMeshDict").write_text(t)

    # controlDict: same iteration count for all, write only at the end (disk is nearly full), add probes
    t = (s / "controlDict").read_text()
    t = re.sub(r"endTime\s+\d+;", f"endTime         {NITER};", t)
    t = re.sub(r"writeInterval\s+\d+;", f"writeInterval   {NITER};", t)
    (s / "controlDict").write_text(t.rstrip() + "\n" + FUNCTIONS)

    t = (s / "decomposeParDict").read_text()
    (s / "decomposeParDict").write_text(re.sub(r"numberOfSubdomains\s+\d+;", f"numberOfSubdomains {NPROC};", t))

    shutil.copy(HERE / "Allrun.case", dst / "Allrun")
    (dst / "Allrun").chmod(0o755)
    print(f"{case:7s}: cluster box level {lev} ({16 / 2**lev:g} m cells)")
print(f"{len(points)} probes -> probes_index.csv")
