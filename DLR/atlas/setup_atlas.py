#!/usr/bin/env python3
"""Build the DLR wind-atlas cases: ONE mesh, many wind directions.

Why one mesh: the mesh study (DLR/meshStudy) showed the base mesh ("coarse": 16 m background,
4 m over campus, 2 m at walls) is within ~1-10 % of the 1 m mesh at drone heights (>= 20 m AGL).
The domain is a cylinder, so a new wind direction only needs a new inlet/outlet split of the
cylinder side and a new inflow direction. No re-meshing.

Why one speed per direction: in neutral flow around sharp-edged buildings the flow pattern
barely depends on speed (Reynolds independence), so U(x)/U_ref is the same at 5 and 10 m/s and
the field can be rescaled to the measured speed. The atlas includes ONE 10 m/s run to test this.

Conventions
    direction theta = meteorological "FROM" angle, degrees clockwise from north (270 = west wind)
    case coordinates: x = east, y = north, z = up (m); case = scene_local - (-335, 279)
    flowDir (where the air goes) = (-sin theta, -cos theta, 0)

Usage: called by run_atlas.sh (`setup_atlas.py base`, then `setup_atlas.py`).
Creates DLR/atlas/cases.noindex/base (mesh with an unsplit cylinder side) and one folder per run.
".noindex" stops Spotlight from indexing hundreds of MB of solver output.
"""
import math, re, shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "meshStudy" / "coarse"        # the base mesh, already meshed and checked
CASES = HERE / "cases.noindex"
NITER = 450      # validation runs and the first directions; later directions were cut to 350 at 01:40 (drift at 240 deg: 0.09 % of Uref per 50 it at drone heights)
NPROC = 8

# Priority order = EDMO 2020-2025 wind climate (share of hours / share of wind energy sum U^3):
# 240: 23 %/58 %   210: 16 %/7 %   60: 8 %/5 %   180: 8 %/2 %   30: 7 %/3 %   150: 7 %/1 %
# 90: 6 %/3 %   300: 4 %/3 %   0: 4 %/1 %   330: 4 %/2 %   120: 3 %/0.4 %   (270: 10 %/15 %, already
# solved in the mesh study). The 10 m/s run tests Reynolds independence on the dominant sector.
RUNS = [(240, 5), (210, 5),
        (255, 5), (240, 10),   # validation runs first: they only need 240 + 270, which are solved first
        (60, 5), (180, 5), (30, 5), (150, 5), (90, 5), (300, 5), (0, 5), (330, 5), (120, 5)]
# 255 = VALIDATION ONLY: tests the real 30 deg blend (240 + 270) where the wind blows most often.


def name(theta, u):
    return f"d{theta:03d}_u{u:02d}"


def flow_dir(theta):
    t = math.radians(theta)
    v = (-math.sin(t), -math.cos(t), 0.0)
    return tuple(0.0 if abs(c) < 1e-12 else round(c, 10) for c in v)


def build_base():
    """Copy the mesh-study mesh. run_atlas.sh then runs createPatchDict.merge and touches BASE_OK."""
    base = CASES / "base"
    if (base / "BASE_OK").exists():
        return base
    shutil.rmtree(base, ignore_errors=True)
    (base / "constant").mkdir(parents=True)
    shutil.copytree(SRC / "system", base / "system")
    shutil.copytree(SRC / "0.orig", base / "0.orig")
    shutil.copytree(SRC / "constant" / "polyMesh", base / "constant" / "polyMesh",
                    ignore=shutil.ignore_patterns("sets"))
    for f in ("transportProperties", "turbulenceProperties"):
        shutil.copy(SRC / "constant" / f, base / "constant" / f)
    # merge the mesh study's west-wind split back into one cylinder side
    (base / "system" / "createPatchDict.merge").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary; object createPatchDict; }\n"
        "pointSync false;\n"
        "patches ( { name side; patchInfo { type patch; } constructFrom patches; patches (inlet outlet); } );\n")
    return base


def make_run(base, theta, u):
    assert (base / "BASE_OK").exists(), "merge the base mesh first (run_atlas.sh does it)"
    c = CASES / name(theta, u)
    if (c / "DONE").exists():
        return c, False
    shutil.rmtree(c, ignore_errors=True)
    shutil.copytree(base, c, ignore=shutil.ignore_patterns("BASE_OK", "log.*"))
    fx, fy, _ = flow_dir(theta)
    up = (-fx, -fy, 0.0)                     # outward normal of the faces that look INTO the wind

    abl = c / "0.orig/include/ABLConditions"
    t = abl.read_text()
    t = re.sub(r"Uref\s+[^;]+;", f"Uref     {float(u)};", t)
    t = re.sub(r"flowDir\s+\([^)]*\);.*", f"flowDir  ({fx} {fy} 0);   // wind FROM {theta} deg", t)
    abl.write_text(t)
    ic = c / "0.orig/include/initialConditions"
    ic.write_text(re.sub(r"flowVelocity\s+\([^)]*\);", f"flowVelocity         ({u*fx:.6g} {u*fy:.6g} 0);",
                         ic.read_text()))

    ts = c / "system/topoSetDict.inletOutlet"
    t = ts.read_text()
    t = re.sub(r"normal \([^)]*\);", f"normal ({up[0]} {up[1]} 0);", t)
    t = t.replace("west, -x", f"from {theta} deg")
    ts.write_text(t)

    cd = c / "system/controlDict"
    t = cd.read_text()
    t = re.sub(r"endTime\s+\d+;", f"endTime         {NITER};", t)
    t = re.sub(r"writeInterval\s+\d+;", f"writeInterval   {NITER};", t, count=1)   # first = solver; the probes keep 1
    t = re.sub(r"startFrom\s+\w+;", "startFrom       startTime;", t)
    cd.write_text(t)
    dp = c / "system/decomposeParDict"
    dp.write_text(re.sub(r"numberOfSubdomains\s+\d+;", f"numberOfSubdomains {NPROC};", dp.read_text()))
    (c / "case.foam").touch()
    (c / "RUN_INFO").write_text(f"theta_from_deg {theta}\nUref_10m {u}\nflowDir {fx} {fy} 0\n")
    return c, True


if __name__ == "__main__":
    import sys
    base = build_base()
    if sys.argv[1:] == ["base"]:
        print("base:", base); sys.exit(0)
    for th, u in RUNS:
        c, new = make_run(base, th, u)
        print(("created " if new else "done    ") + c.name, "flowDir", flow_dir(th))
