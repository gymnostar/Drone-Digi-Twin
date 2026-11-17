#!/usr/bin/env python3
"""Morning script: read the 3 mesh-study cases and write Arthur's tables + 2 plots.
Run from the repo root:   .venv/bin/python DLR/meshStudy/make_table.py
Works with partial results (e.g. fine still running): it uses whatever exists.

Output: DLR/results/mesh_study.md and .html, mesh_study_convergence.png, mesh_study_profiles.png
"""
import csv, html, math, re
from pathlib import Path

import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "results"
CASES = ["coarse", "medium", "fine"]
BOX_CELL = {"coarse": 4, "medium": 2, "fine": 1}   # m, inside the cluster box
RHO = 1.2                                           # kg/m3: simpleFoam's p is kinematic (p/rho)
R = 2.0                                             # refinement ratio between meshes
IDX = list(csv.DictReader(open(HERE / "probes_index.csv")))

def num(pattern, text, cast=float, last=True):
    m = re.findall(pattern, text)
    return cast(m[-1] if last else m[0]) if m else None

def mesh_quality(case):
    t = (HERE / case / "log.checkMesh").read_text()
    return dict(
        cells=num(r"^\s+cells:\s+(\d+)", t.replace("\r", ""), int) or num(r"cells:\s+(\d+)", t, int),
        skew=num(r"Max skewness = ([\d.eE+-]+)", t),
        skew_bad=num(r"(\d+) highly skew faces", t, int) or 0,
        northo_max=num(r"non-orthogonality Max: ([\d.eE+-]+)", t),
        northo_avg=num(r"average: ([\d.eE+-]+)", t),
        aspect=num(r"Max aspect ratio = ([\d.eE+-]+)", t),
        failed=num(r"Failed (\d+) mesh checks", t, int) or 0,
        neg_vol=num(r"Zero or negative cell volume detected.*?(\d+)", t, int) or 0,
    )

def convergence(case):
    t = (HERE / case / "log.simpleFoam").read_text()
    last = t[t.rfind("\nTime = "):]
    res = {f: num(rf"Solving for {f}, Initial residual = ([\d.eE+-]+)", last) for f in ("Ux", "Uy", "Uz", "p", "k", "epsilon")}
    return dict(iters=num(r"\nTime = (\d+)", t, int), res=res,
                cont=num(r"sum local = ([\d.eE+-]+)", last), clock=num(r"ClockTime = ([\d.]+)", t))

def read_probe(case, field):
    """Returns (iterations, values[iter, probe(, 3)]) and the probe indices found by OpenFOAM."""
    f = HERE / case / "postProcessing" / "probes" / "0" / field
    lines = f.read_text().splitlines()
    found = [int(m.group(1)) for l in lines if (m := re.match(r"#\s*Probe\s+(\d+)", l))]
    rows = [l for l in lines if l and not l.startswith("#")]
    it = np.array([float(r.split()[0]) for r in rows])
    if field == "U":
        vals = np.array([[float(v) for v in re.findall(r"[-\d.eE+]+", r.split(None, 1)[1])] for r in rows])
        vals = vals.reshape(len(rows), -1, 3)
    else:
        vals = np.array([[float(v) for v in r.split()[1:]] for r in rows])
    return it, vals, found

data = {}
for c in CASES:
    d = HERE / c
    if not (d / "log.checkMesh").exists():
        continue
    data[c] = {"mesh": mesh_quality(c)}
    if (d / "postProcessing/probes/0/U").exists():
        it, U, found = read_probe(c, "U")
        _, P, _ = read_probe(c, "p")
        spd = np.linalg.norm(U, axis=2)
        data[c].update(conv=convergence(c), it=it, spd=spd, p=P, found=found)

def probe(c, loc, h, what="spd"):
    """final value at location/height; None if not available."""
    if "found" not in data.get(c, {}):
        return None
    for r in IDX:
        if r["location"] == loc and float(r["height_above_ground"]) == h:
            i = int(r["i"])
            if i not in data[c]["found"]:
                return None
            return float(data[c][what][-1, data[c]["found"].index(i)])

def drift(c, loc, h):
    """% change of the probe speed over the last 100 iterations (is it still moving?)."""
    if "found" not in data.get(c, {}):
        return None
    i = next(int(r["i"]) for r in IDX if r["location"] == loc and float(r["height_above_ground"]) == h)
    if i not in data[c]["found"]:
        return None
    s = data[c]["spd"][:, data[c]["found"].index(i)]
    w = s[-100:]
    return 100 * (w.max() - w.min()) / max(abs(s[-1]), 1e-9)

def gci(f3, f2, f1):
    """Roache's 3-grid GCI on the fine grid (f1 = fine). Returns (observed order p, GCI %) or a reason."""
    if None in (f1, f2, f3):
        return None, "needs 3 meshes"
    e32, e21 = f3 - f2, f2 - f1
    if abs(e21) < 1e-12:
        return None, "0 % (identical)"
    if e32 * e21 <= 0:
        return None, "oscillating (not monotonic)"
    p = math.log(abs(e32 / e21)) / math.log(R)
    p_used = min(max(p, 0.5), 2.0)        # common safeguard: clip p to [0.5, 2]
    return p, f"{100 * 1.25 * abs(e21 / f1) / (R ** p_used - 1):.1f} %"

fmt = lambda v, f="{:.2f}": "–" if v is None else f.format(v)
md = ["# DLR mesh-independence study", "",
      "West wind 5 m/s @ 10 m, steady RANS k-ε (simpleFoam). The 3 meshes differ **only** in the cell size inside a "
      "185 × 117 × 48 m box around the DLR weather-station building (4 → 2 → 1 m, ratio 2).", ""]

md += ["## 1. Mesh quality", "",
       "| Mesh | Box cells | Total cells | Max skewness (bad faces) | Non-ortho max / avg | Max aspect ratio | Failed checks |",
       "|---|---|---|---|---|---|---|"]
for c, d in data.items():
    m = d["mesh"]
    md.append(f"| {c} | {BOX_CELL[c]} m | {m['cells']:,} | {fmt(m['skew'])} ({m['skew_bad']}) | "
              f"{fmt(m['northo_max'], '{:.1f}')}° / {fmt(m['northo_avg'], '{:.1f}')}° | {fmt(m['aspect'], '{:.1f}')} | {m['failed']} |")
md += ["", "Usual limits: skewness < 4, non-orthogonality < 70°, no negative volumes. "
       "`checkMesh -writeSets` saved the bad faces to `<case>/postProcessing` / `constant/polyMesh/sets` so you can see where they are.", ""]

md += ["## 2. Convergence", "",
       "| Mesh | Iterations | Final residuals Ux / p / k | Continuity error (sum local) | Max probe drift, last 100 it | Wall time |",
       "|---|---|---|---|---|---|"]
for c, d in data.items():
    if "conv" not in d:
        md.append(f"| {c} | not finished | | | | |"); continue
    v = d["conv"]; r = v["res"]
    drifts = [drift(c, l["location"], float(l["height_above_ground"])) for l in IDX]
    md.append(f"| {c} | {v['iters']} | {fmt(r['Ux'], '{:.0e}')} / {fmt(r['p'], '{:.0e}')} / {fmt(r['k'], '{:.0e}')} | "
              f"{fmt(v['cont'], '{:.0e}')} | {fmt(max([x for x in drifts if x is not None], default=None), '{:.2f}')} % | "
              f"{fmt(v['clock'] / 3600 if v['clock'] else None, '{:.1f}')} h |")
md += ["", "Residuals say how unbalanced the equations still are; the **probe drift** says whether the answer still moves. "
       "Under ~0.5 % = converged.", ""]

KEY = [("upstream", 10), ("gap", 5), ("gap", 20), ("wake", 5), ("wake", 20), ("station", 20), ("station", 40)]
md += ["## 3. Mesh independence (wind speed |U| in m/s at the probes)", "",
       "| Probe (height above ground) | coarse | medium | fine | coarse→medium | medium→fine | GCI fine |",
       "|---|---|---|---|---|---|---|"]
for loc, h in KEY:
    f3, f2, f1 = (probe(c, loc, h) for c in CASES)
    ch = lambda a, b: "–" if None in (a, b) else f"{100 * (b - a) / a:+.1f} %"
    _, g = gci(f3, f2, f1)
    md.append(f"| {loc} {h} m | {fmt(f3)} | {fmt(f2)} | {fmt(f1)} | {ch(f3, f2)} | {ch(f2, f1)} | {g} |")
dp = {c: (None if probe(c, "upstream", 5, "p") is None or probe(c, "wake", 5, "p") is None
          else RHO * (probe(c, "upstream", 5, "p") - probe(c, "wake", 5, "p"))) for c in CASES}
_, g = gci(dp["coarse"], dp["medium"], dp["fine"])
md.append(f"| **Δp upstream→wake, 5 m (Pa)** | {fmt(dp['coarse'])} | {fmt(dp['medium'])} | {fmt(dp['fine'])} | | | {g} |")
md += ["", "**How to read it:** the upstream probe is the control (should be ≈ identical). "
       "If medium→fine is below ~2–5 % at the other probes, the medium mesh is mesh-independent and is good enough. "
       "GCI = Grid Convergence Index (Roache): an error band on the fine-mesh value.", ""]

# plots
OUT.mkdir(exist_ok=True)
locs = ["upstream", "gap", "wake", "station"]
cols = {"coarse": "#9ca3af", "medium": "#2563eb", "fine": "#dc2626"}
fig, ax = plt.subplots(1, 4, figsize=(15, 3.6), sharey=False)
for a, loc in zip(ax, locs):
    h = 20 if loc == "station" else 5
    i = next(int(r["i"]) for r in IDX if r["location"] == loc and float(r["height_above_ground"]) == h)
    for c, d in data.items():
        if "found" in d and i in d["found"]:
            a.plot(d["it"], d["spd"][:, d["found"].index(i)], color=cols[c], label=c)
    a.set_title(f"{loc}, {h} m above ground"); a.set_xlabel("iteration"); a.grid(alpha=.3)
ax[0].set_ylabel("|U| (m/s)"); ax[0].legend()
fig.suptitle("Convergence: probe speed vs iteration (flat at the end = converged)")
fig.tight_layout(); fig.savefig(OUT / "mesh_study_convergence.png", dpi=110)

fig, ax = plt.subplots(1, 4, figsize=(15, 4.2), sharex=True)
for a, loc in zip(ax, locs):
    rows = [r for r in IDX if r["location"] == loc]
    for c, d in data.items():
        if "found" not in d:
            continue
        hs = [float(r["height_above_ground"]) for r in rows if int(r["i"]) in d["found"]]
        vs = [d["spd"][-1, d["found"].index(int(r["i"]))] for r in rows if int(r["i"]) in d["found"]]
        a.plot(vs, hs, "o-", color=cols[c], label=c)
    a.set_title(loc); a.set_xlabel("|U| (m/s)"); a.grid(alpha=.3)
ax[0].set_ylabel("height above ground (m)"); ax[0].legend()
fig.suptitle("Vertical wind profiles at the 4 probe locations, per mesh")
fig.tight_layout(); fig.savefig(OUT / "mesh_study_profiles.png", dpi=110)

md += ["![convergence](mesh_study_convergence.png)", "", "![profiles](mesh_study_profiles.png)", ""]
(OUT / "mesh_study.md").write_text("\n".join(md))

# minimal HTML version (tables + images) for the meeting
def md_table_to_html(lines):
    out, rows = [], [l for l in lines if l.startswith("|") and not l.startswith("|---")]
    for k, l in enumerate(rows):
        cells = [html.escape(x.strip()).replace("**", "") for x in l.strip("|").split("|")]
        tag = "th" if k == 0 else "td"
        out.append("<tr>" + "".join(f"<{tag}>{x}</{tag}>" for x in cells) + "</tr>")
    return "<table>" + "".join(out) + "</table>"
body, buf = [], []
for l in md + [""]:
    if l.startswith("|"):
        buf.append(l); continue
    if buf:
        body.append(md_table_to_html(buf)); buf = []
    if l.startswith("# "): body.append(f"<h1>{html.escape(l[2:])}</h1>")
    elif l.startswith("## "): body.append(f"<h2>{html.escape(l[3:])}</h2>")
    elif l.startswith("!["): body.append(f'<img src="{re.search(r"\((.*?)\)", l).group(1)}">')
    elif l: body.append(f"<p>{html.escape(l).replace('**', '')}</p>")
(OUT / "mesh_study.html").write_text(
    "<!doctype html><meta charset=utf-8><title>DLR mesh study</title><style>"
    "body{font:15px/1.5 -apple-system,sans-serif;max-width:1100px;margin:30px auto;padding:0 20px}"
    "table{border-collapse:collapse;margin:8px 0}td,th{border:1px solid #ccd;padding:5px 9px}th{background:#eef1f8}"
    "img{max-width:100%}</style>" + "\n".join(body))
print("\n".join(md))
print(f"\nwrote {OUT / 'mesh_study.html'}")
