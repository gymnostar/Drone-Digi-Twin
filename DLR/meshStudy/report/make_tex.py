#!/usr/bin/env python3
r"""LaTeX tables and data-driven sentences for the mesh-study report.
Run:  .venv/bin/python DLR/meshStudy/report/make_tex.py  [--include-running]

Reads the same logs and probe files as ../make_table.py. Its parsing logic is COPIED here (with a few
robustness fixes) because make_table.py does all its work at import time. Works with partial results:
a case that is missing or not finished shows a dash.
  --include-running   also show the latest probe values of runs that have not finished (marked with a dagger)

Writes, in this folder:
    tables.tex    the three tables (mesh quality, convergence, mesh independence) as booktabs floats
    tables_appendix.tex   all 19 probes;  tab_quality/convergence/independence.tex: bare tabulars for the slides
    summary.tex   macros used by the .tex files: status line, cell counts, findings sentences.
                  Every number in the report text comes from here, never typed by hand.
    fig/          copies of the pictures (so this folder can be zipped and compiled on Overleaf)
                  + residuals.png (residual and continuity-error history, made here)
"""
import base64, csv, math, re, shutil, sys
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent          # meshStudy/report
STUDY = HERE.parent                              # meshStudy
RES = STUDY.parent / "results"
FIG = HERE / "fig"
CASES = ["coarse", "medium", "fine"]
BOX_CELL = {"coarse": 4, "medium": 2, "fine": 1}            # m, inside the cluster box
BOX_MIN, BOX_MAX = (50.0, 145.0, 12.0), (235.0, 262.0, 60.0)  # = setup.py
NITER = 1000
RHO = 1.2                # kg/m3: simpleFoam's p is kinematic (p/rho)
R = 2.0                  # refinement ratio between meshes
GCI_OK = 5.0             # % : common acceptance level for GCI
INCLUDE_RUNNING = "--include-running" in sys.argv
IDX = list(csv.DictReader(open(STUDY / "probes_index.csv")))
# same key probes as make_table.py
KEY = [("upstream", 10), ("gap", 5), ("gap", 20), ("wake", 5), ("wake", 20), ("station", 20), ("station", 40)]


# ======================================================================= parsing (from make_table.py)
FLT = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"   # a number, without a trailing full stop

def num(pattern, text, cast=float):
    # the patterns below are make_table.py's; its number class [\d.eE+-]+ also swallows a full stop
    # ("Min volume = 0.397."), so it is swapped for the stricter FLT here
    pattern = pattern.replace("([\\d.eE+-]+)", f"({FLT})").replace("[\\d.eE+-]+", FLT)
    m = re.findall(pattern, text, re.M)
    return cast(m[-1]) if m else None

def mesh_quality(c):
    f = STUDY / c / "log.checkMesh"
    if not f.exists():
        return None
    t = f.read_text(errors="replace")
    return dict(
        complete=bool(re.search(r"^End\s*$", t, re.M)),
        cells=num(r"^\s+cells:\s+(\d+)", t, int),
        faces=num(r"^\s+faces:\s+(\d+)", t, int),
        skew=num(r"Max skewness = ([\d.eE+-]+)", t),
        skew_bad=num(r"(\d+) highly skew faces", t, int) or 0,
        northo_max=num(r"non-orthogonality Max: ([\d.eE+-]+)", t),
        northo_avg=num(r"non-orthogonality Max: [\d.eE+-]+ average: ([\d.eE+-]+)", t),
        northo_bad=num(r"severely non-orthogonal \(> [\d.]+ degrees\) faces: (\d+)", t, int) or 0,
        aspect=num(r"Max aspect ratio = ([\d.eE+-]+)", t),
        minvol=num(r"Min volume = ([\d.eE+-]+)", t),
        negvol=num(r"Number of negative volume cells: (\d+)", t, int) or 0,
        pyr_bad=num(r"Error in face pyramids: (\d+) faces", t, int) or 0,
        failed=num(r"Failed (\d+) mesh checks", t, int) or 0,
    )

def skew_face_centres(c):
    """Centres of the highly skewed faces, from the VTK set written by checkMesh -writeSets vtk."""
    files = sorted((STUDY / c / "postProcessing").glob("*/skewFaces/*.vtp"))
    if not files:
        return None
    t = files[-1].read_text(errors="replace")
    blobs = dict(re.findall(r"Name='(\w+)'[^>]*format='binary'>\s*([A-Za-z0-9+/=]+)", t))
    try:
        P = np.frombuffer(base64.b64decode(blobs["Points"])[8:], np.float32).reshape(-1, 3)
        off = np.frombuffer(base64.b64decode(blobs["offsets"])[8:], np.int32)
        con = np.frombuffer(base64.b64decode(blobs["connectivity"])[8:], np.int32)
    except Exception:
        return None
    start = np.r_[0, off[:-1]]
    return np.array([P[con[a:b]].mean(axis=0) for a, b in zip(start, off)])

def ran_all_failed(c):
    """run_all.sh writes '!!! <case> FAILED' when a case's Allrun stops with an error."""
    f = STUDY / "run_all.out"
    return f.exists() and f"!!! {c} FAILED" in f.read_text(errors="replace")

def solver_status(c):
    """'finished' / 'running' / 'failed' / 'mesh checked' / 'meshing' / 'not started', and the iteration."""
    d = STUDY / c
    log = d / "log.simpleFoam"
    if log.exists():
        t = log.read_text(errors="replace")
        it = num(r"^Time = (\d+)\s*$", t, int)
        if re.search(r"^End\s*$", t, re.M):
            return "finished", it
        # (not "Floating point exception": the log header says "Floating point exception trapping enabled")
        if re.search(r"FOAM FATAL|FOAM aborting|sigFpe::sigHandler|core dumped", t):
            return "failed", it
        return ("failed" if ran_all_failed(c) else "running"), it
    if ran_all_failed(c):
        return "failed", None
    if (d / "log.checkMesh").exists():
        return "mesh checked", None
    if (d / "log.blockMesh").exists() or (d / "log.snappyHexMesh").exists():
        return "meshing", None
    return "not started", None

def convergence(c):
    f = STUDY / c / "log.simpleFoam"
    if not f.exists():
        return None
    t = f.read_text(errors="replace")
    res = {k: num(rf"Solving for {k}, Initial residual = ([\d.eE+-]+)", t)
           for k in ("Ux", "Uy", "Uz", "p", "k", "epsilon")}
    return dict(iters=num(r"^Time = (\d+)\s*$", t, int), res=res,
                cont=num(r"sum local = ([\d.eE+-]+)", t),
                cont_global=num(r"sum local = [\d.eE+-]+, global = ([\d.eE+-]+)", t),
                cont_cum=num(r"cumulative = ([\d.eE+-]+)", t),
                clock=num(r"ClockTime = ([\d.]+)", t),
                cont_hist=np.array([float(x) for x in re.findall(r"sum local = ([\d.eE+-]+)", t)]))

NUM = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")

def read_probe(c, field):
    """(iterations, values[iter, probe(, 3)], probe indices found by OpenFOAM) or None.
    Rows with the wrong number of values (e.g. half-written last line of a running case) are skipped."""
    f = STUDY / c / "postProcessing" / "probes" / "0" / field
    if not f.exists():
        return None
    lines = f.read_text(errors="replace").splitlines()
    found = [int(m.group(1)) for l in lines if (m := re.match(r"#\s*Probe\s+(\d+)\s*\(", l))]
    ncomp = 3 if field == "U" else 1
    it, vals = [], []
    for l in lines:
        if not l.strip() or l.startswith("#"):
            continue
        v = NUM.findall(l)
        if len(v) != 1 + ncomp * len(found):
            continue
        it.append(float(v[0])); vals.append([float(x) for x in v[1:]])
    if not it:
        return None
    vals = np.array(vals)
    return np.array(it), (vals.reshape(len(it), -1, 3) if ncomp == 3 else vals), found

def read_residual_history(c):
    f = STUDY / c / "postProcessing" / "residuals" / "0" / "solverInfo.dat"
    if not f.exists():
        return None
    head, rows = None, []
    for l in f.read_text(errors="replace").splitlines():
        if l.startswith("# Time"):
            head = l[1:].split()
        elif l.strip() and not l.startswith("#") and head and len(l.split()) == len(head):
            rows.append(l.split())
    if not rows:
        return None
    out = {"it": np.array([float(r[0]) for r in rows])}
    for name in ("Ux", "Uy", "Uz", "p", "k", "epsilon"):
        col = f"{name}_initial"
        if col in head:
            out[name] = np.array([float(r[head.index(col)]) for r in rows])
    return out


# ======================================================================= collect
data = {}
for c in CASES:
    status, it = solver_status(c)
    d = {"status": status, "iter": it, "mesh": mesh_quality(c), "conv": convergence(c),
         "skew_xyz": skew_face_centres(c), "hist": read_residual_history(c)}
    U, P = read_probe(c, "U"), read_probe(c, "p")
    if U is not None:
        d["it"], Uv, d["found"] = U
        d["spd"] = np.linalg.norm(Uv, axis=2)
        d["p"] = P[1] if P is not None else None
    data[c] = d

def usable(c):
    """Probe values are reported for finished runs (and, with --include-running, for running ones)."""
    d = data[c]
    return "found" in d and (d["status"] == "finished" or (INCLUDE_RUNNING and d["status"] in ("running", "failed")))

def probe_index(loc, h):
    return next(int(r["i"]) for r in IDX if r["location"] == loc and float(r["height_above_ground"]) == h)

def probe(c, loc, h, what="spd"):
    if not usable(c):
        return None
    d, i = data[c], probe_index(loc, h)
    if i not in d["found"] or d.get(what) is None:
        return None
    return float(d[what][-1, d["found"].index(i)])

def drift(c, i):
    d = data[c]
    if "found" not in d or i not in d["found"]:
        return None
    s = d["spd"][:, d["found"].index(i)]
    w = s[-100:]
    return 100 * (w.max() - w.min()) / max(abs(s[-1]), 1e-9)

def max_drift(c):
    v = [x for x in (drift(c, int(r["i"])) for r in IDX) if x is not None]
    return max(v) if v else None

def gci(f3, f2, f1):
    """Roache's 3-grid GCI on the fine grid (f1 = fine), as in make_table.py.
    Returns (GCI %, observed order p, short text for the table)."""
    if None in (f1, f2, f3):
        return None, None, "--"
    e32, e21 = f3 - f2, f2 - f1
    if abs(e21) < 1e-12 * max(1.0, abs(f1)):
        return 0.0, None, "0.0"
    if abs(e32) < 1e-12 * max(1.0, abs(f1)) or e32 * e21 < 0:
        return None, None, "osc."
    p = math.log(abs(e32 / e21)) / math.log(R)
    p_used = min(max(p, 0.5), 2.0)          # common safeguard: clip p to [0.5, 2]
    g = 100 * 1.25 * abs(e21 / f1) / (R ** p_used - 1)
    return g, p, f"{g:.1f}"

def change(a, b):
    return None if None in (a, b) or abs(a) < 1e-12 else 100 * (b - a) / a

def dp(c):
    a, b = probe(c, "upstream", 5, "p"), probe(c, "wake", 5, "p")
    return None if None in (a, b) else RHO * (a - b)


# ======================================================================= LaTeX formatting
DASH = "--"
def f_(v, fmt="{:.2f}"):
    return DASH if v is None else fmt.format(v)

def sci(v):
    if v is None:
        return DASH
    m, e = f"{v:.1e}".split("e")
    return rf"${m}\times10^{{{int(e)}}}$"

def drift_txt(v):
    return DASH if v is None else ("$<0.01$" if v < 0.005 else f"{v:.2f}")

def pct(v):
    return DASH if v is None else rf"${v:+.1f}$"

def big(n):
    return DASH if n is None else f"{n:,}".replace(",", r"\,")

def bold_if(s, cond):
    return rf"\textbf{{{s}}}" if cond and s != DASH else s

def row(*cells):
    return " & ".join(cells) + r" \\"

DAGGER = r"$^\dagger$"
def head(c):
    mark = DAGGER if INCLUDE_RUNNING and data[c]["status"] in ("running", "failed") and usable(c) else ""
    return rf"{c}{mark} ({BOX_CELL[c]}\,m)"

def status_word(c):
    s, it = data[c]["status"], data[c]["iter"]
    if s == "running":
        return f"running (iteration {it} of {NITER})" if it else "running"
    if s == "failed":
        return f"failed at iteration {it}" if it else "failed"
    return {"mesh checked": "meshed, solver not started", "meshing": "meshing"}.get(s, s)


# ---------------------------------------------------------------- Table 1: mesh quality
def table_mesh():
    M = {c: data[c]["mesh"] for c in CASES}
    g = lambda c, k: None if M[c] is None else M[c][k]
    rows = [
        ("Cell size in the box (m)", DASH, lambda c: str(BOX_CELL[c]) if M[c] else DASH),
        ("Total cells", DASH, lambda c: big(g(c, "cells"))),
        ("Max skewness", r"$<4$", lambda c: bold_if(f_(g(c, "skew")), (g(c, "skew") or 0) > 4)),
        ("Faces with skewness $>4$", "0", lambda c: bold_if(f_(g(c, "skew_bad"), "{:d}"), (g(c, "skew_bad") or 0) > 0)),
        (r"Max non-orthogonality (\textdegree)", r"$<70$",
         lambda c: bold_if(f_(g(c, "northo_max"), "{:.1f}"), (g(c, "northo_max") or 0) > 70)),
        (r"Average non-orthogonality (\textdegree)", DASH, lambda c: f_(g(c, "northo_avg"), "{:.1f}")),
        (r"Faces with non-orthogonality $>70$\textdegree", "0",
         lambda c: bold_if(f_(g(c, "northo_bad"), "{:d}"), (g(c, "northo_bad") or 0) > 0)),
        ("Max aspect ratio", r"$<1000$", lambda c: bold_if(f_(g(c, "aspect"), "{:.1f}"), (g(c, "aspect") or 0) > 1000)),
        (r"Min cell volume (m$^3$)", r"$>0$", lambda c: bold_if(f_(g(c, "minvol"), "{:.3g}"), g(c, "minvol") is not None and g(c, "minvol") <= 0)),
        ("Invalid cells (zero or negative volume)", "0", lambda c: bold_if(f_(g(c, "negvol"), "{:d}"), (g(c, "negvol") or 0) > 0)),
        ("Wrongly oriented faces (face pyramids)", "0", lambda c: bold_if(f_(g(c, "pyr_bad"), "{:d}"), (g(c, "pyr_bad") or 0) > 0)),
        ("Failed checkMesh checks", "0", lambda c: bold_if(f_(g(c, "failed"), "{:d}"), (g(c, "failed") or 0) > 0)),
    ]
    L = [r"\begin{table}[htbp]", r"\centering", r"\small",
         r"\caption{Mesh quality from \texttt{checkMesh}. Values outside the usual limit are in bold. "
         r"A dash means the mesh does not exist yet.}", r"\label{tab:quality}",
         r"\begin{tabular}{llrrr}", r"\toprule",
         row("Quantity", "Usual limit", *[head(c) for c in CASES]), r"\midrule"]
    L += [row(name, lim, *[fn(c) for c in CASES]) for name, lim, fn in rows]
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return L


# ---------------------------------------------------------------- Table 2: convergence
def table_conv():
    C = {c: data[c]["conv"] for c in CASES}
    r_ = lambda c, k: None if C[c] is None else C[c]["res"][k]
    rows = [
        ("Status", lambda c: status_word(c)),
        ("Iterations done", lambda c: f_(None if C[c] is None else C[c]["iters"], "{:d}")),
        (r"Initial residual $U_x$", lambda c: sci(r_(c, "Ux"))),
        (r"Initial residual $U_y$", lambda c: sci(r_(c, "Uy"))),
        (r"Initial residual $U_z$", lambda c: sci(r_(c, "Uz"))),
        (r"Initial residual $p$", lambda c: sci(r_(c, "p"))),
        (r"Initial residual $k$", lambda c: sci(r_(c, "k"))),
        (r"Initial residual $\varepsilon$", lambda c: sci(r_(c, "epsilon"))),
        ("Continuity error, sum local", lambda c: sci(None if C[c] is None else C[c]["cont"])),
        ("Continuity error, cumulative", lambda c: sci(None if C[c] is None else C[c]["cont_cum"])),
        (r"Max probe drift, last 100 it.\ (\%)", lambda c: drift_txt(max_drift(c))),
        ("Wall-clock time (h)", lambda c: f_(None if C[c] is None or C[c]["clock"] is None else C[c]["clock"] / 3600, "{:.1f}")),
    ]
    L = [r"\begin{table}[htbp]", r"\centering", r"\small",
         r"\caption{Iterative convergence at the last iteration of each run (values of a running case are "
         r"the latest ones). Residuals are the initial residuals of the last iteration. "
         r"Probe drift is the largest change of any probe's $|U|$ over the last 100 iterations, relative to its final value.}",
         r"\label{tab:convergence}",
         r"\begin{tabular}{lrrr}", r"\toprule",
         row("Quantity", *[head(c) for c in CASES]), r"\midrule"]
    L += [row(name, *[fn(c) for c in CASES]) for name, fn in rows]
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return L


# ---------------------------------------------------------------- Table 3 / A1: mesh independence
def indep_rows(pairs):
    out = []
    for loc, h in pairs:
        f3, f2, f1 = (probe(c, loc, h) for c in CASES)
        g, p, gtxt = gci(f3, f2, f1)
        name = f"{loc}, {h:g}\\,m" + (" (control)" if loc == "upstream" else "")
        out.append(row(name, f_(f3), f_(f2), f_(f1), pct(change(f3, f2)), pct(change(f2, f1)),
                       f_(p, "{:.2f}"), bold_if(gtxt, g is not None and g > GCI_OK)))
    return out

INDEP_HEAD = [r"\begin{tabular}{lrrrrrcr}", r"\toprule",
              r" & \multicolumn{3}{c}{$|U|$ (m/s)} & \multicolumn{2}{c}{Change (\%)} & & \\",
              r"\cmidrule(lr){2-4}\cmidrule(lr){5-6}",
              None,  # filled in below
              r"\midrule"]

def indep_head():
    h = list(INDEP_HEAD)
    h[4] = row("Probe (height above ground)", *[head(c) for c in CASES],
               r"c.$\to$m.", r"m.$\to$f.", r"$p$", r"GCI$_\mathrm{fine}$ (\%)")
    return h

def table_indep():
    dps = [dp(c) for c in CASES]
    g, p, gtxt = gci(*dps)
    L = [r"\begin{table}[htbp]", r"\centering", r"\small",
         r"\caption{Mesh independence: wind speed $|U|$ at the key probes at the last iteration, its change from one "
         r"mesh to the next (relative to the coarser mesh), the observed order $p$ (as computed; the GCI uses $p$ "
         r"limited to 0.5--2) and the grid convergence index on the fine mesh. The last row is the pressure drop $\Delta p = \rho\,(p_\mathrm{upstream}-p_\mathrm{wake})$ "
         r"at 5\,m above ground. c.$\to$m.: coarse to medium; m.$\to$f.: medium to fine. "
         r"``osc.'': the three values do not change monotonically, so no GCI. "
         rf"GCI above {GCI_OK:g}\,\% in bold. A dash: case not finished yet."
         + (r" $^\dagger$: run not finished, latest values." if INCLUDE_RUNNING else "") + "}",
         r"\label{tab:independence}"] + indep_head()
    L += indep_rows(KEY)
    L += [r"\midrule",
          row(r"$\Delta p$ upstream$\to$wake, 5\,m (Pa)", *[f_(v, "{:.3f}") for v in dps],
              pct(change(dps[0], dps[1])), pct(change(dps[1], dps[2])), f_(p, "{:.2f}"),
              bold_if(gtxt, g is not None and g > GCI_OK))]
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return L

def table_all_probes():
    L = [r"\begin{table}[htbp]", r"\centering", r"\small",
         r"\caption{All 19 probes (same columns as Table~\ref{tab:independence}). Heights are above the local ground.}",
         r"\label{tab:allprobes}"] + indep_head()
    locs = list(dict.fromkeys(r["location"] for r in IDX))
    for k, loc in enumerate(locs):
        if k:
            L.append(r"\midrule")
        L += indep_rows([(loc, float(r["height_above_ground"])) for r in IDX if r["location"] == loc])
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return L


# ======================================================================= sentences (summary.tex)
def esc(s):
    return s.replace("%", r"\%")

def quality_sentence(c):
    m = data[c]["mesh"]
    if m is None:
        return None
    if not m["complete"]:
        return f"{c}: checkMesh log incomplete."
    if m["failed"] == 0:
        return (f"{c} mesh ({big(m['cells'])} cells): passes all checkMesh checks "
                f"(max skewness {f_(m['skew'])}, max non-orthogonality {f_(m['northo_max'], '{:.1f}')}\\textdegree).")
    probs = []
    if (m["skew"] or 0) > 4:
        probs.append(f"max skewness {f_(m['skew'])} on {m['skew_bad']} face{'s' if m['skew_bad'] != 1 else ''} (limit 4)")
    if m["northo_bad"]:
        probs.append(f"{m['northo_bad']} faces with non-orthogonality above 70\\textdegree")
    if m["negvol"]:
        probs.append(f"{m['negvol']} cells with zero or negative volume")
    if m["pyr_bad"]:
        probs.append(f"{m['pyr_bad']} wrongly oriented faces")
    s = (f"{c} mesh ({big(m['cells'])} cells): fails {m['failed']} checkMesh check{'s' if m['failed'] != 1 else ''}"
         + (": " + "; ".join(probs) if probs else "") + ".")
    xyz = data[c]["skew_xyz"]
    if xyz is not None and len(xyz):
        inside = int(np.sum(np.all((xyz >= BOX_MIN) & (xyz <= BOX_MAX), axis=1)))
        P = np.array([[float(r["x"]), float(r["y"]), float(r["z"])] for r in IDX])
        dmin = float(np.min(np.linalg.norm(xyz[:, None, :] - P[None, :, :], axis=2)))
        where = ("all outside the refinement box" if inside == 0 else
                 f"{inside} of them inside the refinement box")
        s += f" The skewed faces are {where}; the nearest is {dmin:.0f}\\,m from a probe."
    return s

def conv_sentence(c):
    d, v = data[c], data[c]["conv"]
    if v is None:
        return None
    r = {k: x for k, x in v["res"].items() if x is not None}
    if not r:
        return None
    lo, hi = min(r, key=r.get), max(r, key=r.get)
    nice = lambda k: {"epsilon": r"$\varepsilon$", "Ux": "$U_x$", "Uy": "$U_y$", "Uz": "$U_z$"}.get(k, f"${k}$")
    when = f"after {v['iters']} iterations" if d["status"] == "finished" else f"at iteration {v['iters']} (still {d['status']})"
    s = (f"{c}: initial residuals {when} between {sci(r[lo])} ({nice(lo)}) and {sci(r[hi])} ({nice(hi)}); "
         f"continuity error (sum local) {sci(v['cont'])}")
    md = max_drift(c)
    if md is not None:
        s += (f"; the probes moved by {'less than 0.01' if md < 0.005 else f'at most {md:.2f}'}\\,\\% "
              "over the last 100 iterations")
    return s + "."

def indep_sentences():
    out = []
    inbox = [(l, h) for l, h in KEY if l != "upstream"]
    for a, b in (("coarse", "medium"), ("medium", "fine")):
        ch = [(abs(change(probe(a, l, h), probe(b, l, h))), l, h) for l, h in inbox
              if change(probe(a, l, h), probe(b, l, h)) is not None]
        if not ch:
            continue
        m, l, h = max(ch)
        ctrl = change(probe(a, "upstream", 10), probe(b, "upstream", 10))
        s = (f"From {a} to {b}, $|U|$ at the key probes in the campus changes by at most {m:.1f}\\,\\% "
             f"(largest at {l}, {h:g}\\,m) and by {f_(None if ctrl is None else abs(ctrl), '{:.1f}')}\\,\\% at the upstream control")
        dch = change(dp(a), dp(b))
        if dch is not None:
            s += f"; the pressure drop changes by {dch:+.1f}\\,\\%"
        out.append(s + ".")
    gs = [gci(*(probe(c, l, h) for c in CASES)) for l, h in KEY]
    valid = [g for g, _, t in gs if g is not None]
    osc = sum(1 for _, _, t in gs if t == "osc.")
    if valid or osc:
        out.append(f"GCI$_\\mathrm{{fine}}$ is below {GCI_OK:g}\\,\\% at {sum(g < GCI_OK for g in valid)} of the "
                   f"{len(KEY)} key probes" + (f"; at {osc} probe{'s' if osc != 1 else ''} the three values are not "
                   f"monotonic, so no GCI is defined there" if osc else "") + ".")
    return out

def summary_tex():
    now = datetime.now().strftime("%H:%M on %d %b %Y")
    st = "; ".join(f"{c}: {status_word(c)}" for c in CASES)
    L = ["% generated by make_tex.py -- do not edit, rerun build.sh instead",
         rf"\newcommand{{\StudyStatus}}{{Status of the runs when this document was built ({now}): {st}.}}"]
    for c in CASES:
        m = data[c]["mesh"]
        L.append(rf"\newcommand{{\Cells{c.capitalize()}}}{{{big(m['cells']) if m else DASH}}}")
    qs = [s for s in (quality_sentence(c) for c in CASES) if s]
    cs = [s for s in (conv_sentence(c) for c in CASES) if s]
    ins = indep_sentences()
    def items(lst, empty):
        return (r"\begin{itemize}" + "".join(rf"\item {x}" for x in lst) + r"\end{itemize}") if lst else rf"\textit{{{empty}}}"
    L.append(rf"\newcommand{{\FindingsQuality}}{{{items(qs, 'No mesh has been checked yet.')}}}")
    L.append(rf"\newcommand{{\FindingsConvergence}}{{{items(cs, 'No solver run has started yet.')}}}")
    L.append(rf"\newcommand{{\FindingsIndependence}}{{{items(ins, 'At least two finished runs are needed to compare meshes; not available yet.')}}}")
    m = data["coarse"]["mesh"]
    if m is None:
        sk = "The quality of the coarse mesh is not available yet."
    elif (m["skew"] or 0) > 4:
        sk = (f"The coarse mesh fails the skewness check: the maximum skewness is {f_(m['skew'])} on "
              f"{m['skew_bad']} of its {big(m['faces'])} faces, against a usual limit of 4.")
    else:
        sk = f"The coarse mesh passes the skewness check (maximum {f_(m['skew'])}, usual limit 4)."
    L.append(rf"\newcommand{{\SkewSentence}}{{{sk}}}")
    n_fin = sum(data[c]["status"] == "finished" for c in CASES)
    L.append(rf"\newcommand{{\NFinished}}{{{n_fin}}}")
    return "\n".join(L) + "\n"


# ======================================================================= residual plot
SERIES = [("Ux", "$U_x$", "#2a78d6"), ("Uy", "$U_y$", "#eb6834"), ("Uz", "$U_z$", "#1baf7a"),
          ("p", "$p$", "#eda100"), ("k", "$k$", "#e87ba4"), ("epsilon", r"$\varepsilon$", "#008300")]

def residual_plot(path):
    fig, axs = plt.subplots(1, 3, figsize=(10, 3.4), dpi=160, sharey=True, layout="constrained")
    for ax, c in zip(axs, CASES):
        h, v = data[c]["hist"], data[c]["conv"]
        ax.set_title(f"{c} ({BOX_CELL[c]} m in the box)", fontsize=10)
        if h is None:
            ax.text(0.5, 0.5, "no residuals yet", transform=ax.transAxes, ha="center", va="center", color="0.5")
        else:
            for k, lab, col in SERIES:
                if k in h:
                    ax.semilogy(h["it"], h[k], color=col, lw=1.2, label=lab)
        if v is not None and len(v["cont_hist"]):
            ax.semilogy(np.arange(1, len(v["cont_hist"]) + 1), v["cont_hist"], color="k", lw=1.0, ls="--",
                        label="continuity\n(sum local)")
        ax.set_xlim(0, NITER); ax.set_xlabel("iteration"); ax.grid(alpha=0.3, which="major")
    axs[0].set_ylabel("initial residual,\ncontinuity error")
    handles, labels = [], []
    for ax in axs:
        for hh, ll in zip(*ax.get_legend_handles_labels()):
            if ll not in labels:
                handles.append(hh); labels.append(ll)
    if handles:
        fig.legend(handles, labels, loc="outside right center", fontsize=9, frameon=False)
    fig.savefig(path); plt.close(fig)


# ======================================================================= main
def bare_tabular(lines):
    """Only the tabular of a table float (for the slides)."""
    a = next(i for i, l in enumerate(lines) if l.startswith(r"\begin{tabular}"))
    b = next(i for i, l in enumerate(lines) if l.startswith(r"\end{tabular}"))
    return lines[a:b + 1]

def main():
    t_mesh, t_conv, t_indep = table_mesh(), table_conv(), table_indep()
    tables = (["% generated by make_tex.py -- do not edit, rerun build.sh instead"]
              + t_mesh + [""] + t_conv + [""] + t_indep + [""])
    (HERE / "tables.tex").write_text("\n".join(tables) + "\n")
    for name, lines in (("quality", t_mesh), ("convergence", t_conv), ("independence", t_indep)):
        (HERE / f"tab_{name}.tex").write_text("\n".join(["% generated by make_tex.py"] + bare_tabular(lines)) + "\n")
    (HERE / "tables_appendix.tex").write_text("\n".join(["% generated by make_tex.py"] + table_all_probes()) + "\n")
    (HERE / "summary.tex").write_text(summary_tex())
    # pictures: fresh copies, so the folder is self-contained (Overleaf: zip report/ and upload)
    FIG.mkdir(exist_ok=True)
    for old in FIG.glob("*.png"):
        old.unlink()
    n = 0
    for src in list((RES / "mesh_fig").glob("*.png")) + [RES / "mesh_study_convergence.png", RES / "mesh_study_profiles.png"]:
        if src.exists():
            shutil.copy2(src, FIG / src.name); n += 1
    residual_plot(FIG / "residuals.png")
    print("status: " + "; ".join(f"{c}: {status_word(c)}" for c in CASES))
    print(f"wrote {HERE / 'tables.tex'}, tables_appendix.tex, tab_*.tex (slides), summary.tex; "
          f"{n} pictures copied to fig/ + residuals.png")


if __name__ == "__main__":
    main()
