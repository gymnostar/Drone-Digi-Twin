"""Planner table + one-paragraph reading for slides.tex, from DLR/planner/results/plan_summary.json.
Uses the 'fair' comparison (straight-line shortest vs energy-aware-smoothed energy path) when present,
otherwise falls back to the raw graph paths."""


def _mode(v, base):
    for prefix in ("fair_", "raw_"):
        if prefix + base in v:
            return prefix + base, prefix.rstrip("_")
    return None, None


def build(plan):
    per = plan.get("per_direction", {})
    if not per:
        return "\\textit{planner results not available yet}\n", "--"
    rows, label = [], None
    for key in sorted(per, key=lambda k: (k.split("_")[1] != "edmo", int(k.split("_")[0]))):
        v = per[key]
        th, kind = key.split("_")
        mb, label = _mode(v, "best")
        mf, _ = _mode(v, "fix10")
        if not mb:
            continue
        b, f = v[mb], v[mf]
        rows.append(f"{th}$^\\circ$ & {v['wind_speed']:.1f} & {f['save_mean']:.1f} & {b['save_mean']:.1f} & "
                    f"{b['save_median']:.1f} & {b['save_p90']:.1f} & {b['extra_len_mean']:.1f} \\\\")
    w = plan.get("weighted", {})
    t = ["\\footnotesize\\begin{tabular}{@{}lrrrrrr@{}}", "\\toprule",
         " & & \\multicolumn{4}{c}{energy saved vs.\\ shortest path (\\%)} & \\\\",
         "\\cmidrule(lr){3-6}",
         "Wind from & $U_{10}$ (m/s) & fixed 10\\,m/s & best speed & median & p90 & extra length (\\%) \\\\",
         "\\midrule"]
    t += rows
    t.append("\\midrule")
    for kind, name in (("edmo", "climate-weighted, EDMO speeds"), ("strong", "climate-weighted, 8\\,m/s")):
        e = w.get(kind, {}).get("expected_save_pct", {})
        kb = f"{label}_best" if f"{label}_best" in e else None
        kf = f"{label}_fix10" if f"{label}_fix10" in e else None
        if kb:
            t.append(f"\\multicolumn{{2}}{{@{{}}l}}{{{name}}} & {e[kf]:.1f} & {e[kb]:.1f} & & & \\\\")
    t += ["\\bottomrule", "\\end{tabular}"]
    ndir = len(plan.get("directions_available", []))
    cmp = ("straight line vs.\\ energy path with energy-aware smoothing" if label == "fair"
           else "raw graph paths (both on the 26-neighbour grid)")
    s = (f"{plan.get('n_pairs', '--')} random start--goal pairs over the campus ($\\geq$\\,250\\,m apart), both legs, "
         f"{ndir} wind directions. Comparison: {cmp}.\n\n\\medskip\n"
         "\\textbf{Reading:} at typical winds the wind-aware \\emph{route} saves a few per cent; at 8\\,m/s it saves more, "
         "with a long tail of trips that save 15--25\\,\\%. Choosing the \\emph{speed} for the local wind adds more "
         "than choosing the route. The direction of travel matters most: the same straight path costs far more upwind than downwind.")
    return "\n".join(t) + "\n", s
