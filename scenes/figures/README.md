# Figures for the cheaper-mesh wind study

These are the presentation figures for the Niigata buildings plus SWiFT turbine study: how much cheaper a wind map gets when interior mesh cells are merged, how close it stays to the strict reference, and what the released lookup returns. Every figure is a 1920 px wide PNG on a white background, built only from files in `scenes/` (`joint/results/`, `joint/lookup/`, `review-export/`, `buildings/results/`, `turbine/results/`). No solver, mesher or Docker job runs when they are rebuilt.

One palette is used throughout: blue, orange, aqua and yellow mark categories, a light-to-dark blue ramp shows wind speed or cell size, and grey is buildings or "unchanged". Labels and marker shapes repeat what the colour says, so no figure relies on colour alone.

**Main caveat.** The merged (production) mesh in fig03 and fig04 is **reconstructed, not loaded**. Only the 422,894-cell strict mesh is stored (`review-export/joint-baseline.vtu`). `src/recon.py` re-applies the merge rules of `joint/agglomeration/agglomerate.C` to it and reproduces the recorded result exactly (27,428 groups, 230,898 cells). No production flow volume is stored, so no figure shows the merged-mesh velocity field in 3D.

## Rebuild

From the repository root, with [uv](https://docs.astral.sh/uv/) (any Python 3.12 virtual environment with these four packages works):

```sh
uv venv -p 3.12 .venv-figures
uv pip install --python .venv-figures numpy matplotlib pyvista vtk
cd scenes/figures/src
for s in render3d fig_turbine fig_data fig_nne fig_fail fig_ne fig_flow fig_query; do ../../../.venv-figures/bin/python $s.py; done
../../../.venv-figures/bin/python recon.py /tmp/recon.npz && ../../../.venv-figures/bin/python fig_mesh.py /tmp
```

The figures are written into this folder; set `FIG_OUT=/some/folder` to write them elsewhere. `recon.py` takes about a minute and prints the reconstructed merge counts next to the recorded ones. With numpy 2.5, matplotlib 3.11, pyvista 0.49 and vtk 9.7 on macOS, the scripts reproduce the stored PNGs pixel for pixel. Another font (the scripts prefer Helvetica, then DejaVu Sans) shifts text slightly. `01-turbine-result.png`, `03-city-limit.png`, `niigata-map.png` and `swift-sketch.png` come from earlier slides and have no script here.

| Script | Figures |
|---|---|
| `render3d.py` | fig01 |
| `fig_turbine.py` | fig02 |
| `recon.py`, then `fig_mesh.py` | fig03, fig04 |
| `fig_data.py` | fig05 to fig08 |
| `fig_nne.py` | fig09 |
| `fig_fail.py` | fig10 |
| `fig_ne.py` | figNE |
| `fig_flow.py` | figflow |
| `fig_query.py` | figquery |
| `style.py` | shared palette, fonts and paths |

---

## figflow-pipeline.png: from geometry to a wind value

- **Purpose:** show the whole pipeline once, so every later figure has a place in it.
- **How to read it:** eight steps, top row left to right, then the bottom row. Step 3 (green) is the only change in the cheap "production" recipe; the dashed arrow shows the strict reference skipping it. Step 6 (orange) is the accuracy gate. The text on the right is the cold time of steps 2–5 for each released heading.
- **What to point at:** the strict and production runs share geometry, meshing, solver and sampling; they differ in the merge and in looser stopping rules (639 vs 308 iterations for N, 426 vs 136 for NE). Only a field that passes step 6 becomes a lookup table.
- **Takeaway:** a wind map is one CFD solve per heading, checked against a strict solve, then stored as 74,277 points that answer a query in well under a microsecond.
- **Caveats:** every number is read from `joint/results/*/metrics.json`, `joint/lookup/index.json` and `joint/lookup/query-timing.json` when the figure is built. Times are one cold run each on Studio (3 CPUs, 5 GiB). The comparison (step 6) is not part of the timed pipeline. Exporting a lookup table took 9.6 s (N) and 9.2 s (NE) and is not timed either. Query times are warm in-process Python calls on one Mac, not a drone-in-the-loop latency.

## figNE-heading-difference.png: why each heading needs its own solve

- **Purpose:** show how much the field changes when the wind turns by 45°, and where.
- **How to read it:** four plan views at hub height (32.1 m) on the 5 m lookup grid, east to the right and north up. Panels 1 and 2 are wind speed for wind from N and from NE, on one shared 0–10 m/s scale. Panel 3 is the length of the vector difference, |u(NE) − u(N)|, in m/s. Panel 4 is the speed difference, NE minus N (blue: faster from NE, red: faster from N). Grey cells have no value because they are inside a building; at this height only the 60 m tower near (−60, 0) reaches up (24 masked cells in each table). The black bar is the rotor, turned to face each wind; panels 3 and 4 show both orientations.
- **What to point at:** the tower wake points south for N and south-west for NE, and the rotor wake swings the same way. In panel 3 almost every point differs by about 6.5 m/s (median 6.54, P95 6.89, maximum 10.90 m/s). That is mostly the turn itself: two 8.7 m/s winds 45° apart differ by 6.66 m/s (the tick on the colour bar). Panel 4 removes direction and shows where the wakes and building shadows moved: speed changes by more than 1 m/s at 21% of hub-height points, by up to 8.4 m/s.
- **Takeaway:** a field for one heading cannot stand in for another. The wakes move with the wind, and the nonlinear flow around buildings changes shape, not just direction. That is why every heading needs its own solve and why the in-between heading fallback is marked unvalidated.
- **Caveats:** one height only; near the ground (5 m) the masks differ slightly between headings (152 cells), at 32.1 m they are identical. These are the released production fields, not strict references. The NE table comes from the `table-NE-simplec` run (136 iterations); the 1.99× timing is its same-recipe retime (`table-NE-simplec-retimed`), which has the same iteration count and the same zone errors.

## figquery-lookup-example.png: what the lookup returns

- **Purpose:** make the product tangible: a position goes in, a wind vector comes out.
- **How to read it:** every value is a real call to `scenes/joint/lookup.py` on the released tables. Left: wind speed along an east–west line 70 m south of the hub at 32.1 m. Markers are grid nodes (`at()`), the line is `interpolate()` every 0.5 m. Right: vertical profile 70 m downwind of the hub for each heading: (0, −350) for N and (−50, −330) for NE, with the 7 stored heights as markers and the rotor height range shaded. Bottom: four example calls and exactly what they returned (east, north, up velocity in m/s).
- **What to point at:** for wind from N the line crosses two dips: the rotor wake at east 0 (4.94 m/s) and the 60 m tower wake near east −55 (5.17 m/s). For NE the rotor wake has moved to east −70 (6.10 m/s at this line, which is about 99 m downwind for NE). A position inside the tower returns `None`, not a made-up value. A heading of 22.5° returns the N field but labels it `nearest_released_unvalidated`.
- **Takeaway:** the lookup answers in about 0.14 µs per grid-node call and 0.71 µs interpolated, and it says clearly when it has no value or when the heading is not validated.
- **Caveats:** query times are the medians of 5 warm batches of 10,000 Python calls in `joint/lookup/query-timing.json` (one Mac). Interpolation is trilinear between nodes; its separate error against CFD is in `joint/lookup/N-offgrid-check.json` and `NE-offgrid-check.json`. Values are a numerical solution, not measured wind.

## fig01-coupled-scene-3d.png: the coupled scene in 3D

- **Purpose:** show what "buildings and turbine in one air volume" looks like before any numbers.
- **How to read it:** grey blocks are the Niigata buildings at full scale. The large translucent sheet is a horizontal slice through the air at 32.1 m (the rotor hub height). The upright sheet is a vertical slice along x = 0, running through the rotor and its wake. The orange ring is the actuator disk: the region where the solver brakes the air to represent the turbine. Colour is wind speed from 0 m/s (near white) to 10 m/s (dark blue). Wind comes from the north (arrow).
- **What to point at:** the pale streak behind the ring, which is the turbine wake (slower air), and the pale patch behind the 60 m tower, which is the building wake. Most other buildings sit below the 32.1 m slice, so the slice there shows nearly undisturbed wind.
- **Takeaway:** one solve captures both wakes and their interaction. Separate maps are never added together.
- **Caveats:** this is the saved `joint-baseline` field: the 422,894-cell common mesh solved with the cheaper stopping rules (364 iterations). It is not the strict reference and not the merged production field, because neither of those volumes is stored. The view is cropped to x ±230 m, north −470 to 230 m and heights up to 130 m. The ring is a marker, not blade geometry.

## fig02-turbine-alone-mesh.png: turbine alone, reference vs cheap mesh

- **Purpose:** show the change that worked for the isolated turbine: bigger cells in the wake, with the force cells kept fine.
- **How to read it:** each panel is a vertical slice through the hub. The x axis is distance downstream of the rotor in metres and the y axis is height in metres. Every outlined polygon is one mesh cell, and its colour is the wind speed in that cell (4.5 to 10 m/s). The orange bar is the 27 m rotor.
- **What to point at:** the wake region between 0 and about 180 m. Top panel: about 1.7 m cells. Bottom panel: about 6.75 m cells, four times wider. In both panels the cells at the rotor stay equally fine, so the same 6,496 force cells receive the force. The slow-air plume has the same shape and length in both panels.
- **Takeaway:** 359,334 cells became 112,206 and map generation dropped from 131.1 s to 22.1 s (5.94× in development, 5.97–6.09× on three held-out conditions), while wake error stayed at 1.38–1.63%.
- **Caveats:** the panels show the development condition, not the three held-out conditions. The timings are one run each. The thin curved lines are slice artefacts of non-convex polyhedral cells, not flow features.

## fig03-mesh-before-after.png: where the coupled mesh got coarser (side view)

- **Purpose:** explain "merge groups of 8 interior cells" on the real coupled mesh, and show what stays fine.
- **How to read it:** both panels are the same vertical north–south section through the rotor (x = 0). The x axis is the north coordinate in metres, and the wind blows from right to left. The y axis is height. Each polygon is one cell, coloured by edge length from 1.25 m (darkest) to 40 m (lightest). The orange dashed box is the protected zone around the rotor, and the orange bar is the disk. White gaps at the bottom are buildings.
- **What to point at:** (1) Above 40 m over the city, 10 m and 20 m cells became 20 m and 40 m cells. (2) In the wake region on the left (north −460 to −250 m), the 2.5 m cells became 5 m cells, except inside the protected box. (3) Cells touching buildings and the ground stay unchanged. (4) The rotor cells are identical in both panels.
- **Takeaway:** the merge removes 191,996 cells (422,894 to 230,898) without moving any surface or rotor-force cell.
- **Caveats:** the merged mesh is reconstructed, not loaded (see the main caveat above). It reproduces the recorded result exactly: 27,428 groups, split by original refinement level as [7772, 2395, 1304, 15957], giving 230,898 cells. Cell sizes for non-box (snapped) cells are estimated from their bounding boxes. Thin slanted lines are slice artefacts of concave cells.

## fig04-mesh-what-changed.png: which cells were merged (near ground) and how many by size

- **Purpose:** show that merging happens in open air and stops at walls, and give the whole-mesh count.
- **How to read it:** left panel: a plan view at 6 m height. Aqua cells were merged (8 became 1), grey cells are unchanged, and white areas are buildings. The orange dashed box and bar mark the rotor and its protected zone. Right panel: the number of cells of each size in the whole domain, before (blue) and after (aqua), with exact counts printed.
- **What to point at:** the grey band one or two cells wide that hugs every building, and the aqua streets, open ground and wake strip. On the right, the 2.5 m class drops from 229,552 to 101,896 cells, and 7,772 new 40 m cells appear aloft. The 1.25 m class, which holds the rotor-force cells, is unchanged at 17,043.
- **Takeaway:** cells next to surfaces, and the force cells, keep their exact geometry. Only interior cells merge, so the counts fall where the flow is smooth.
- **Caveats:** the same reconstruction caveat as fig03 applies. Size classes are approximate for snapped cells. A 6 m slice shows only one height. At 32.1 m fewer cells change: about 22% of cells whose centres lie within 0.7 m of that height. The reason is that most merging happens near the ground band and above 40 m.

## fig05-timing-breakdown.png: where the time goes

- **Purpose:** show which stage the savings come from, and why NE stops just short of 2×.
- **How to read it:** each bar is one complete cold map generation, from scratch, timed once on Studio with 3 CPUs and 5 GiB. Segments are mesh (blue), solve (orange), extracting the comparison samples (aqua), and setup plus automated host checks (yellow). Numbers inside the bars are seconds. The dashed line in each group is half the strict reference time, which a cheap run must beat for 2×.
- **What to point at:** meshing costs 118–126 s in **every** run, because production rebuilds the reference's snapped mesh and then merges it. All savings come from the solve: 393 s to 106 s for N, and 267 s to 59 s for NE with SIMPLEC. NE SIMPLEC finishes at 200.4 s against a 199.4 s target, about 1 s late.
- **Takeaway:** with meshing fixed at about 122 s, the solve is the only remaining lever. NE misses 2× by 0.5%.
- **Caveats:** one timing per recipe. Earlier runs show mesh times varying between 121.9 and 123.0 s, so a 1 s margin is inside the noise. The N production mesh segment includes a failed first audit and its recovery (113.2 + 12.7 s). The "same mesh, looser stopping" row uses the 10× looser residual limits of the production recipe, which alone gives 1.65×. The NE SIMPLE row was originally compared with a 401.14 s reference from an earlier session (1.69×). The dashed line here uses the same-session 398.86 s retime.

## fig06-verdict-speedup-vs-error.png: one-slide verdict

- **Purpose:** put every recipe on one chart against both targets.
- **How to read it:** the x axis is speed-up (strict reference time divided by cheap time). The y axis is vector RMSE in the worst of the four zones, as a percentage of the 8.7 m/s inflow. The shaded box is the target: at least 2× faster and at most 5% error. Blue dots: turbine alone (three held-out conditions). Hollow orange ring: city alone. Aqua markers: coupled recipes, each marker shape labelled.
- **What to point at:** the turbine sits deep inside the target. Coupled N is just inside (2.11×, 2.22%). Both NE versions sit at the 2× line with about 4.5% error, the largest coupled error. The square shows that loosening the stopping rules alone gives 1.65× with almost no error.
- **Takeaway:** the turbine clearly passes, coupled N passes, coupled NE narrowly misses on cost, and the city alone was never certified.
- **Caveats:** the city point compares against an **unconverged** reference, so its 8.6× is descriptive only and no speed-up is claimed. The turbine RMSE covers the wake zone only. NE SIMPLE uses the 401.14 s reference and NE SIMPLEC the 398.86 s same-session reference. All errors are against a numerical reference, not real wind.

## fig07-errors-by-zone.png: accuracy per zone

- **Purpose:** show where the cheap fields differ most from the strict reference.
- **How to read it:** left panel: vector RMSE. Right panel: P95, the error that 95% of points stay below. Both are percentages of the 8.7 m/s inflow, shown for the four gated zones (city, building wake, rotor wake, flight region). Bars: N production (blue), NE first release (orange), NE released SIMPLEC (aqua). Dashed lines mark the limits: 5% for RMSE, 15% for P95.
- **What to point at:** the NE building wake at 4.5–4.6% RMSE and 10.1–10.3% P95, against N's 2.2% and 4.5%. Switching NE to SIMPLEC barely changes the errors (4.49% to 4.56%), because relaxation changes the iteration path, not the answer.
- **Takeaway:** every zone passes with full coverage (100%). The weak spot is the building wake when the wind comes from NE.
- **Caveats:** zone sizes differ (rotor wake: 24,200 points; flight region: 2,961). This is not a map; see "Not built" below for why no point-by-point error map exists.

## fig08-strict-residuals-vs-limit.png: why strict references converge or not

- **Purpose:** show the convergence gate for every strict (reference-quality) solve, and why the tighter inner tolerance was needed.
- **How to read it:** each row is one strict solve. The x axis (log scale) is the final residual divided by its strict limit, where the residual is the leftover imbalance in the equations. Left of the dashed line passes. Blue circles: pressure (limit 1e-4). Orange diamonds: the largest velocity component (limit 1e-5). The right column gives iterations and whether the run converged or hit the 600 s cap.
- **What to point at:** for N, the original run stopped at 2.4× the pressure limit. Damping momentum made both residuals worse (9–13×). Tightening the inner solver tolerance (1e-7 to 1e-10) converged in 639 iterations. NE converged at its first strict attempt. Every intermediate or extra heading stayed 1.5–27× above a limit: the four NNE settings, NbE, and NW. The NNE SIMPLEC 0.9 run lowered pressure but pushed crosswind velocity to 27×.
- **Takeaway:** for N the problem was the inner linear solves skipping velocity updates, not the relaxation. The other headings stay above the limits within 600 s, even with that fix included.
- **Caveats:** these are final values only. Per-iteration residual histories are not stored; one trial records only three points. "N strict, original" is the value at iteration 995. Its exported field is from iteration 900. NNE, NbE and NW already use the tight inner tolerance.

## fig09-nne-residual-location.png: where the 22.5° solve is stuck

- **Purpose:** show that the NNE convergence failure sits at specific building edges, not at the turbine.
- **How to read it:** left panel: plan view of the district (grey footprints) and the rotor (black bar, turned to face 22.5° wind). Orange squares are 25 m blocks ranked highest for crosswind-momentum residual. The label gives each block column's share of the domain total and its height range. Orange diamonds are the three single cells with the largest crosswind residual. Blue circles are the two cells with the largest pressure residual, at a building–ground corner near (−168, 16, 1) m. Right panel: the share of the total crosswind residual in each group of cells.
- **What to point at:** two stacked blocks at one elevated building edge (25–75 m high) hold about 10.2% of the whole domain's imbalance. The 2,184 rotor-force cells hold 0.02%.
- **Takeaway:** the solver is stuck at one building edge and one building–ground corner. The turbine is not the cause.
- **Caveats:** the data come from short replays (three iterations) of the archived SIMPLEC run at iteration 902. The SIMPLE replay at 1002–1003 agrees on the pressure corner and the same edge. These are unfinished solves, so the figure locates the imbalance, not error in a final answer. Block shares sum the top-5 blocks by column only.

## fig10-what-failed-and-why.png: what was tried and the measured reason

- **Purpose:** one summary panel of all mesh and solver changes, each with its measured outcome.
- **How to read it:** one row per attempt: what changed, the measured result, and the outcome (✓ accepted or fixed, ✗ rejected or failed).
- **What to point at:** most mesh ideas failed on the strict "same rotor-force cells" rule. That includes a 0.65 µm vertex shift. The rule is deliberately strict, so the comparison isolates the mesh effect.
- **Takeaway:** two hypotheses were accepted (tighter tolerance for the N reference, post-snap merging for production). SIMPLEC relaxation sped NE up but left it at 1.99×.
- **Caveats:** a text panel, not a chart; the rows are typed into `src/fig_fail.py` from `joint/results/*` and `joint/results/optimization-summary.json`. The 0.55 mm audit row was a diagnostic false alarm: two different centroid estimators were compared, and the exact common-integral check passed to roundoff.

## Earlier slide images

- `niigata-map.png`: the Niigata district and the 80 wind-tunnel measurement points (district at 1:250 scale in the tunnel).
- `swift-sketch.png`: the SWiFT site: upstream met tower for the inflow, rear-facing lidar for the wake.
- `01-turbine-result.png`: the isolated-turbine result on three held-out conditions: 5.97–6.09× faster, 1.38–1.63% wake RMSE against the 5% limit. fig02 shows the same change on the mesh itself.
- `03-city-limit.png`: why the city alone was not certified: 6.42% (whole city), 7.24% (below 40 m) and 7.40% (downstream buildings) against a provisional reference that never converged.

These four predate the shared palette and are not 1920 px wide (they are 1800–2872 px).

---

## Not built, and why

- **Error map, |production − strict|:** the strict-reference grid samples (`north-reference.csv.gz`, `northeast-reference.csv.gz`) used by `joint/lookup/*-grid-check.json` are machine-local on Studio and not in the repository. Only the production lookup tables are stored. fig07 gives the per-zone numbers instead. City error by height for N is in `joint/results/poc-N-agglomerated-siblings/metrics.json` (largest at 20 m: 0.256 m/s, 2.9%).
- **Residual history vs iteration:** no per-iteration logs are stored. fig08 shows final residuals against their limits.
- **Merged-mesh velocity or a cheap-mesh 3D render:** no production volume is stored. fig03 and fig04 reconstruct the merged geometry only.
- **Per-failure figures** (city-floor coarsening, undeformed-only merge, pre-snap refinement, post-snap split): the records contain hashes, counts and log facts, not geometry, so fig10 summarises them.
