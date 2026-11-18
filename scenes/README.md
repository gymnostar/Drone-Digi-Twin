# Buildings and turbine wind study

**Question:** can we generate wind maps more cheaply while keeping their numerical error small? We tested the full Niigata building district and a SWiFT V27 turbine separately, then solved their interacting flow in one domain.

| Study | Result | Interpretation |
|---|---|---|
| Isolated turbine | About **6× faster** complete map generation; wake vector RMSE **1.38–1.63%** on three reserved numerical conditions | Successful numerical cost reduction within the tested setup |
| Buildings | **112,358 cells, 81.85 s**; city RMSE **6.42%**, above the 5% target | Cheap solve works, but accuracy target is unmet and stricter references did not converge |
| Buildings + turbine | **230,898 cells, 249.42 s**, versus **422,894 cells, 525.27 s** for the checked north-wind reference | **2.106×** cold speedup; all four regional velocity and coverage gates passed |

Errors compare calculations with a numerical reference, **not with real wind**. These are research results, not flight-safety or physical-validation claims.

### What we did

A wind map answers a practical question: **at this position, which way is the air moving and how fast?** We want that answer throughout the space a drone might use, not only at one weather station. A building can create shelter, speed up flow around an edge and leave a disturbed region downstream. A turbine removes momentum from the air and creates a wake (the region of altered flow behind it).

We cannot get those interactions by copying the incoming wind speed everywhere. We solve equations across the whole air volume. The pipeline is **geometry → mesh → solve → wind map → comparison**. A mesh divides the air into small cells. In each cell, the solver estimates velocity, pressure and turbulence. Neighbouring cells exchange information until pressure and velocity approximately satisfy conservation of mass (air does not disappear) and momentum (forces explain its motion).

An iteration is one numerical update of this solution. It is not one second of simulated flight. More cells mean more unknown values and more work per update. More iterations mean repeating that work more times. A detailed mesh may capture stronger changes near a wall or in a wake, but costs more; a coarse mesh may erase those changes. The central experiment is to find where we can reduce that detail without changing the outputs too much. We measured preparation, meshing, solving and extraction together (extraction means writing the predicted wind at the fixed positions used for comparison). Reading an existing map is a different operation and is not the reported speedup.

Before testing, we fixed provisional targets: at least **2× faster**, vector RMSE ≤**5%** of incoming speed, P95 ≤**15%**, and ≥**95%** coverage. RMSE (root mean square error) measures typical differences and penalises large ones. Vector error includes direction and speed. P95 is the error below which 95% of sampled errors fall, not a maximum or confidence interval. Coverage counts common valid samples against reference-valid points. These are proof-of-concept targets, not agreed drone tolerances.

Runs use three CPUs, 5 GiB, ten-minute mesh/solve limits and bounded repairs. The [frozen protocol](poc-protocol.json) records conditions, thresholds and budgets; it has not been edited retrospectively. Failed trials remain in the structured evidence.

### What the future lookup table needs

The planner will ask for a wind vector at a position and inflow condition. For the coupled scene, choosing a production mesh therefore requires two comparisons. First, check that each solve converges and that the cheaper recipe meets the four regional gates for the city, building wake, rotor wake and flight region. Then compare across mesh resolutions at the same physical positions. Changing solver settings on one mesh measures sensitivity to the solve and leaves spatial discretisation error open.

The coupled comparison also reports error and coverage at 5, 10, 20, 40 and 80 m above ground and in a narrow-gap proxy. The proxy selects valid city samples with a missing reference sample within 10 m on each side along either horizontal grid axis at the same height. It indicates where the sampling grid sees a constriction. The masks come from the reference, and these additional breakdowns leave the four original acceptance gates intact.

Once a mesh is supported, each wind direction needs its own field with world-frame velocity components in m/s, inflow metadata, rotor orientation and a valid-air mask. Reserve a whole direction to check interpolation or regression. Measure offline field-generation cost separately from repeated query latency. The current joint protocol covers one north-wind condition and sets no planner latency or vehicle-level error tolerance.

### Main implementation choices

- **OpenFOAM v2306, steady RANS:** the solver predicts average turbulent flow. We chose this affordable starting point to study mean wind maps before attempting expensive time-resolved turbulence. RANS models unresolved fluctuations rather than resolving individual gusts. SIMPLE/SIMPLEC are iteration methods coupling pressure and velocity.
- **Full building geometry:** Niigata after-construction geometry is scaled from its 1:250 experiment into full-size metres, without cropping. The numerical scenario has its own declared inflow; it does not recreate the tunnel experiment.
- **Neutral atmospheric boundary layer:** wind increases with height using a coherent inlet/ground treatment. The common development condition is 8.7 m/s at 32.1 m, roughness length 0.01158 m (a ground-drag parameter). Empty-domain runs check profile transport before adding obstacles.
- **One global realizableKE closure:** a closure supplies unresolved turbulent mixing. This stock model was selected from turbine development evidence and used across the coupled domain. Its physical superiority around buildings is unproven.
- **Actuator disk:** a 27 m rotor at 32.1 m hub height slows the air through a momentum sink (opposing force), without moving blades. Locked v3 rotor curves supply coefficients. Coupled forcing uses solved local inflow. Tower/nacelle drag, swirl and tilt are omitted.
- **Fixed force support and queries:** the same physical cells must receive rotor force across mesh comparisons; equal total thrust is insufficient. Building vicinity, rotor wake and proposed flight region are evaluated separately. Missing samples are not moved or filled.

### Turbine: successful bounded comparison

**The approximately sixfold speedup came from a cheaper CFD discretisation, not a neural network or a different computer.** We kept the same physical model and solved it on a mesh with less detail away from the force application region.

There are two different comparisons in this work. The early physical pilot asks whether a calculation resembles measurements. The cost experiment asks whether a cheaper calculation resembles a more detailed, checked calculation under the same declared inputs. The sixfold result belongs to the second comparison. Its denominator is the more detailed **same-physics numerical reference**, not the original experimental pilot. Earlier pilot timings used different settings and forcing support, so the result is not a claim that every new case is six times faster than every old case. In particular, the originally very small turbine pilot was already fast, but did not provide this controlled accuracy comparison. Changing the turbulence model during development is not what is being counted as the speedup: both timed recipes use realizableKE.

#### Exactly what changed

The turbine has a diameter of 27 m. The reference uses approximately **1.69 m cells in the refined wake region**, while the cheap recipe uses approximately **6.75 m cells there**. Refinement means dividing larger mesh cells into smaller ones. Here the cheap wake cells are four times wider in each direction. A single fully refined volume could therefore have many fewer cells, but the whole mesh does not shrink by that factor because other regions retain their resolution.

The total cell counts fall from **359,334 to 112,206**, about **3.2 times fewer cells**. We did not remove the rotor or simply remove difficult sample points.

Two important regions stay fine:

- The force is applied to the same **6,496 rotor-source cells**, approximately **0.844 m** wide, in both meshes. Their actual spatial geometry is checked, not just their number. Otherwise a supposedly cheaper mesh might change how the turbine slows the air, making the comparison unfair.
- A common refined near-ground band controls ground-related resolution and sampling differences. The ground treatment is retained, rather than making everything uniformly coarse.

The domain size, rotor coefficients, incoming wind, ground roughness, turbulence model and velocity transport scheme are matched. The fixed wake sample locations are matched too. Within this comparison, we changed the numerical resolution around the wake while protecting the forcing and near-ground treatment.

#### Why fewer cells made the whole calculation faster

With fewer cells there are fewer unknowns to update. In the development comparison, the reference also needed **464 iterations**, while the final cheap recipe needed **191** to meet the stopping checks. A coarser set of numerical equations can settle differently from a finer one; iteration count is an observed result, not a fixed guarantee for every condition.

Crucially, the final recipes have **the same residual stopping thresholds**. We did not obtain the accepted result merely by stopping the cheap run much earlier at a looser standard. The first cheap attempt did use looser stopping, and its wind at the monitored points was still changing too much even though its field error looked acceptable. We rejected that attempt, tightened stopping, and only then froze the recipe for the reserved conditions.

Fewer cells and fewer iterations both contribute. They should not be multiplied as independent theoretical speedup promises: iteration cost includes mesh connectivity, numerical behaviour and stage overhead. The evidence is the measured end-to-end timing:

| Development stage | Reference | Final cheap recipe | What the stage does |
|---|---:|---:|---|
| Preparation | 0.19 s | 0.19 s | Verify inputs and create case/query files |
| Mesh | 6.29 s | 2.65 s | Build and check the cells |
| Solve | 120.26 s | 17.58 s | Iterate the flow equations, including parallel setup and reconstruction |
| Extraction | 4.36 s | 1.64 s | Export wind at the fixed comparison positions |
| **Complete map generation** | **131.10 s** | **22.06 s** | **All four stages together** |

The speedup is **reference time / cheap time = 131.10 / 22.06 ≈ 5.94**. Most saved time comes from the solve, with additional savings in meshing and extraction. This is one measured timing per recipe, not a statistical timing study. The absolute seconds can change with the execution environment; the controlled comparison uses consistent resources.

#### Why we can accept the cheaper map within this study

Making a calculation faster is easy if its answer is allowed to change arbitrarily. We therefore compare the velocity vectors at fixed locations, particularly in the wake. A vector has direction and magnitude: two winds of equal speed but opposite direction are very different answers.

The reference is checked against its own numerical convergence, flow-balance and output-settling requirements. We then compare the cheap and reference maps on the same valid positions. A valid position is one that is actually sampled in fluid, not inside a solid or missing from the output. Points are not moved to make errors look smaller.

For example, an RMSE of 1.63% with an incoming speed of 8 m/s corresponds to roughly **0.13 m/s vector RMSE against that numerical reference**. It does not mean every point is within 0.13 m/s, and it does not establish that the real wind is predicted to that accuracy. P95 and coverage provide complementary checks on the error distribution and missingness.

The final recipe was chosen on development evidence. Three entire wind/roughness conditions were reserved (held out) before final testing. We applied the frozen recipe to them once without tuning it afterwards. This is stronger than testing neighbouring points from the same field, because the incoming conditions themselves change.

| Reserved condition: speed / roughness | Reference / cheap map time | Speedup | Wake RMSE / incoming speed |
|---|---|---|---|
| 8.0 m/s / 0.008 m | 130.82 / 21.47 s | 6.09× | 1.63% |
| 9.5 m/s / 0.015 m | 135.72 / 22.73 s | 5.97× | 1.38% |
| 8.7 m/s / 0.025 m | 140.82 / 23.38 s | 6.02× | 1.39% |

Each comparison has 19,800/19,800 common wake samples and passes convergence and spatial-force checks. P95 is 3.26–3.90%. An enlarged-domain check changes wake RMSE by 0.506% of incoming speed. This is limited boundary-sensitivity evidence, not proof of mesh independence. Authoritative reports in [turbine results](turbine/results/) end in `verified-comparison.json`.

### Why the same idea was harder for buildings

A turbine wake in steady RANS is relatively smooth away from the force region. A district has many walls, gaps, rooftops and separated flows (flow detaching from a surface). A cell size that is adequate in open air may be inadequate in a narrow passage or behind a building. Coarsening can change the local geometry represented by the mesh as well as the flow gradients.

There are two distinct gates. **Numerical convergence** asks whether a particular mesh's equations have been solved sufficiently. **Accuracy against a reference** asks whether that converged answer is close enough to the comparison field. A cheap solution can pass the first gate and fail the second, which is what the city results show. A strict reference can also fail to converge, in which case it cannot certify the cheaper result.

### Buildings: accuracy target not established

Three strict-reference attempts missed convergence limits. Residuals (remaining numerical equation imbalance) plateaued or oscillated even when overall flow balance looked good. Two cheap meshes also failed geometry checks. The final cheap mesh/solve passes its own checks, but RMSE against the provisional reference is 6.42% across the city, 7.24% in the lower city and 7.40% downstream. Whole-city coverage is 30,301/30,505 points.

No certified city speedup is reported. Unused direction holdouts remain untouched because development did not establish a valid reference/method. See the [city summary](buildings/results/poc-study-summary.json).

### What combining the scenes actually means

We place the buildings and rotor in one air domain and solve one velocity and pressure field. The buildings alter the velocity reaching the turbine, and the rotor force alters the downstream field. Adding a standalone turbine wake to a city wind map would miss this mutual interaction and can break conservation; we do not do that.

The rotor is still an actuator disk, not a blade-resolved turbine. This lets us investigate the mean momentum effect at manageable cost. The visible ring is an annotation of where the force acts, not a solid object that was added as a wall.

### Coupled flow and the selected north-wind recipe

The hypothetical rotor is downstream at (0, −280, 32.1) m, with a local-inflow monitor at (0, −215, 32.1) m and rotor-bottom clearance of 18.6 m. One solve computes building/turbine interaction; separate wind fields are never added.

Different mesh recipes slightly moved rotor-source vertices during snapping (fitting cells to surfaces). We rejected their unequal force support and used the same 422,894-cell mesh for the solver comparisons. All retained runs share the same 2,184 rotor-source cells, verified by their vertices. The mesh exceeds the original 150,000-cell cheap target and fits the declared 500,000-cell joint allowance.

The baseline converged in 364 iterations and passed its own cheaper stopping checks. Its recorded pipeline cost of 306.76 s excludes the original host force-support audit. The original strict reference timed out, leaving a field at iteration 900. The previously deferred SIMPLEC candidate completed in 377.65 s. Its field agreed closely with that provisional reference across all four regions. However, the four rotor-wake monitors changed by 0.247% of inflow over iterations 350–375, above the 0.20% settling limit, so we rejected it.

An invalid early trial caused by a dictionary-substitution bug was stopped and excluded; generated solver settings are now checked before execution. Failed implementation/support trials and costs remain in the [joint summary](joint/results/poc-study-summary.json).

The first **checked strict north-wind reference** took 525.27 s for complete map generation on Studio. The original strict solver timed out after 600 s of solving. Reducing momentum relaxation from 0.60 to 0.45 also timed out. Restoring 0.60 and tightening the inner linear-solver absolute tolerance for velocity and turbulence from `1e-7` to `1e-10` produced convergence in 639 iterations and 393.09 s of solving. The first changed velocity update appeared at iteration 436, where the tighter solver took an additional sweep. The final pressure residual was `9.93e-5` against a limit of `1e-4`. The final field and force were sampled at the same iteration and passed the flux, force and settling checks. [Strict reference metrics](joint/results/poc-N-reference-tight-transport/metrics.json) · [Experiment record](joint/results/optimization-summary.json)

A fresh baseline produced byte-identical samples to the historical run. With both host audits included, it took **317.55 s**, giving **1.65×** against the checked reference. At the 24,200 common rotor-wake points, vector RMSE was **0.0457 m/s** (0.525% of inflow) and P95 was **1.19%**. Coverage was 100% in each of the four gated regions. All 657 narrow-gap proxy points were also present, with RMSE of 0.0099 m/s. Each recipe has one complete timing. At the measured reference time, meeting the 2× goal requires a complete map in **262.64 s** or less. [Retimed baseline and regional errors](joint/results/poc-N-baseline-retimed/metrics.json)

The working recipe builds the same snapped mesh from the original geometry, then merges 27,428 interior groups of eight cells. Boundary polygons, the 2,184 rotor-source cells and every surviving face retain their coordinates. The resulting **230,898-cell** mesh reached the stopping criteria in **308 iterations**. Complete map generation took **249.42 s**, including compilation, a failed centroid diagnostic and its recovery. There was one timed run per recipe. This establishes a checked north-wind numerical comparison. The original 150,000-cell aspiration and physical validation remain open.

| Region | Common / reference-valid points | Vector RMSE | P95 vector error |
|---|---:|---:|---:|
| City | 30,504 / 30,504 | 0.194 m/s | 0.362 m/s |
| Building wake | 6,354 / 6,354 | 0.193 m/s | 0.391 m/s |
| Rotor wake | 24,200 / 24,200 | 0.117 m/s | 0.233 m/s |
| Flight region | 2,961 / 2,961 | 0.140 m/s | 0.293 m/s |

City error was largest at 20 m above ground among the five sampled heights (0.256 m/s RMSE). The 657 narrow-gap proxy points had 0.162 m/s RMSE. All regional percentages use the prescribed 8.7 m/s inflow. Standard and custom mesh checks passed. The extended geometry check flagged 41,904 concave cells under the existing acceptance rule. [Full field and mesh metrics](joint/results/poc-N-agglomerated-siblings/metrics.json)

Two earlier mesh trials stopped before solving. Coarsening the city floor moved the rotor-source vertices. Merging only undeformed Cartesian cells failed boundary verification. For the selected mesh, the initial centroid audit compared OpenFOAM estimates calculated from different pyramid decompositions. Measuring both geometries with the same face triangulation verified volume and centroid conservation to roundoff. The record retains the original discrepancy and charges both verification attempts. [Trial history, versions and commands](joint/results/optimization-summary.json)

To reproduce the selected pair on Studio with the pinned raw inputs and Docker image, run from the repository root with fresh case names.

```sh
set -e
study=scenes/joint/study.py
ref=north-strict
candidate=north-coarse
python3 "$study" build --case "$ref" --role reference --common-mesh --tight-transport
for stage in mesh solve sample evaluate; do python3 "$study" "$stage" --case "$ref"; done
python3 "$study" build --case "$candidate" --role baseline --common-mesh --agglomerate
for stage in mesh solve; do python3 "$study" "$stage" --case "$candidate" --reference "$ref"; done
python3 "$study" sample --case "$candidate"
python3 "$study" evaluate --case "$candidate" --reference "$ref"
```

### Physical evidence and limits

The original Niigata tunnel pilot has speed-ratio RMSE about 0.198 on 77/80 samples, with unresolved inlet/turbulence metadata and profile drift. The new full-scale neutral-inflow toy is a different scenario.

The initial SWiFT standard-k–epsilon comparison has measured wake RMSE about 1.11 m/s. Its observed periods already entered development, so historical date splits are exploratory. The three reserved conditions above test numerical generalisation, not independent physical observations. There are no joint interaction observations.

A numerical reference is a more demanding calculation, not physical truth. We have not demonstrated instantaneous gust accuracy, flight safety or a validated integrated wind/RF/drone system.

### What the current result is useful for

The turbine result demonstrates that a carefully chosen mesh can reduce repeated mean-flow map generation cost while preserving agreement with a checked numerical reference over the tested condition range. It gives us an economical recipe for a broader numerical campaign, subject to further validation.

It does not give an instant model for arbitrary buildings, turbine placements or weather. New geometry can change the required mesh and convergence behaviour. The city and coupled results are therefore useful boundary cases: they show where a successful isolated-object recipe stops transferring automatically.

Nor have we trained a surrogate (a compact model fitted to previously computed flow fields). Such a model could make predictions very quickly later, but generating its training database would still cost CFD time. That offline cost and genuinely new test conditions would have to be included in a separate experiment.

### From the north-wind case to a direction table

1. Retain the north-wind pair as the baseline for directional development. The two accepted resolutions quantify sensitivity at the sampled locations. A three-level refinement trend and directional transfer remain untested.
2. After choosing a working mesh, define the multi-direction campaign and reserve a whole direction before fitting. Save a valid-air mask and east, north and vertical velocity components for every field. The frozen coupled experiment covers a single north wind.
3. Compare lookup or regression predictions with the held-out CFD field, including the cost of generating the table and the latency of repeated queries. Physical validation still needs measurements under matched conditions. The standalone building study also needs its own checked reference.


### Open the results

Double-click `review-export/Open3D.command` for the coupled baseline. To select another scene, run that launcher with `turbine-cheap` or `turbine-reference`. No simulation is started by opening these files.

### Are these different files or different views?

There are **three separate saved simulation scenes**, not one run with three alternatives overlaid:

| Scene | What it contains | Role in the study |
|---|---|---|
| `joint-baseline` | Buildings and turbine in one solved flow | Coupled baseline; stricter reference remains provisional |
| `turbine-cheap` | Isolated turbine with the verified cheaper wake mesh | Accepted cheap development run |
| `turbine-reference` | Isolated turbine with the finer wake mesh | Checked numerical comparison for that development condition |

Each scene has a `.vtu` file (the actual mesh and computed fields) and a `.pvsm` file (ParaView's saved display state: camera, filters and colours). The coupled scene also uses `buildings.stl` for the surfaces. The state does not solve the equations; it arranges the saved data for viewing. The exported turbine scenes show the development condition, not all three reserved-condition results. The other runs and their metrics remain separate evidence.

**Within one scene, the pipeline items are different ways of showing that same result.** The two slices are sections through one three-dimensional velocity field. Streamlines are generated from that field. They are not extra CFD runs.

### What to show and hide in ParaView

Use the eye icon beside an item in **Pipeline Browser**. Selecting its name changes the settings shown in **Properties**; it does not necessarily make the item visible. After changing a filter setting, click **Apply** when it is enabled. Drag to rotate and scroll to zoom.

Start with **Simulated building surfaces**, **Ground**, **Rotor source marker**, **Simulation description**, and **Velocity at rotor height** visible. Hide the other items initially. For isolated-turbine scenes, there is no building-surface item.

| Pipeline item | Show/hide during presentation | What it does |
|---|---|---|
| **Saved CFD volume** | Hide | The full computed air mesh and its velocity `U`, pressure `p` and turbulence `k`. Keep the source loaded; showing its outer surface often hides everything inside. |
| **Velocity interpolation for display** | Hide | Converts cell values to point values so slices and streamlines can be displayed smoothly. This is display processing, not another solve or an accuracy improvement. |
| **Display region** | Hide | Clips the displayed data to a useful area near the obstacles. It does not crop the original simulation domain or change boundary conditions. |
| **Velocity at rotor height** | Show for a top/oblique overview | Horizontal section, initially at 32.1 m. Colours show mean speed in m/s. Select it and change **SliceType → Origin → Z** to inspect another height. This moves the viewing plane, not the turbine or inlet. |
| **Vertical wake section** | Show instead of the horizontal slice for a side view | Section through the wake centre: x=0 in the coupled case, y=0 in the isolated-turbine frame. It shows how the predicted wake varies with height and downstream distance. Hide the horizontal slice if the overlap is confusing. |
| **Steady mean-flow streamlines** | Usually hide; show **Flow paths** instead | Computes curves following the saved mean velocity, starting from seed points. Moving the seeds changes which curves are drawn; it does not change the wind. |
| **Flow paths** | Show when explaining direction; hide for a clear slice | Turns the streamline curves into visible tubes. They show the direction of steady mean flow, not instantaneous gusts, actual particle tracks or drone routes. |
| **Simulated building surfaces** | Show | Original building geometry used for this coupled case. The surfaces help locate shelter and wakes. Visibility does not remove walls from the already solved case. |
| **Ground** | Show for orientation; hide if it obstructs the slice | A display plane at ground level. The actual ground boundary condition belongs to the solver case, not this display plane. |
| **Actuator disk location** | Hide | Untransformed ring source used to construct the rotor marker. Showing it as well can put an extra ring at the wrong location. |
| **Rotor source marker** | Show | The transformed ring at the modelled rotor position. It marks where momentum was removed; it is not rotating blade geometry or a new force applied by ParaView. |
| **Simulation description** | Show, or hide for a cleaner view | A text label identifying the selected scene. |

Suggested explanation sequence: **geometry first** (buildings, ground and rotor marker), **horizontal slice next** (sheltered and faster regions), **vertical slice next** (wake through height), then **Flow paths** (mean-flow direction). Show one main diagnostic at a time so overlapping transparent surfaces do not obscure it. All colours use the same 0–12 m/s display range; changing the colour scale alters appearance, not the predicted wind.

### Which settings actually define a run?

The ParaView filters above are post-processing (looking at a calculation after it has finished). They do not set up the CFD run. The case builder and OpenFOAM dictionaries define the domain, geometry, mesh refinement, incoming wind, ground roughness, turbulence closure, rotor force and convergence limits. Changing those requires building/running a new case. Moving the camera, hiding a building, changing a slice or moving streamline seeds does not require CFD and does not change those physical assumptions.

The coupled scene uses a wind from north flowing toward negative y; the isolated turbine uses a wind-aligned x axis pointing downstream. Their camera orientations and coordinate axes therefore differ. Compare the numerical values and the downstream direction, rather than assuming the same screen direction means the same world direction.

[Scripts and metrics](./) · [3D viewer](review-export/Open3D.command) · [Frozen protocol](poc-protocol.json)
