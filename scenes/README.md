# Buildings and turbine wind study

**Question:** can we generate wind maps more cheaply while keeping their numerical error small? We tested the full Niigata building district and a SWiFT V27 turbine separately, then solved their interacting flow in one domain.

| Study | Result | Interpretation |
|---|---|---|
| Isolated turbine | About **6× faster** complete map generation; wake vector RMSE **1.38–1.63%** on three reserved numerical conditions | Successful numerical cost reduction within the tested setup |
| Buildings | **112,358 cells, 81.85 s**; city RMSE **6.42%**, above the 5% target | Cheap solve works, but accuracy target is unmet and stricter references did not converge |
| Buildings + turbine | **422,894 cells**; baseline passes its numerical checks; stricter reference timed out | Actual interacting flow, but no certified joint speedup; candidate mesh ready, solve pending |

Errors compare calculations with a numerical reference, **not with real wind**. These are research results, not flight-safety or physical-validation claims.

### What we did

A wind map answers a practical question: **at this position, which way is the air moving and how fast?** We want that answer throughout the space a drone might use, not only at one weather station. A building can create shelter, speed up flow around an edge and leave a disturbed region downstream. A turbine removes momentum from the air and creates a wake (the region of altered flow behind it).

We cannot get those interactions by copying the incoming wind speed everywhere. We solve equations across the whole air volume. The pipeline is **geometry → mesh → solve → wind map → comparison**. A mesh divides the air into small cells. In each cell, the solver estimates velocity, pressure and turbulence. Neighbouring cells exchange information until pressure and velocity approximately satisfy conservation of mass (air does not disappear) and momentum (forces explain its motion).

An iteration is one numerical update of this solution. It is not one second of simulated flight. More cells mean more unknown values and more work per update. More iterations mean repeating that work more times. A detailed mesh may capture stronger changes near a wall or in a wake, but costs more; a coarse mesh may erase those changes. The central experiment is to find where we can reduce that detail without changing the outputs too much. We measured preparation, meshing, solving and extraction together (extraction means writing the predicted wind at the fixed positions used for comparison). Reading an existing map is a different operation and is not the reported speedup.

Before testing, we fixed provisional targets: at least **2× faster**, vector RMSE ≤**5%** of incoming speed, P95 ≤**15%**, and ≥**95%** coverage. RMSE (root mean square error) measures typical differences and penalises large ones. Vector error includes direction and speed. P95 is the error below which 95% of sampled errors fall, not a maximum or confidence interval. Coverage counts common valid samples against reference-valid points. These are proof-of-concept targets, not agreed drone tolerances.

Runs use three CPUs, 5 GiB, ten-minute mesh/solve limits and bounded repairs. The [frozen protocol](../scene-selection/poc-protocol.json) records conditions, thresholds and budgets; it has not been edited retrospectively. Failed trials remain in the structured evidence.

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

### Coupled flow: actual interaction, unfinished optimization

The hypothetical rotor is downstream at (0, −280, 32.1) m, with a local-inflow monitor at (0, −215, 32.1) m and rotor-bottom clearance of 18.6 m. One solve computes building/turbine interaction; separate wind fields are never added.

Different mesh recipes slightly moved rotor-source vertices during snapping (fitting cells to surfaces). We rejected their unequal force support and used the identical reference mesh recipe for baseline/candidate. **Mesh-cost reduction failed to transfer; the remaining comparison tests solver/stopping settings only.** The full joint mesh exceeds the initial cheap-cell target; its declared cap adjustment is recorded.

The baseline converged in 364 iterations. Its recorded pipeline cost is 306.76 s, a lower bound because original host force-support verification overhead was not timed. The strict reference timed out at 600.20 s; its saved iteration-900 field is provisional. Baseline/reference differences are small (rotor-wake RMSE 0.512%, full common coverage), but cannot certify reference accuracy or speedup. The candidate mesh has matching force support; its solve is pending.

An invalid early trial caused by a dictionary-substitution bug was stopped and excluded; generated solver settings are now checked before execution. Failed implementation/support trials and costs remain in the [joint summary](joint/results/poc-study-summary.json).

### Physical evidence and limits

The original Niigata tunnel pilot has speed-ratio RMSE about 0.198 on 77/80 samples, with unresolved inlet/turbulence metadata and profile drift. The new full-scale neutral-inflow toy is a different scenario.

The initial SWiFT standard-k–epsilon comparison has measured wake RMSE about 1.11 m/s. Its observed periods already entered development, so historical date splits are exploratory. The three reserved conditions above test numerical generalisation, not independent physical observations. There are no joint interaction observations.

A numerical reference is a more demanding calculation, not physical truth. We have not demonstrated instantaneous gust accuracy, flight safety or a validated integrated wind/RF/drone system.

### What the current result is useful for

The turbine result demonstrates that a carefully chosen mesh can reduce repeated mean-flow map generation cost while preserving agreement with a checked numerical reference over the tested condition range. It gives us an economical recipe for a broader numerical campaign, subject to further validation.

It does not give an instant model for arbitrary buildings, turbine placements or weather. New geometry can change the required mesh and convergence behaviour. The city and coupled results are therefore useful boundary cases: they show where a successful isolated-object recipe stops transferring automatically.

Nor have we trained a surrogate (a compact model fitted to previously computed flow fields). Such a model could make predictions very quickly later, but generating its training database would still cost CFD time. That offline cost and genuinely new test conditions would have to be included in a separate experiment.

### Next steps and pace

1. Finish the one remaining capped candidate solve. This alone will not fix the failed reference.
2. Establish a checked city/joint reference through a focused numerical study: revised steady setup, justified longer budget, or time-dependent flow. Do not loosen thresholds just to pass.
3. Obtain better matched inflow metadata or independent observations. Agree which wind quantities and errors matter for drone planning.
4. Expand turbine conditions first: comparable cheap cases took 21–23 s versus 131–141 s for references. This is measured throughput within this setup, not guaranteed scaling for new geometry. City/joint expansion is gated by the reference problem.
5. Then expand terrain, turbine count or fidelity, or build a surrogate (a smaller model approximating simulations). Include database/training cost. More compute increases throughput but does not fix modelling gaps. Wind/RF exchange and drone planning follow; reinforcement learning is optional.

The project has a reproducible pipeline and a successful isolated-turbine result. City reference convergence and physical validation remain the main blockers. No reliable completion date for those repairs follows from these measurements.


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

[Scripts and metrics](./) · [3D viewer](../review-export/Open3D.command) · [Frozen protocol](../scene-selection/poc-protocol.json)
