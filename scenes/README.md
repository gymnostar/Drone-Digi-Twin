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

### Query a coupled wind field

The lookup table supplies steady mean wind at fixed positions. Each compressed CSV includes the inflow speed and direction, rotor orientation, source case, evaluation checksum and an explicit valid-air column. Positions are local east and north coordinates in metres, with height above flat ground. Velocity components point east, north and up, in m/s. Wind direction is the meteorological **from** angle. [Campaign settings](joint/lookup/campaign.json)

The common grid covers east −200 to 200 m and north −450 to 200 m at 5 m spacing. It samples heights of 5, 10, 20, 32.1, 40, 60 and 80 m. Every position has a row, including locations inside a building or missing from CFD sampling. The prescribed speed is 8.7 m/s at 32.1 m and the ground roughness length is 0.01158 m.

Run from the repository root with Python. This reads an existing field and requires no OpenFOAM installation.

```sh
python3 scenes/joint/lookup.py query N 0 -350 20
```

Here the north-wind wake has velocity **(0.015, −4.839, 0.056) m/s**. An exact grid point without a CFD sample returns `masked`. A between-grid or outside position returns `unsampled`. For repeated queries, load the field once.

```python
from scenes.joint.lookup import WindField

wind = WindField.for_direction("N")
u_east, u_north, u_up = wind.at(0, -350, 20)
```

The [field index](joint/lookup/index.json) lists released directions, their sampling roles, file checksums and measured export times. The Python interface returns `None` for a masked point and raises `KeyError` for an unsampled point. The masks describe sampled fluid cells at this resolution. A grid lookup does not test whether a path between two points crosses a wall.

For a short integration check, query the four waypoints below in order for each heading. All are at east −195 m and height 5 m. At the two endpoints the entries are (east, north, up) velocities in m/s, rounded to three decimals.

| North coordinate (m) | North wind | Northeast wind |
|---:|---|---|
| 5 | (−0.085, −1.374, −0.172) | (−0.392, −0.126, −0.103) |
| 10 | masked | masked |
| 11 | unsampled | unsampled |
| 25 | (−1.030, −3.677, −0.295) | (−2.933, −0.870, −0.316) |

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

The coupled directional experiment below uses two independently generated fields to fit a per-position wind estimate. A simulation at a third heading tests how much of the change between directions that estimate captures.

### What transferred to northeast wind

The cheaper mesh passed the original numerical gates and the four regional checks on the exported table for both training directions. The northeast building wake was closest to the 5% RMSE limit, at **4.49% of inflow** (0.391 m/s). All 5,069 reference-valid points there had a candidate value. On the common 74,277-point grid, the masks agreed with the corresponding strict reference at every position for each direction. [North grid comparison](joint/lookup/N-grid-check.json) · [Northeast grid comparison](joint/lookup/NE-grid-check.json)

| Wind from | Strict / production cells | Strict / production cold CFD time | Speedup |
|---|---:|---:|---:|
| N | 422,894 / 230,898 | 525.27 / 249.42 s | 2.11× |
| NE | 425,006 / 228,432 | 401.14 / 237.59 s | 1.69× |

Each entry is one complete generation timing on Studio using three CPUs and 5 GiB. The northeast case missed the original twofold cost target. Preparing the additional common-grid CSV from the finished production field took 9.59 s for north and 9.06 s for northeast in one export each. For northeast, 34,253 valid overlapping rows agreed exactly with the original regional samples. [North CFD comparison](joint/results/poc-N-agglomerated-siblings/metrics.json) · [Northeast CFD comparison](joint/results/table-NE-production/metrics.json)

The first NE mesh exceeded the 500,000-cell limit because the solver-axis bounds enlarged the refined footprint. Rotating the refinement boxes with the city brought the strict mesh to 425,006 cells. The first NNE mesh failed the extended checks at two building corners. Local refinement removed the low-determinant cells and low-weight face. Exact rotor-source geometry matches across the coarse and strict meshes **within** each training direction. Between headings, the disk selects 2,184 cells for N and 1,810 for NE as its position shifts relative to the solver lattice. This contributes numerical uncertainty to directional wake predictions. [Rejected NE mesh](joint/results/table-NE-strict/mesh-verdict.json) · [Rejected NNE mesh](joint/results/table-NNE-oriented/mesh-verdict.json)

### Estimating wind between two headings

The experimental model fits a straight line through the N and NE vectors at each fixed world position. It expresses each vector along and across its incoming wind, blends those components and the vertical velocity, then rotates the result into the requested world heading. This preserves a uniform inflow's speed through a turn. A direct average of equal-speed world vectors 45° apart would lose 7.6% of their speed at the midpoint. The estimate requires valid samples at both endpoints and accepts angles from 0° through 45°. For an intermediate heading the command returns `estimated` or `masked_endpoint`. Building edges can fall in different sampled cells as the mesh turns.

```sh
python3 scenes/joint/lookup.py estimate scenes/joint/lookup/N.csv.gz scenes/joint/lookup/NE.csv.gz 22.5 0 -350 20
```

We reserved the entire **22.5° NNE direction** before fitting. Its strict SIMPLEC solve reached the ten-minute limit with pressure and velocity residuals above the frozen thresholds. A second cold run changed only the pressure–velocity coupling to standard SIMPLE on a byte-identical mesh and also timed out. Both attempts and their costs remain recorded. NNE has no released CFD field or qualifying accuracy score. [SIMPLEC run](joint/results/table-NNE-corners/solve-verdict.json) · [Standard SIMPLE run](joint/results/table-NNE-simple/solve-verdict.json)

The saved checkpoints still help locate the discrepancies. At iteration 900, the two solvers differed by 0.014 m/s vector RMSE in the building wake. The estimate differed from each checkpoint by 1.76 m/s there. Corresponding differences in the flight region were 0.003 m/s between solvers and 0.66 m/s between the estimate and either checkpoint. In the narrow-gap proxy only 564 of 626 checkpoint-valid locations had valid samples in both training fields. These figures describe unfinished iterates. The formal RMSE target is 0.435 m/s against a checked reference. [Checkpoint diagnostics by region and height](joint/lookup/NNE-diagnostic.json)

After the midpoint timeouts, we reserved 11.25° (north by east) as a secondary test and kept the fit fixed. One local refinement repaired its mesh. The strict solve completed 1,021 iterations within the time limit, with final pressure and x-velocity residuals of 2.29 × 10⁻⁴ and 1.48 × 10⁻⁵, above their limits of 10⁻⁴ and 10⁻⁵. Neither reserved heading produced a qualifying field. [Initial mesh](joint/results/table-NbE-strict/mesh-verdict.json) · [Secondary solve](joint/results/table-NbE-corner/solve-verdict.json)

For repeated queries in the first planner experiment, use the two released CFD headings directly. A warm direct query took 0.15 µs and a wind-aligned estimate took 0.72 µs on Studio (median batch mean over five batches of 10,000 calls). Loading each field took about 0.10–0.12 s. The campaign reserves 16.875° to compare the original fit with a split interval if a 22.5° reference qualifies. Each method will be scored on the same reference-valid locations in each region. [Query workload and timings](joint/lookup/query-timing.json)

Studio can reproduce the diagnostic from the two saved parallel cases.

```sh
python3 scenes/joint/checkpoint_diagnostic.py scenes/joint/runs/table-NNE-corners \
  scenes/joint/runs/table-NNE-simple --common-iteration 900
```

To locate the remaining equation imbalance, we replayed three iterations of each archived NNE solve and recorded the per-cell pressure and crosswind-momentum residuals. In the SIMPLEC replay, two adjoining 25 m cubes at one elevated building edge contained about 10% of the domain's absolute crosswind-momentum residual. Both solvers placed their largest individual pressure residual at the same building–ground junction, near world (−168, 16, 1) m. Across the four sampled iterates, cells in the rotor force region contributed at most 0.021% of the crosswind-momentum residual. These fractions measure equation imbalance in unfinished solves. [Per-cell residual results](joint/results/nne-residual-diagnostic.json)

```sh
python3 scenes/joint/residual_diagnostic.py scenes/joint/runs/table-NNE-corners --checkpoint 900
python3 scenes/joint/residual_diagnostic.py scenes/joint/runs/table-NNE-simple --checkpoint 1000
```

### Directional simulation settings

A new build defaults to the cheaper production recipe for N and NE, and to the strict reference for the reserved NNE and NbE directions. Use `--recipe strict` or `--recipe production` to select a comparison explicitly. Both use the same rotor support for a given direction. The building geometry and city refinement boxes rotate together into the solver frame. The world hub remains at (0, −280, 32.1) m and the rotor yaws with the inflow. Fixed world positions make direct field comparisons possible.

```sh
set -e
study=scenes/joint/study.py
python3 "$study" build --case ne-reference --direction NE --recipe strict
for stage in mesh solve sample evaluate; do python3 "$study" "$stage" --case ne-reference; done
python3 "$study" build --case ne-candidate --direction NE --recipe production
for stage in mesh solve sample evaluate; do
  python3 "$study" "$stage" --case ne-candidate --reference ne-reference
done
python3 scenes/joint/lookup.py export --case scenes/joint/runs/ne-candidate \
  --results scenes/joint/results/ne-candidate --output /tmp/NE.csv.gz
```

The exporter samples the final reconstructed field on the common grid. It requires a checked strict reference or a candidate that passes the matched-reference velocity, coverage and support checks. Re-extraction verifies overlapping points against the original regional samples. The two original north-wind resolutions and the directional pilot leave a three-level refinement trend and validation against joint observations open.

### Open the results

Double-click `review-export/Open3D.command` for the coupled baseline. To select another scene, run that launcher with `turbine-cheap` or `turbine-reference`. No simulation is started by opening these files.

### Are these different files or different views?

There are **three separate saved simulation scenes**, not one run with three alternatives overlaid:

| Scene | What it contains | Role in the study |
|---|---|---|
| `joint-baseline` | Buildings and turbine in one solved flow | Historical coupled baseline used for the saved 3D view |
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
