# Fröttmaning–Freimann: Project 3

## Project outline: goal and route to completion

**Goal:** build a useful wind model around a representative turbine, check how accurate it is, then make a cheaper model reproduce the important wake behaviour (the slower, disturbed airflow behind the rotor). Deliver reusable geometry, repeatable simulations and a clear accuracy-versus-computing-cost comparison. Start with one turbine under known conditions; add site complexity only when it helps answer that question.

The professor's **root README is the main reference**: one turbine, terrain and buildings, OpenFOAM with an actuator disk first, and a wider digital twin for drones later. The older **Project 3 proposal** adds the research sequence: choose a trusted reference, measure wake-model errors, correct a cheaper model, and test new conditions. Its requirements are recorded below; this document does not need the PDF to be understood.

The current phase is **Project 3 only**. RF propagation, swarm coordination, autonomous route planning and Unreal/Project AirSim remain later integration work from the README. Reinforcement learning (an agent learning through rewards) is also later; it is not needed to establish a valid wind model or a first correction. QGIS currently performs the geographic assembly that the original plan assigned to Fusion 360.

| Order | What needs to be done | Programs and data | Result / decision before continuing |
|---|---|---|---|
| 1 | Agree on the comparison and required accuracy | Professor; measurements, a published turbine experiment, or already validated CFD | A named reference with known wind, turbulence, turbine and sampling conditions. Munich has no measured reference in this repo. A simpler published case may be needed first. |
| 2 | Check the geometry needed by the chosen case | Benchmark geometry, or the existing Munich terrain/structures for a site case; Python and optional QGIS inspection | Correct dimensions, position and height conventions. The full Munich scene is reusable context; multiple QGIS views are not a prerequisite for validation. |
| 3 | Define one single-turbine case and its wind conditions | OpenFOAM case settings; the chosen reference data and geometry | Choose a computational domain (the volume of air being simulated), incoming wind and rotor orientation. No separate QGIS crop is needed. Use a simpler benchmark first if that provides the trustworthy comparison. |
| 4 | Prepare a separate CFD geometry copy | Geometry scripts; terrain/building surfaces | A surface mesh (triangles describing solid boundaries), with recorded repairs and no unintended holes. This preparation still needs rebuilding after the old archive was deleted. |
| 5 | Define the turbine and incoming atmosphere | OpenFOAM dictionaries (text settings); rotor dimensions, thrust/power coefficients, wind profile, roughness and turbulence | An actuator disk (a region applying a braking force to air). Use documented coefficients or benchmark values; the printable blades are not needed for this model. |
| 6 | Build and check the air mesh | OpenFOAM `blockMesh`, `snappyHexMesh`, `checkMesh`; local refinement near the rotor, wake and obstacles | A volume mesh (cells filling the air), checked for shape and connectivity. Choose cell sizes from accuracy needs, not the terrain's 1 m spacing. |
| 7 | Calculate a reference-quality mean flow | OpenFOAM `simpleFoam` with a justified turbulence model; ParaView for inspection | Stable residuals (equation error measures), stable probes (sample points), reasonable mass balance, and acceptable mesh checks. These numerical checks do not replace physical validation. |
| 8 | Check mesh and domain dependence | Repeat with finer cells and farther boundaries, holding the physical conditions fixed | Important wake results change little under the agreed tolerance. A more expensive calculation is not automatically a trustworthy reference. Compare with the external reference from step 1. |
| 9 | Establish a cheaper baseline | Proposed first choice: coarser steady CFD with the same geometry/conditions. An analytical wake model (a formula) is an alternative if its terrain/building limits are acceptable. | Measure accuracy, runtime and memory before adding a correction. Choose the baseline with the professor; neither alternative is implemented as a new study yet. |
| 10 | Fit the smallest useful correction | Python; matched reference/baseline cases | Begin with calibration factors or regression (fitting a relationship). Try a more complex model only if a measured error needs it. Keep the original uncorrected baseline. |
| 11 | Test generalization and report | Separate wind/yaw/speed conditions; ParaView extraction; Python/Matplotlib plots | Test whole conditions not used for fitting. Report where correction helps or fails, its cost, and its limits. Keep configurations, versions and any random seeds. |
| 12 | Package the result for later drone integration | Geometry/data, case-generation/run scripts, wind arrays, comparison plots and notes | A repeatable Project 3 result: validated/benchmarked wake model, cheaper baseline, correction, accuracy/cost comparison and limits. Drone planning can then use it as an input; it is a separate task. |

### What to measure

Compare **velocity deficit** (loss of wind speed behind the rotor), **wake width**, **turbulence intensity** (fluctuation level relative to a stated mean wind speed), and **downstream recovery** (how the wind returns towards its upstream state). Compare the same positions and heights under the same conditions. Positions such as 2D, 5D and 10D downstream are possible starting points; D means rotor diameter, and the reference dataset should determine the actual comparison points.

For deficit, state the upstream reference speed and use the streamwise velocity (component along the incoming wind), rather than silently substituting total speed. State error units and normalization. Match turbulence definitions across sources. Record mesh cell count, meshing time, solver time, post-processing time and peak memory separately. Choose accuracy thresholds before selecting a preferred result.

Keep fitting and test conditions separate; randomly splitting neighbouring samples from one wind field would not demonstrate performance on new weather. A steady RANS model (average flow with modelled turbulence) does not supply individual gusts or a time series of blade motion.

### Implementation choices still to agree

The README's actuator-disk approach is the starting point. `kOmegaSST` (the pilot's turbulence model) is a candidate, not a validated choice for this site. Atmospheric inflow (wind entering the domain), ground roughness and turbine coefficients need sources. Do not add detailed moving blades, LES (a more expensive time-dependent turbulence method), or learning before identifying a concrete need.

Position, wind direction, **yaw** (the nacelle/rotor's facing direction), turbine operating state and blade motion are distinct. A future case must set them consistently. For Project 3 they can be controlled scenario inputs; live measurement and drone state estimation are outside the current phase. Fixed buildings and hills mean that an existing west-wind field cannot simply be rotated to represent another wind direction.

## What has been done

- Assembled Fröttmaning–Freimann from Bavarian DGM1 (ground heights), LoD2 (simplified 3D structures) and local official building outlines. No DLR geometry was mixed into Munich.
- Expanded the terrain to 4 km × 3 km, preserving original coordinates and heights. Packaged a QGIS scene independent of the raw downloads.
- Retained the source towers and added illustrative rotors using recorded turbine dimensions. Styled roofs, bridges, edges and hillshade (ground shading) for inspection.
- Ran one preliminary E-66 wind calculation, then retired it while reviewing the scope. Its inspection copy remains in `wind-review/`.
- Kept one full-area QGIS view and one static site-reference image. Interactive crop variants were removed because they were optional presentation work, not a Project 3 requirement.

The old external archive was deleted by Fausto. Raw downloads, the complete old OpenFOAM case, old preparation scripts and detailed logs are no longer available here. The retained QGIS data and wind-review fields still open, but the old solve cannot be resumed or fully reproduced from those review files alone. No new wind solve has been run during this scope review. The restored full-area project was checked from a separate folder: all eight layers and the original 3D renderers load, with one saved 4 km × 3 km view. This is a file/display check, not wind validation.

## Geographic scene: dimensions and assumptions

| Item | Current representation |
|---|---|
| Ground | 12 km², original 1 m samples, no vertical exaggeration |
| Structures | 2,983 complete displayed objects, including 53 bridges; 52 objects crossing the ground boundary are hidden |
| Coordinates | EPSG:25832 (east/north in metres); DHHN2016 (national height reference) |
| Fröttmaning E-66 | Rotor diameter 66 m; specified hub height 66.8 m above ground |
| Freimann E-138 | Rotor diameter 138 m; specified hub height 80 m above ground |

Turbine XY positions come from the source tower models, not a new survey. QGIS rotor shapes and their fixed north-facing pose are illustrative. They are not the professor's printable model or verified manufacturer blades. The old CFD disk used a different orientation, aligned with its west-to-east flow.

The display uses LoD2 tower-base heights. At Fröttmaning the recorded LoD2/DGM1 bases differ by about 0.13 m; at Freimann by about 2.09 m. Reconcile the base and hub-height convention before a new simulation; do not silently lift all source objects. There are also 26 terrain/structure contacts to review. These flags are not 26 proven errors. The scene is not a validated collision model.

## Computational domain, cost and optional site context

Project 3 needs a single-turbine case, a trusted comparison, a cheaper baseline and a tested correction. It does not require multiple geographic scenes or an interactive QGIS crop system. Comparing different wind speeds or yaw states changes the simulation conditions; it does not require new QGIS visualizations.

A computational domain (the air volume being solved) is still necessary. Comparing smaller/larger domains can reveal boundary effects (distortion caused by nearby box edges). Configure that comparison in the CFD case, keeping local cell sizes and physics fixed. It is separate from making presentation views in QGIS.

Cost depends on air-cell count, refinement (smaller cells where needed), physical model and solver settings. The old case had 196,568 cells; its 6 GB/four-CPU container cap was a chosen budget, not the machine's maximum. Halving cell width throughout 3D space gives roughly eight times as many cells. Measure time and memory on the first checked case before making laptop-capacity claims.

The [static site-reference image](scene/crop-reference.png) shows possible areas around each turbine and both together. Its scene1–scene6 labels belong only to that image; they are not active QGIS scenes, approved CFD domains or a current work plan. Missing northern coverage is marked for possible future use. No extra data is needed now just to maintain these presentation ideas.

Future drone work can reuse terrain, obstacle geometry, turbine positions and numeric wind fields. QGIS is a useful preparation/inspection tool, not a required part of the drone simulator or controller. Keep the geographic data without maintaining several interactive crop views.

## What the retired wind pilot established

It demonstrated a geometry-to-solver-to-viewer workflow for one E-66 case: 229 candidate structures, 196,568 air cells, west wind at 5.4 m/s, steady `simpleFoam` / `kOmegaSST`, and an actuator disk of diameter 66 m and thickness 10 m. The domain top was 1,000 m in the national height system. Add `[696000, 5343500, 500]` to local output coordinates to recover the geographic frame. Ct=0.75 and Cp=0.5625 (thrust and power coefficients) were provisional test values, not a measured Enercon performance curve.

### Why a mesh is needed

A **mesh** divides geometry or space into small connected pieces that a computer can work with. We used two different kinds:

- **Surface mesh:** triangles covering terrain, buildings and the tower, telling the solver where air meets solids. The derived solids were joined and foundations extended for meshing; the original QGIS geometry stayed separate.
- **Volume mesh:** cells filling the surrounding air. OpenFOAM uses the finite-volume method (balancing mass and momentum entering/leaving each cell) to calculate wind. These cells are mainly hexahedral (box-like), with more complex shapes near surfaces; surface triangles do not mean the whole air volume consists of triangles.

`blockMesh` made the initial grid, `snappyHexMesh` refined it and fitted it to the solid boundaries, and `checkMesh` inspected cell quality. The pilot's background spacing was about 20 m, with smaller cells locally. DGM1's 1 m ground spacing and the README's proposed 2 m CFD spacing are separate from the actual air mesh. A finer mesh can capture sharper changes but costs more; refinement and validation determine the useful resolution.

### What the turbulence model does

**CFD** means computational fluid dynamics (calculating fluid motion). The pilot used OpenCFD OpenFOAM **v2306**, with `simpleFoam` for steady, incompressible flow (a time-averaged solution with constant air density). **SIMPLE** repeatedly adjusts pressure and velocity until the equations agree sufficiently. The saved iteration 750 is an iteration count, not 750 seconds of weather. [Solver reference](https://doc.openfoam.com/2306/tools/processing/solvers/rtm/incompressible/simpleFoam/).

Turbulence contains many swirling motions at different scales. Resolving all of them is too expensive here. **RANS** (Reynolds-averaged Navier–Stokes) solves the average flow and models their effect. Our **`kOmegaSST`** model solves two extra quantities: **k** (turbulent kinetic energy, describing fluctuation strength, in m²/s²) and **omega** (a turbulence decay-rate scale, in 1/s). **SST** means shear stress transport; it blends near-wall `k–omega` behaviour with `k–epsilon` behaviour away from walls. This estimates turbulent mixing; it does not generate a sequence of individual gusts. Suitability for our wake still needs checking. [v2306 model reference](https://api.openfoam.com/2306/classFoam_1_1RASModels_1_1kOmegaSST.html).

| Other pilot choice | Meaning and limitation |
|---|---|
| Actuator disk | Applies a braking force across the rotor region, producing a slower wake. The turbine extracts energy; it does not blow a jet of air. No rotating blades or blade-tip vortices are resolved. Ct/Cp and disk orientation must match the intended operating condition. |
| Uniform inflow; smooth-wall functions; no added boundary layers | Boundary conditions specify how air enters/exits and meets solids. Wind entered at one speed rather than a height-dependent atmospheric profile. Wall functions approximate near-surface flow; no thin cell layers were added there. Future ground roughness, inflow and **y+** (dimensionless distance of the first cell from a wall) need a consistent treatment. |
| Bounded upwind discretization | A robust first-order way to transport quantities between cells. It introduces numerical diffusion (artificial smoothing), which can blur wake gradients. Finer cells or higher-order schemes need checks rather than an automatic accuracy claim. |
| Fields and viewers | **U** is the 3D mean velocity vector; **p** is kinematic pressure (pressure divided by density, in m²/s²); **k/omega** describe modelled turbulence. ParaView inspects the full field with slices, arrows and probes; Matplotlib makes comparison charts from extracted numbers. Neither viewer computes the wind. |

The old figure showed one slice of a 3D result. Airflow arrows are not drone paths, and low mean speed behind a building does not establish safe shelter: turbulence, clearance and drone limits also matter. The historical extended mesh check flagged 54 underdetermined and 6,704 concave cells (poor cell shapes). Numerical settling did not establish physical accuracy: there was no measurement/benchmark comparison, mesh-independence study or domain-independence study. The pilot is therefore an inspection example, not a validated reference for fitting a correction. No planner or live yaw estimation was connected.

### Options for the detailed and cheaper wind fields

These are possible study choices, not additional methods already implemented. **Detailed** and **cheap** describe a comparison under stated conditions; neither label guarantees accuracy.

| Method | What it could contribute / main limitation |
|---|---|
| Refined steady RANS + actuator disk | Practical candidate for the detailed mean-flow case, retaining terrain/building interactions. Test mesh/domain dependence and benchmark agreement first. |
| Coarser steady RANS + actuator disk | Proposed cheaper CFD baseline: fewer cells, initially the same geometry, physical inputs and turbulence model. This makes discretization errors easier to compare; excessive coarsening can remove important obstacles. |
| Another RANS closure, such as realizable `k–epsilon` | An alternative turbulence assumption to benchmark, not an automatic speed or accuracy improvement. Changing closure and grid together makes causes of error harder to separate. |
| URANS / LES | URANS (unsteady RANS, e.g. with `pimpleFoam`) allows evolving average flow while still modelling turbulence. LES (large eddy simulation) resolves larger turbulent motions and models smaller ones; it needs finer cells, time steps and sufficient averaging time. Neither was used in the pilot. |
| Analytical wake formulas | Jensen (uniform deficit across the wake) or Gaussian (smooth deficit profile), available in tools such as FLORIS, can give very cheap baselines. These are not CFD and do not automatically resolve flow around our hill and individual buildings. [Model examples](https://natlabrockies.github.io/floris/wake_models.html). |
| Calibration / surrogate | Correct a measured baseline error, or fit a model that approximates reference outputs. Cheap evaluation follows an initial data-generation cost; test complete unseen wind/yaw conditions. No correction or surrogate has been trained here. |

Turbine detail is a separate choice from turbulence modelling: an **actuator line** applies forces along moving blade lines using blade/aerodynamic data; **blade-resolved CFD** meshes the actual rotating blades and needs much more geometry and computation. Either can be paired with different flow models. They are later options if the required wake behaviour cannot be captured with a disk, not prerequisites for starting Project 3.

The proposed first comparison remains one checked steady RANS/disk case versus a coarser version, followed by the smallest useful correction. Extra realism should address a demonstrated error. A simpler published benchmark may establish this method before using the Munich site.

## Adding RF propagation maps later

**RF propagation** (how radio signals travel) is part of the wider README goal, separate from the current Project 3 wind comparison. Wind maps describe airflow; RF maps describe the radio channel between specified transmitters and possible receiver locations. No Munich RF simulation has been run.

A possible sequence is:

1. **Define the radio task:** communication coverage, RF positioning, or both. Specify transmitter positions, frequency, power, antenna patterns (directions of stronger/weaker transmission), orientation, polarization, bandwidth and receiver assumptions. Example frequencies must be chosen for the intended system, not inferred from the geographic scene.
2. **Prepare an RF surface scene:** reuse terrain/buildings in metres with the same local origin and height reference as CFD. Export triangles in a Sionna-compatible scene format and assign electromagnetic materials: permittivity (response to an electric field), conductivity, thickness and any scattering parameters. LoD2 shapes and QGIS colours do not provide these properties. Start with documented material assumptions and test their influence. This uses surface geometry, not OpenFOAM's air-cell mesh. [Material reference](https://nvlabs.github.io/sionna/rt/api/radio_materials.html).
3. **Compute a baseline and then ray tracing:** free-space loss gives a simple unobstructed reference; a line-of-sight/obstruction model adds coarse blockage. Sionna RT traces paths interacting with surfaces. Record enabled reflection, scattering and diffraction (bending around edges), interaction limits, ray samples and seeds; support/settings depend on the chosen version. Increasing samples reduces sampling noise but costs time/memory.
4. **Export maps at several heights:** `RadioMapSolver` can produce path gain per transmitter on horizontal measurement grids. Stack these into a 3D dataset, recording coordinates and whether height means national elevation or height above local ground. **Path gain g** is a linear received/transmitted power ratio: `received dBm = transmitted dBm + 10*log10(g)` (dBm measures power relative to 1 milliwatt), under the map's antenna assumptions. Account for receiver noise to obtain SNR (signal-to-noise ratio); include competing transmissions for SINR (signal-to-interference-plus-noise ratio). Standard maps assume an isotropic receiver (equal sensitivity in all directions); an oriented drone antenna may need separate path/receiver evaluation. [Map definitions](https://nvlabs.github.io/sionna/rt/api/radio_maps.html).
5. **Align and check before integration:** wind and RF may use different resolutions. Resample to explicit positions, mask solid/outside locations, and verify axes and units. The original wind helper stores `(x,y,z)` arrays; raw RF stacks use `(transmitter,z,y,x)`. Cell-centre versus grid-point positions also need checking. A 1 m RF map is a spatial sampling choice, not a guarantee of 1 m positioning accuracy. Validate against a known case or measurements before claiming site accuracy.

The professor's `openfoam_cost_map.py` and `sionna_cost_map.py` suggest a future connection. The first converts wind into planning penalties; its `sqrt(2k/3)` estimate assumes isotropic turbulence (equal fluctuation strength in every direction), not measured gusts or a proven drone energy model. The second uses radio strength and transmitter geometry to estimate a **CRLB** (Cramér–Rao lower bound: an ideal lower bound on unbiased positioning error under assumed measurement noise), not achieved localization accuracy. It assumes time-of-arrival measurements (estimating distance from signal travel time). Multipath (several reflected arrivals) can bias those distances, so a strong signal alone does not guarantee good positioning. This is not a completed positioning system. Sionna's installed version is unconfirmed; check the original helper against the selected API before use. Current online Sionna documentation is 2.2.0, not a version tested here.

Later, a planner could query wind, RF and obstacles at the same position to compare exposure, connectivity and clearance. Combining such costs requires explicit scales and weights; it does not couple radio physics into the airflow solver. Begin with static maps; moving-blade radio effects and live weather would be separate extensions. None of this requires additional QGIS crop views.

## Conclusions and next actions

The master scene is a reusable geographic starting point. The next research milestone is one checked single-turbine wind case and a trusted comparison. Extra QGIS views would not establish physical accuracy or improve the correction study. Detailed printable blades are not required by the actuator-disk method; correct dimensions, operating inputs and a trusted comparison matter first.

The original professor STEP/STLs and `windCase.zip` remain unchanged as an alternative example. Munich does not depend on them, but ask whether an original example is a required benchmark before disregarding it. Check its units and coordinate frame before reuse. The original cost-map scripts are later conversion tools; they are not completed Project 3 results.

Next, ask the professor:

1. Which trusted reference and accuracy measures should anchor Project 3? Should the first validation use a known simple turbine case before Munich?
2. Which wind speeds/directions, yaw states, atmospheric inputs and turbine coefficients are available or acceptable?
3. Should the first case use a simple benchmark or E-66, and is a two-turbine result ultimately required? Must the supplied geometry be used as a benchmark?
4. Which cheaper baseline should be tested first, and what accuracy/runtime improvement would make a correction useful?

Then rebuild one explicit, repeatable CFD case, fix numerical quality issues, validate it, and proceed through the outline above. Keep the uncorrected baseline and report unsuccessful cases as well as improvements.

## Source basis

The professor's current `README.md` takes priority. Project 3 from the older `Topic Proposals.pdf` (previously read, pages 3–4) supplies the validation/correction/generalization tasks transcribed here; the PDF is no longer needed alongside these notes. Site dimensions and source fingerprints are in `scene/qgis/scene_manifest.json`; retained pilot settings are in `wind-review/assumptions.json`.

Useful references: [OpenFOAM surface/volume meshing](https://www.openfoam.com/documentation/user-guide/4-mesh-generation-and-conversion/4.4-mesh-generation-with-the-snappyhexmesh-utility), [v2306 actuator-disk source](https://api.openfoam.com/2306/classFoam_1_1fv_1_1actuationDiskSource.html), [ParaView tutorials](https://docs.paraview.org/en/latest/Tutorials/index.html). The separate personal run guide explains the actual installed programs and commands; it is not part of this repository.
