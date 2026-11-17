# Buildings and turbine wind study

**Question:** can we generate wind maps more cheaply while keeping their numerical error small? We tested the full Niigata building district and a SWiFT V27 turbine separately, then solved their interacting flow in one domain.

| Study | Result | Interpretation |
|---|---|---|
| Isolated turbine | About **6× faster** complete map generation; wake vector RMSE **1.38–1.63%** on three reserved numerical conditions | Successful numerical cost reduction within the tested setup |
| Buildings | **112,358 cells, 81.85 s**; city RMSE **6.42%**, above the 5% target | Cheap solve works, but accuracy target is unmet and stricter references did not converge |
| Buildings + turbine | **422,894 cells**; baseline passes its numerical checks; stricter reference timed out | Actual interacting flow, but no certified joint speedup; candidate mesh ready, solve pending |

Errors compare calculations with a numerical reference, **not with real wind**. These are research results, not flight-safety or physical-validation claims.

### Show the actual simulation in 3D

Open [Open3D.command](../review-export/Open3D.command) with ParaView installed. It opens the coupled baseline. To switch views from the repository root:

```sh
./review-export/Open3D.command joint-baseline
./review-export/Open3D.command turbine-cheap
./review-export/Open3D.command turbine-reference
```

The [viewer script](../review-export/view.py) opens the selected view without running a solver. Keep the exported data and state files together. The generated states use the actual saved CFD volume, not reconstructed plots.

In ParaView, drag to rotate, scroll to zoom, and use the eye icons in **Pipeline Browser** to show or hide buildings, slices and flow paths. Select **Velocity at rotor height**, change **SliceType → Origin → Z**, then click **Apply** to move the horizontal section. Select **Vertical wake section** for the wake through its centre. The source named **Saved CFD volume** retains the full computed domain; **Display region** only clips the view.

Rotate the buildings, follow the flow downstream, then inspect the rotor wake. Change the slice height to show how shelter changes with height. Colours show predicted mean speed in m/s. Streamlines follow the computed mean velocity: they are not gust trajectories or drone routes. The rotor marker locates the modelled force; it is not moving blade geometry. The coupled view is the baseline, whose stricter reference remains provisional.

<details>
<summary>Explanation for presenting: methods, results, limitations and next steps</summary>

### What we did

The pipeline is **geometry → mesh → solve → wind map → comparison**. A mesh divides the air into cells; the solver balances mass and momentum in them. We measured preparation, meshing, solving and extraction together. Reading an existing map is a different operation and is not the reported speedup.

Before testing, we fixed provisional targets: at least **2× faster**, vector RMSE ≤**5%** of incoming speed, P95 ≤**15%**, and ≥**95%** coverage. RMSE (root mean square error) measures typical differences and penalises large ones. Vector error includes direction and speed. P95 is the error below which 95% of sampled errors fall, not a maximum or confidence interval. Coverage counts common valid samples against reference-valid points. These are proof-of-concept targets, not agreed drone tolerances.

Runs use three CPUs, 5 GiB, ten-minute mesh/solve limits and bounded repairs. The [frozen protocol](../scene-selection/poc-protocol.json) records conditions, thresholds and budgets; it has not been edited retrospectively. Failed trials remain in the structured evidence.

### Main implementation choices

- **OpenFOAM v2306, steady RANS:** the solver predicts average turbulent flow. RANS models unresolved fluctuations rather than resolving individual gusts. SIMPLE/SIMPLEC are iteration methods coupling pressure and velocity.
- **Full building geometry:** Niigata after-construction geometry is scaled from its 1:250 experiment into full-size metres, without cropping. The numerical scenario has its own declared inflow; it does not recreate the tunnel experiment.
- **Neutral atmospheric boundary layer:** wind increases with height using a coherent inlet/ground treatment. The common development condition is 8.7 m/s at 32.1 m, roughness length 0.01158 m (a ground-drag parameter). Empty-domain runs check profile transport before adding obstacles.
- **One global realizableKE closure:** a closure supplies unresolved turbulent mixing. This stock model was selected from turbine development evidence and used across the coupled domain. Its physical superiority around buildings is unproven.
- **Actuator disk:** a 27 m rotor at 32.1 m hub height slows the air through a momentum sink (opposing force), without moving blades. Locked v3 rotor curves supply coefficients. Coupled forcing uses solved local inflow. Tower/nacelle drag, swirl and tilt are omitted.
- **Fixed force support and queries:** the same physical cells must receive rotor force across mesh comparisons; equal total thrust is insufficient. Building vicinity, rotor wake and proposed flight region are evaluated separately. Missing samples are not moved or filled.

### Turbine: successful bounded comparison

The cheap recipe reduces wake refinement while retaining fixed fine rotor-force cells and a common near-ground band. The first stopping rule failed settling (the solution was still changing); tighter stopping was frozen before testing whole reserved conditions, with no later tuning.

| Reserved condition: speed / roughness | Reference / cheap map time | Speedup | Wake RMSE / incoming speed |
|---|---|---|---|
| 8.0 m/s / 0.008 m | 130.82 / 21.47 s | 6.09× | 1.63% |
| 9.5 m/s / 0.015 m | 135.72 / 22.73 s | 5.97× | 1.38% |
| 8.7 m/s / 0.025 m | 140.82 / 23.38 s | 6.02× | 1.39% |

Each comparison has 19,800/19,800 common wake samples and passes convergence and spatial-force checks. P95 is 3.26–3.90%. An enlarged-domain check changes wake RMSE by 0.506% of incoming speed. This is limited boundary-sensitivity evidence, not proof of mesh independence. Authoritative reports in [turbine results](turbine/results/) end in `verified-comparison.json`.

### Buildings: accuracy target not established

Three strict-reference attempts missed convergence limits. Residuals (remaining numerical equation imbalance) plateaued or oscillated even when overall flow balance looked good. Two cheap meshes also failed geometry checks. The final cheap mesh/solve passes its own checks, but RMSE against the provisional reference is 6.42% across the city, 7.24% in the lower city and 7.40% downstream. Whole-city coverage is 30,301/30,505 points.

No certified city speedup is reported. Unused direction holdouts remain untouched because development did not establish a valid reference/method. See the [city summary](buildings/results/poc-study-summary.json).

### Coupled flow: actual interaction, unfinished optimization

The hypothetical rotor is downstream at (0, −280, 32.1) m, with a local-inflow monitor at (0, −215, 32.1) m and rotor-bottom clearance of 18.6 m. One solve computes building/turbine interaction; separate wind fields are never added.

Different mesh recipes slightly moved rotor-source vertices during snapping (fitting cells to surfaces). We rejected their unequal force support and used the identical reference mesh recipe for baseline/candidate. **Mesh-cost reduction failed to transfer; the remaining comparison tests solver/stopping settings only.** The full joint mesh exceeds the initial cheap-cell target; its declared cap adjustment is recorded.

The baseline converged in 364 iterations. Its recorded pipeline cost is 306.76 s, a lower bound because original host force-support verification overhead was not timed. The strict reference timed out at 600.20 s; its saved iteration-900 field is provisional. Baseline/reference differences are small (rotor-wake RMSE 0.512%, full common coverage), but cannot certify reference accuracy or speedup. The candidate mesh has matching force support; its solve is pending.

An invalid early trial caused by a dictionary-substitution bug was stopped and excluded; generated solver settings are now checked before execution. Failed implementation/support trials and costs remain in the [joint summary](joint/results/poc-study-summary.json).

### Physical evidence and limits

The original Niigata tunnel pilot has speed-ratio RMSE about 0.198 on 77/80 samples, with unresolved inlet/turbulence metadata and profile drift. The new full-scale neutral-inflow toy is a different scenario.

The initial SWiFT standard-k–epsilon comparison has measured wake RMSE about 1.11 m/s. Its observed periods already entered development, so historical date splits are exploratory. The three reserved conditions above test numerical generalisation, not independent physical observations. There are no joint interaction observations.

A numerical reference is a more demanding calculation, not physical truth. We have not demonstrated instantaneous gust accuracy, flight safety or a validated integrated wind/RF/drone system.

### Next steps and pace

1. Finish the one remaining capped candidate solve. This alone will not fix the failed reference.
2. Establish a checked city/joint reference through a focused numerical study: revised steady setup, justified longer budget, or time-dependent flow. Do not loosen thresholds just to pass.
3. Obtain better matched inflow metadata or independent observations. Agree which wind quantities and errors matter for drone planning.
4. Expand turbine conditions first: comparable cheap cases took 21–23 s versus 131–141 s for references. This is measured throughput within this setup, not guaranteed scaling for new geometry. City/joint expansion is gated by the reference problem.
5. Then expand terrain, turbine count or fidelity, or build a surrogate (a smaller model approximating simulations). Include database/training cost. More compute increases throughput but does not fix modelling gaps. Wind/RF exchange and drone planning follow; reinforcement learning is optional.

The project has a reproducible pipeline and a successful isolated-turbine result. City reference convergence and physical validation remain the main blockers. No reliable completion date for those repairs follows from these measurements.

</details>

### Files and reproduction

- [scenes/buildings](buildings/): input lock, case builder, runner, diagnostics and metrics.
- [scenes/turbine](turbine/): benchmark input locks, actuator-disk preparation, evaluation and comparisons.
- [scenes/joint/study.py](joint/study.py): coupled preparation and checks, reusing separate-scene helpers.
- Downloaded inputs and generated `runs/` are ignored; structured metrics and source hashes retain scientific evidence. Exported 3D data are generated deliverables.

Scripts provide `--help`. For example, `python3 scenes/joint/study.py --help` lists build, mesh, solve, sample and evaluate stages. Existing cases are preserved: use a fresh case name for a new experiment. Viewing saved results does not run CFD.
