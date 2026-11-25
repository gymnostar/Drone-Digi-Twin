# DLR campus wind atlas

Wind at any point of the DLR Oberpfaffenhofen campus, for **any wind direction and speed**, from 12 steady CFD runs. This is the wind input for the cost function and the planner (`DLR/planner/`).

```python
import sys; sys.path.insert(0, "DLR/atlas")
from wind_atlas import WindAtlas
atlas = WindAtlas()                          # loads grid.noindex/*.npz (~0.5 s)
U = atlas.field(theta=250, speed=6.0)        # (169, 161, 49, 3) m/s on the shared grid, NaN in solids
u = atlas.query([[0, 0, 60]], 250, 6.0)      # (N, 3) at arbitrary points, trilinear
```

`theta` is the meteorological **from** direction (deg clockwise from north; 270 = west wind). `speed` is the wind at 10 m above open ground upstream, which is what the EDMO airfield reports. Coordinates are the case frame: x east, y north, z up (m). Case = scene_local − (−335, 279), and the scene origin is UTM32 E670000 N5328000.

## How it works (three shortcuts, each one tested)

| Step | What | Why it is allowed | Test |
|---|---|---|---|
| 1. One mesh | All runs use the mesh-study base mesh: 16 m background, 4 m over campus, 2 m at walls; 1.19 M cells. Only the inlet/outlet split of the cylinder and the inflow direction change. | Mesh study (`DLR/meshStudy`), coarse vs. 1 m mesh: ≤ 2 % at most probes 20–40 m above ground, 9 % in a street canyon (not mesh-converged there: GCI 15 %), and up to 15 % at 5 m in a building wake. | `DLR/meshStudy/report` |
| 2. One speed | Every direction is solved at 5 m/s, then rescaled: U ∝ speed, k ∝ speed². | Reynolds independence: around sharp-edged buildings, with rough walls (z0), the flow pattern doesn't depend on speed. | One extra 10 m/s run at 240° |
| 3. Direction blending | Between two solved directions (30° apart): rotate each field's horizontal vectors to the requested direction, then blend linearly. | The far field is then exact (blending raw vectors shrinks it by cos(gap/2)). Near buildings it can smear wakes, so it's the weakest step. | Leave-one-out (60° gaps) + one extra 255° run (the real 30° gap, in the dominant sector) |

Each run is resampled **once** onto a regular grid: x −400..440 and y −360..440 every 5 m, z 10..130 m (absolute) every 2.5 m. Points inside buildings or below the terrain are NaN, and `valid` is the air mask. A regular grid gives every consumer (planner, RF side, AirSim) an O(1) lookup in the same frame. The CFD mesh is unstructured, so it can't do that.

## Runs

The order follows the EDMO 2020–2025 wind climate (share of reports with a defined direction / share of wind energy ΣU³; calm and variable reports, about 14 %, excluded), so the most useful directions finish first:

| from | 240 | 270 | 210 | 60 | 180 | 30 | 150 | 90 | 300 | 0 | 330 | 120 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| directional reports % | 22.9 | 10.3 | 15.8 | 8.3 | 8.2 | 7.1 | 6.8 | 6.1 | 4.1 | 3.6 | 3.7 | 3.1 |
| wind energy % | 57.7 | 15.1 | 7.3 | 5.1 | 1.5 | 3.2 | 1.2 | 2.6 | 2.8 | 1.4 | 1.7 | 0.4 |

- 270° is the mesh study's coarse run (1000 iterations).
- **210° oscillates** (a limit cycle of steady RANS around the tall blocks, amplitude ≈ 1 % of Uref); the stored field is one snapshot. The other runs converge monotonically (probe drift over the last 50 iterations < 0.1 % of Uref at drone heights; those probes sit near the weather station and were laid out for west wind).
- 240, 210, 255 and the 10 m/s run use 450 SIMPLE iterations; the remaining directions use 350 (cut overnight when the machine was busy: the 240 run's drift at drone heights was 0.09 % of Uref per 50 iterations). The mesh study showed every probe within 0.5 % of its final value by iteration ~360, and `validate_atlas.py` re-checks the probe drift of every run.
- Validation-only runs: 240° at 10 m/s (Reynolds check) and 255° at 5 m/s (blending check).

## Files

| File | Role |
|---|---|
| `setup_atlas.py` | Builds `cases.noindex/base` (mesh study mesh, cylinder side merged back into one patch) and one case per run (`dTTT_uSS`). |
| `run_atlas.sh` | Overnight runner: topoSet/createPatch for the direction → simpleFoam → resample → `grid.noindex/<run>.npz`. Skips runs that have a `DONE` file. |
| `sample_grid.py` | pvpython: resamples one solved case onto the shared grid. |
| `wind_atlas.py` | The lookup: loading, speed scaling, direction blending, trilinear queries, ground height. |
| `validate_atlas.py` | Convergence, Reynolds and blending checks; lookup timing. Writes `results/validation.json`. |
| `make_figures.py` | `results/atlas_maps_*.png`, `atlas_reynolds.png`, `atlas_loo.png`. |

`cases.noindex/` (~4 GB) and `grid.noindex/` (~25 MB per run) are git-ignored; rebuild them with `nohup caffeinate -i sh DLR/atlas/run_atlas.sh > DLR/atlas/run_atlas.out 2>&1 &`. The `.noindex` suffix keeps Spotlight from indexing solver output.

The 270° grid comes from the mesh study and has its own command:

```sh
printf "theta_from_deg 270\nUref_10m 5\nflowDir 1 0 0\n" > DLR/meshStudy/coarse/RUN_INFO
/Applications/ParaView-6.2.0.app/Contents/bin/pvpython DLR/atlas/sample_grid.py DLR/meshStudy/coarse DLR/atlas/grid.noindex/d270_u05.npz
```

## Limits

- **Heights above ground are +0 to 2.5 m:** `ground()` returns the first air node of the 2.5 m grid, and `terrain()` smooths it with a 60 m minimum filter.
- **Extrapolation is refused:** a direction lying in a gap of more than 60° between solved runs raises an error instead of blending across half the compass.


- **Model, not reality.** Steady RANS k-ε with a neutral log-law inflow (z0 = 0.1 m). No gusts, no thermal stability, no validation against measurements yet (the DLR weather-station data is pending IT clearance).
- **Terrain beyond the original strip is invented flat ground** (see `DLR/README.md`).
- **Mean wind only.** The turbulence k is stored (`atlas.field(..., var="k")`) but the planner doesn't use it yet.
- **Low altitude near buildings is the least accurate.** The mesh study shows up to 13 % differences at 5 m in building wakes.
