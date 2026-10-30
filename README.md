# Drone Digi Twin

## Project Overview

Develop a simulation environment for training and validating autonomous drones in complex environments. The digital twin combines terrain, buildings, wind turbine infrastructure, CFD wind simulations, and RF propagation.

**Initial scope:** One wind turbine, surrounding buildings, terrain, and simplified wind modelling. Expand the environment progressively to limit computational costs.

**Later:** Implement RL workflow
## 1. Environment and Geometry

### Terrain and Buildings

Geographical data sources:
- [Bavarian LoD2 3D Building Models](https://geodaten.bayern.de/opengeodata/OpenDataDetail.html?pn=lod2)
- DGM1 digital terrain model (1 m resolution).
- DLR building geometry.

### Wind Turbine

- [Wind Turbine Geometry – Energy Encyclopedia](https://www.energyencyclopedia.com/en/free-downloads/3d-printing/wind-turbine-14)

Start with one turbine. Multiple turbines can be introduced later to investigate wake interactions, downstream velocity reduction, and increased turbulence.

### Geometry Preparation

Use **Autodesk Fusion 360** to assemble terrain, buildings, and turbine geometry. Prepare separate representations for CFD and RF simulations, as their geometric requirements differ.

## 2. CFD – OpenFOAM

### Turbine Modelling

Initially use OpenFOAM's `actuationDiskSource` to represent the rotor as a momentum sink without explicitly modelling blades.

Future extensions:
- Actuator line modelling.
- Explicit blade geometry.
- Large Eddy Simulation (LES).

### Mesh Generation

Use OpenFOAM's `blockMesh` and `snappyHexMesh` workflow:

1. **Background mesh:** Generate the computational domain using `blockMesh`.
2. **Castellation:** Refine cells around STL geometry.
3. **Snapping:** Adjust mesh points to match geometry surfaces.
4. **Boundary layers:** Optional; not initially required for the large-scale environmental model.

Reference: [Wind Flow Predictions Around Offshore Wind Turbines](https://www.emerald.com/ec/article/43/7/2807/1298297/Development-of-a-fast-numerical-tool-for-wind-flow)

The study uses poly-hexcore meshing as an alternative approach.

## 3. RF Propagation – Sionna

Sionna requires triangulated surface geometry and appropriate electromagnetic material properties.

Unlike CFD, RF ray tracing is sensitive to:
- Surface materials and conductivity.
- Reflections and scattering.
- Sharp edges and geometric details.
- Metal structures and window frames.

Therefore, simplify geometry for CFD while retaining relevant details for RF propagation.

**Initial target resolutions (to be validated):**
- CFD: approximately 2 m.
- RF (Sionna): approximately 1 m.

## 4. Development Workflow

| Stage | Development | Software |
|---|---|---|
| 1 | Terrain and building geometry | Fusion 360 |
| 2 | Wind turbine integration | Fusion 360 |
| 3 | CFD mesh and wind simulation | OpenFOAM |
| 4 | RF geometry and propagation maps | Sionna |
| 5 | Wind and RF map integration | Python |
| 6 | Autonomous mission planning | Python |
| 7 | Drone simulation and validation | Unreal Engine + Project AirSim |

## 5. Software and Versions

| Software | Version / Configuration | Purpose |
|---|---|---|
| OpenFOAM | v2306 (DLR HPC, to be confirmed) | CFD and wind modelling |
| Fusion 360 | N/A | Geometry preparation |
| Sionna | To be confirmed | RF ray tracing |
| Python | 3.9.13 (local Windows environment) | Simulation integration and mission planning |
| Project AirSim Python Client | 1.0.2 (`projectairsim==1.0.2`) | Python API for drone control and simulation |
| NumPy | 2.0.2 | Numerical computing |
| OpenCV | 5.0.0.93 | Image processing |
| Matplotlib | 3.9.4 | Data visualization |
| Shapely | 2.0.7 | Geometric operations |
| pynng | 0.9.0 | Communication between Python client and simulator |
| msgpack | 1.1.2 | Data serialization |
| Azure Core | 1.39.0 | Azure service integration |
| Azure Storage Blob | 12.30.3 | Cloud storage integration |
| Unreal Engine | 5.7 | 3D simulation environment |
| Project AirSim | 1.0.1 (local installation) | Autonomous drone simulation |
| ParaView | To be confirmed | CFD post-processing |

*Note: Project AirSim simulator (1.0.1) and Python client (1.0.2) have separate version numbers.*
**Long-term objective:** Integrate CFD wind fields, RF propagation maps, and autonomous mission planning into a digital twin for energy-aware drone navigation in complex environments.
