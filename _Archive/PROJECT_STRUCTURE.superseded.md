# Project Structure

Reorganised on 2026-05-06, flattened on 2026-08-05. Every entry-point script
lives at the root; everything else is grouped by role under a code subfolder.
Anything in `_Archive/` is intentionally off the active MATLAB path.

## How to start

From the workspace root (`MSc/MATLAB/`), run:

```matlab
START_HERE
```

That is the whole setup. It puts the RVC3 toolbox on the path (so
`CentralCamera` and `se3` resolve), calls `addProjectPaths`, changes into this
folder and lists the entry points. Opening `rvc3setup.prj` by hand is no longer
needed, and neither is navigating down through subfolders — this project used to
sit five levels deep under `Simulation/Optimising Camera Placement/Genetic
Algorithm/Real-coded GGA/Optimal Camera Placement Genetic Algorithm/` and is now
the root of its own repository.

## ⚠️ `Results/` is empty on this machine

The current (end-of-July 2026) GA runs live on the **simulation machine** and
have never been committed — `Results/` is in `.gitignore`, so they were never
picked up. Until they are force-added there and pushed, this machine has no
valid run data and `viewGALog` / `analyseConfiguration` will not find anything.

See the "Getting the results off the simulation machine" section of the
workspace `README.md` for the exact sequence.

`_Archive/Results_pre-bugfix/` holds 514 files from Mar–May 2026. Those predate
the FOV and normalisation fixes and are **not valid results** — neither
`viewGALog` nor `resolveRunPath` will fall back to them, deliberately.

## Layout

```
<root>/
├── runCameraOptimiser.m       Main entry script — single optimisation run
├── batchRunGA.m               Batch sweep across cam counts / cost functions / etc.
├── app1-Shannons-PC.m         Older standalone entry script
├── plotGARuns.m               Generates every figure for the thesis from logs
├── analyseConfiguration.m     Interactive analysis of the best CF3 result
├── ConfigAnalyser.m           Class form of the same analysis tool
├── viewGALog.m                Print a tabular summary of the master log
├── addProjectPaths.m          Bootstrap — adds every code subfolder to the path
├── PROJECT_STRUCTURE.md       This file
│
├── GA_Core/                   The genetic algorithm engine
│   ├── RunGA.m
│   ├── Tournament.m
│   ├── RouletteWheelSelection.m
│   ├── DoublePointCrossover.m
│   ├── Mutate.m
│   ├── SortPopulation.m
│   ├── initialPopulation.m
│   └── fixPoorCameras.m
│
├── CostFunctions/             Cost evaluation + visibility / triangulability
│   ├── combinedCostFunction.m
│   ├── resUncertainty.m
│   ├── resUncertaintyCost.m
│   ├── dynamicOcclusion.m
│   ├── dynamicOcclusionCost.m
│   ├── computePointUncertainty.m
│   ├── calculateOccludedSections.m
│   ├── calculatePointOcclusion.m
│   ├── checkSectionTriangulability.m
│   ├── checkTriangulability.m
│   ├── findFrontCameras.m
│   └── findVisibleCameras.m
│
├── Geometry/                  Pyramids, vertices, target space, ellipsoid, camera setup
│   ├── buildPyramidSurf.m
│   ├── calcVertices.m
│   ├── isInsidePlanes.m
│   ├── quantToWorld.m
│   ├── optimiseEllipsoid.m
│   ├── generateSectionCentres.m
│   ├── generateTargetSpace.m
│   └── setupCameras.m
│
├── Setup/                     Parameter / problem setup
│   ├── setupCostParams.m
│   ├── setupGAparams.m
│   ├── setupHardwareSpecs.m
│   └── setupProblem.m
│
├── Plotting/                  Plot helpers & styling
│   ├── plotResults.m
│   ├── plotCoverageHeatmap.m
│   ├── plotGA_ComputationTime.m
│   ├── plotGA_Convergence.m
│   ├── plotGA_CostBoxPlots.m
│   ├── plotGA_FactorEffects.m
│   ├── plotGA_PopulationDiversity.m
│   ├── plotGA_WarmColdEffect.m
│   ├── visualizeCameraCoverage.m
│   ├── gaPlotStyle.m
│   ├── gaStatsHelpers.m
│   ├── applyThesisStyle.m
│   └── getModalityLabel.m
│
├── Analysis/                  Save / load / resolve helpers
│   ├── saveResults.m          Writes per-run mat+txt into Results/<N>Cams/
│   ├── loadGARuns.m           Loads + filters the master log
│   └── resolveRunPath.m       Maps a bare RunFilename to its full path
│                              under Results/<N>Cams/ (back-compat with
│                              existing log entries that store basenames)
│
├── Results/                   EMPTY ON THIS MACHINE — see warning above
│   ├── 6Cams/  7Cams/  8Cams/  Logs/    (awaiting push from the sim machine)
│   ├── Sensitivity/           spacing sweep .mat files (UAV + UGV)
│   └── Sweep_PopVsGen/        population-vs-generations sweeps (Jul 2026)
│
├── figures/                   Output PDFs (from plotGARuns) + PNG plots
│
├── Reference/
│   └── Rahimian-P2I-python/   Rahimian et al. Python implementation (reference)
│
└── _Archive/                  OFF the active MATLAB path — nothing deleted
    ├── Unused/                UniformCrossover.m, CameraConfigPlot.m,
    │                          computeOcclusionAngle.m, old autosaves
    ├── pre-bugfix-snapshots/  Dated snapshots (pre_FOVfix, pre_normfix)
    ├── figures_pre-bugfix/    Figures generated before the fixes
    ├── Results_UGV_fine_grid_archive/
    ├── EarlyRootScripts/      optiTrackConfig.m, sectionCentres.m (superseded)
    ├── autosaves/             *.asv MATLAB editor backups
    └── build_cache/           slprj Simulink build artefacts (regenerable)
```

## Where things live (quick lookup)

| Need to…                              | Path                                                   |
|---------------------------------------|--------------------------------------------------------|
| Run a single optimisation             | `runCameraOptimiser.m`                                 |
| Run the full batch sweep              | `batchRunGA.m`                                         |
| Generate every thesis figure          | `plotGARuns.m`                                         |
| Inspect / tweak the best result       | `analyseConfiguration.m` or `ConfigAnalyser.m`         |
| List runs in the master log           | `viewGALog.m`                                          |
| Set the path and get started          | `START_HERE` from the workspace root                   |
| Find a per-run `.mat` / `.txt`        | `Results/<N>Cams/<N>Cams_Run_<timestamp>.<ext>`        |
| Load the master log                   | `Results/Logs/GGA_RunsLog.mat`                         |
| Find batch-sweep state                | `Results/Logs/BatchLog_<timestamp>.mat`                |
| Compare GA vs the ad-hoc OptiTrack rig| `saveOptiTrackAsRun.m` (stay in GA mode — see below)   |
| Add a new cost function               | drop in `CostFunctions/`                               |
| Add a new plotting helper             | drop in `Plotting/`                                    |

## How file resolution works

Every entry-point script and every interactive function calls
`addProjectPaths()` at the top. That walks the project root and adds
`GA_Core/`, `CostFunctions/`, `Geometry/`, `Setup/`, `Plotting/` and
`Analysis/` to the MATLAB path, but **not** `_Archive/`. So any function in
those subfolders is reachable by name from anywhere.

Note that `START_HERE` deliberately does **not** add
`1_CameraPlacement/OptiTrackConfig/` at the same time as this project. Those two
folders define 15 same-named functions (`buildPyramidSurf`, `calcVertices`,
`resUncertainty`, `combinedCostFunction`, …) with different bodies, so having
both on the path would silently resolve to the wrong one. Use
`START_HERE optitrack` to switch.

Per-run `.mat` files are referenced from log entries by their **basename
only** (e.g. `7Cams_Run_20260331_001929.mat`). The helper
`resolveRunPath(filename, numCams)` in `Analysis/` maps that basename to
the full path under `Results/<N>Cams/`. This means existing log entries —
written before the move — continue to resolve without rewriting any data.

`saveResults.m` (and the inline save block in `app1-Shannons-PC.m`) write
new runs straight into `Results/<N>Cams/` and append to
`Results/Logs/GGA_RunsLog.mat`, keeping the layout self-consistent going
forward.

## Why some things were marked unused

* **`UniformCrossover.m`** — defined but never called from `RunGA` or
  anywhere else. The active crossover operator is `DoublePointCrossover`.
* **`CameraConfigPlot.m`** — top-level script with hard-coded camera
  positions, never invoked by any other file. Looks like a one-off
  visualisation snippet.
* **`computeOcclusionAngle.m`** — the file declares a function called
  `calculatePointOcclusion`, which is a duplicate of the standalone
  `CostFunctions/calculatePointOcclusion.m`. Calling
  `computeOcclusionAngle(...)` in MATLAB would actually fail (function
  name doesn't match the file name). Moved out so it can't be picked up.
* **`_Archive/autosaves/*.asv`** — MATLAB editor backups. Kept (rather than
  deleted) so they're available if you ever need to recover an in-progress
  edit, but out of the way.
