# Project Structure

Reorganised on 2026-05-06. Every entry-point script lives at the root; everything
else is grouped by role under a code subfolder. Run results and logs live under
`Results/`. Anything in `Unused/` is intentionally off the active MATLAB path.

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
├── Results/
│   ├── 6Cams/                 Per-run .mat + .txt for 6-camera runs (146 files)
│   ├── 7Cams/                 Per-run .mat + .txt for 7-camera runs (224 files)
│   ├── 8Cams/                 Per-run .mat + .txt for 8-camera runs (144 files)
│   ├── Logs/                  GGA_RunsLog.mat, GA_RunsLog.mat, BatchLog_*.mat
│   └── 7C_UAV_setup.txt       Curated camera setup sheet (kept here for visibility)
│
├── figures/                   Output PDFs (from plotGARuns) + PNG plots
│
└── Unused/                    OFF the active MATLAB path
    ├── UniformCrossover.m     Alternative crossover not wired into RunGA
    ├── CameraConfigPlot.m     One-off plot script with hard-coded positions
    ├── computeOcclusionAngle.m  Duplicate of calculatePointOcclusion.m
    └── Autosaves/             *.asv MATLAB autosaves
```

## Where things live (quick lookup)

| Need to…                              | Path                                                   |
|---------------------------------------|--------------------------------------------------------|
| Run a single optimisation             | `runCameraOptimiser.m`                                 |
| Run the full batch sweep              | `batchRunGA.m`                                         |
| Generate every thesis figure          | `plotGARuns.m`                                         |
| Inspect / tweak the best result       | `analyseConfiguration.m` or `ConfigAnalyser.m`         |
| List runs in the master log           | `viewGALog.m`                                          |
| Find a per-run `.mat` / `.txt`        | `Results/<N>Cams/<N>Cams_Run_<timestamp>.<ext>`        |
| Load the master log                   | `Results/Logs/GGA_RunsLog.mat`                         |
| Find batch-sweep state                | `Results/Logs/BatchLog_<timestamp>.mat`                |
| Add a new cost function               | drop in `CostFunctions/`                               |
| Add a new plotting helper             | drop in `Plotting/`                                    |

## How file resolution works

Every entry-point script and every interactive function calls
`addProjectPaths()` at the top. That walks the project root and adds
`GA_Core/`, `CostFunctions/`, `Geometry/`, `Setup/`, `Plotting/` and
`Analysis/` to the MATLAB path, but **not** `Unused/`. So any function in
those subfolders is reachable by name from anywhere.

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
* **`Unused/Autosaves/*.asv`** — MATLAB editor backups. Kept (rather than
  deleted) so they're available if you ever need to recover an in-progress
  edit, but out of the way.
