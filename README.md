# Optimal Camera Placement — Guided Genetic Algorithm

A guided genetic algorithm (GGA) that finds where to mount the cameras of an
optical motion-capture system so that the tracked volume is covered as well as
possible.

Each candidate solution is a **chromosome** describing every camera's pose —
6 real-valued genes per camera `[Xc, Yc, Zc, alpha, beta, gamma]` (position in
metres, orientation in radians). For 7 cameras that is a 42-gene chromosome.
The GA evolves a population of these and returns the lowest-cost arrangement.

---

## 1. Quick start

**From the workspace root (`MSc/MATLAB/`), in MATLAB:**

```matlab
START_HERE
```

That is the only setup step. It puts Peter Corke's RVC3 toolbox on the path
(this project needs `CentralCamera` and `se3` from it), adds every code folder
here, and moves you into this project.

Then, to reproduce a single optimisation run:

```matlab
edit runCameraOptimiser    % set your inputs at the top (see §3)
runCameraOptimiser         % run it
```

To look at results that already exist:

```matlab
viewGALog                  % table of every logged run
analyseConfiguration       % inspect the best 7-camera combined-cost result
plotGARuns                 % regenerate every thesis figure into figures/
```

You do not need to know where any of these files live — `START_HERE` puts them
all on the MATLAB path, so calling them by name works from anywhere.

> **Run `START_HERE` once per MATLAB session**, not once per script. The path
> persists for the whole session.

---

## 2. Requirements

- MATLAB (developed against R2026a) with the Image Processing and Computer
  Vision toolboxes.
- **RVC3-MATLAB** (Peter Corke's *Robotics, Vision & Control* toolbox), at
  `../../3_External/RVC3-MATLAB`. `START_HERE` adds it automatically; you do
  **not** need to open `rvc3setup.prj`.

---

## 3. Running the GA

### Single run — `runCameraOptimiser`

This is a **script**, not a function: you edit the values at the top and run it.
The settings that matter, all in the `%% User Inputs` block:

| Variable | Meaning | Typical |
|---|---|---|
| `numCams` | number of cameras to place | `7` |
| `volume` | tracked volume `[xmin xmax; ymin ymax; zmin zmax]` in metres | `[-4 4; -4 4; 0 4]` |
| `cameraLowerBounds` / `cameraUpperBounds` | per-camera search bounds (mounting constraints), 6 values each | wall-mount limits |
| `targetType` | `1` = UAV (whole volume), `2` = UGV (floor slab) | `2` |
| `targetMode` | `1` = uniform grid, `2` = centre-concentrated grid | `1` |
| `spacing` | evaluation grid spacing in metres (x–y) | `1` |
| `UGV_maxHeight`, `UGV_zSpacing` | slab height and z-step, UGV only | `0.5`, `0.25` |
| `costFunctionType` | `1` = resolution uncertainty, `2` = dynamic occlusion, `3` = combined | `3` |
| `maxGenerations` | GA generations | `100` |
| `populationSize` | auto-scales as `numCams × 6 × 10` | `420` for 7 cams |
| `warmStartUsed` / `warmStartBestSol` | seed from a previous chromosome | `false` |

For `costFunctionType = 3` the two weights must sum to 1:

```matlab
specs.WeightUncertainty = 0.5;
specs.WeightOcclusion   = 0.5;
```

Finer `spacing` means many more evaluation points and a much slower run — it is
the main cost/runtime lever. Runs of a few hours are normal.

When it finishes it prints the best cost and elapsed time, draws the coverage
and convergence plots, and **saves automatically** (see §5).

> **Warm start:** only ever seed from a *post-bugfix* chromosome. Earlier runs
> were optimised against a corrupted cost surface and will mislead the GA.

### Batch sweep — `batchRunGA`

Sweeps combinations of conditions and repeats each one. All options are
name–value pairs; defaults in brackets:

```matlab
batchRunGA()                                  % full default sweep
batchRunGA('CameraRange', 7, 'NumRepeats', 3) % just 7 cameras, 3 repeats
batchRunGA('DryRun', true)                    % list what would run, run nothing
batchRunGA('ResumeFrom', 12)                  % resume an interrupted sweep
```

Key parameters: `CameraRange` `[6 7 8]`, `CostFunctions` `[1 2 3]`,
`TargetTypes` `[1 2]`, `GridModes` `[1 2]`, `Spacings` `[1.0]`,
`NumRepeats` `5`, `MaxGenerations` `100`, `SkipWarmStart` `false`.

**Start with `'DryRun', true`.** The full default sweep is hundreds of runs.

`batchRun7Cameras` is a narrower preset for the 7-camera case.

---

## 4. The three cost functions

| Type | Name | Measures |
|---|---|---|
| 1 | Resolution uncertainty | how precisely a point can be triangulated |
| 2 | Dynamic occlusion | how robust coverage is when cameras get blocked |
| 3 | Combined (CF3) | weighted sum of both — the one used for results |

CF3 is not a raw sum. `cf3Terms.m` is the single source of truth: each objective
is min–max scaled as `(raw − utopia)/norm` so utopia → 0 and nadir → 1, then
weighted and summed. The utopia/nadir constants come from
`ParameterTesting/getNormConstants.m`.

**Lower cost is better** for all three.

---

## 5. Where results go

Saving is automatic. `Analysis/saveResults.m` writes, per run:

```
Results/<N>Cams/<N>Cams_Run_<timestamp>.mat    full result + specs + chromosome
Results/<N>Cams/<N>Cams_Run_<timestamp>.txt    human-readable summary
Results/Logs/GGA_RunsLog.mat                   master log, one row per run
```

`Results/` is tracked in git so runs done on the simulation machine can be
pulled onto a laptop for analysis.

The master log stores run files by **basename only**;
`Analysis/resolveRunPath.m` maps a basename back to its full path, so log
entries keep working even if the folders move.

---

## 6. Reading your results

| Task | Command |
|---|---|
| List every logged run | `viewGALog` |
| Filter the log | `viewGALog('NumCameras', 7)` |
| Inspect / hand-tweak the best config | `analyseConfiguration` |
| Same thing, object form | `ConfigAnalyser` |
| Regenerate all thesis figures | `plotGARuns` |
| Per-component cost breakdown | `reportCostBreakdown` |
| Full 7-camera analysis suite | `analyseBest7CamConfigs` |

`analyseConfiguration` returns a handle you can poke at interactively:

```matlab
cfg = analyseConfiguration('NumCameras', 7);
cfg.printCameras()                              % list poses
cfg.moveCamera(3, [4.0 -2.0 2.5], [])           % move camera 3, keep orientation
cfg.roundAll(0.10, 5)                           % snap to 10 cm / 5 degrees
```

That last one matters in practice: it snaps an optimised layout to values a
person can actually mount, and re-costs it so you can see what the rounding cost
you.

`plotGARuns` writes PDFs into `figures/`.

### Comparing against the real (ad-hoc) OptiTrack rig

Use **`saveOptiTrackAsRun`** — and stay in GA mode.

It builds the OptiTrack rig's chromosome, evaluates it with *this project's*
CF3 cost function on the same grid, and writes it in `saveResults` format so it
loads alongside GA runs and is directly comparable.

Do **not** use the cost functions in `../OptiTrackConfig/`. They were last
touched 2025-12-09 and have diverged: `combinedCostFunction` there takes
`(cameras, specs)` and divides by `specs.PreComputed` norms, while this
project's takes `(cameraChromosome, specs)` and applies `cf3Terms` utopia
scaling. The two do not produce comparable numbers. Treat `OptiTrackConfig/` as
rig visualisation and calibration tooling only.

---

## 7. Folder layout

Only `addProjectPaths.m` and the docs sit at the root; everything else is
grouped by what it is for.

```
GeneticAlgorithm/
├── README.md               this file
├── addProjectPaths.m       defines the project root — must stay here
│
├── Run/                    ► scripts you execute to perform a run
│   ├── runCameraOptimiser.m    single run (main entry point)
│   ├── batchRunGA.m            batch sweep across conditions
│   └── batchRun7Cameras.m      7-camera preset
│
├── ParameterTesting/       ► calibration, validation, diagnostics
│   ├── getNormConstants.m           CF3 utopia/nadir constants
│   ├── buildNormalisationSchedule.m build them via a calibration pass
│   ├── buildNormFromBatch.m         derive them from existing CF1/CF2 runs
│   ├── sweepPopVsGenerations.m      population vs generations sweep
│   ├── feasibilityMap.m             config-independent coverage feasibility
│   ├── checkCameraFOV.m             does simulated FOV match the real lens?
│   ├── sanityCheckCoverage.m        per-camera / per-point coverage check
│   └── test*.m                      equivalence + range regression tests
│
├── GA_Core/                ► the genetic algorithm engine
│   ├── RunGA.m                  main loop
│   ├── initialPopulation.m      guided seeding from section centres
│   ├── Tournament.m, RouletteWheelSelection.m
│   ├── DoublePointCrossover.m, Mutate.m
│   ├── SortPopulation.m, fixPoorCameras.m
│
├── CostFunctions/          ► the objective
│   ├── combinedCostFunction.m   CF3 entry point
│   ├── cf3Terms.m               scaling + weighting (single source of truth)
│   ├── resUncertainty*.m        CF1
│   ├── dynamicOcclusion*.m      CF2
│   └── visibility / triangulability helpers
│
├── Geometry/               ► cameras, pyramids, target space, projection
├── Setup/                  ► setupProblem, setupGAparams,
│                             setupHardwareSpecs, setupCostParams
├── Plotting/               ► plotGARuns + every plot helper
├── Analysis/               ► saving, loading, and the analysis entry points
├── Sensitivity/            ► spacing sensitivity + OptiTrack chromosome
│
├── Results/                ► run outputs (TRACKED in git)
├── figures/                ► generated PDFs
├── Reference/              ► Rahimian et al. Python implementation
└── _Archive/               ► superseded work, OFF the MATLAB path
```

### How a run fits together

```
runCameraOptimiser  (Run/)
   ├─ setupProblem / setupHardwareSpecs / setupCostParams / setupGAparams   (Setup/)
   ├─ generateTargetSpace, generateSectionCentres, setupCameras             (Geometry/)
   ├─ RunGA                                                                (GA_Core/)
   │     └─ each candidate scored by combinedCostFunction → cf3Terms       (CostFunctions/)
   ├─ visualizeCameraCoverage, plotResults                                 (Plotting/)
   └─ saveResults → Results/<N>Cams/ + Results/Logs/                       (Analysis/)
```

---

## 8. Conventions worth knowing

- **Never use `fileparts(mfilename('fullpath'))` to find the project root.**
  From inside a subfolder that gives you the subfolder. Call
  `projectRoot = addProjectPaths();` instead — it returns the root regardless of
  where the calling file lives, so files can be moved without breaking paths.
- `addProjectPaths()` is safe to call repeatedly and is called by every entry
  point, so scripts work even if you forgot `START_HERE`.
- `_Archive/` is deliberately **off** the path. Nothing in it is used, and
  `_Archive/Results_pre-bugfix/` in particular predates the FOV and
  normalisation fixes — those runs are **not valid results** and nothing falls
  back to them.
- The GA is stochastic. Repeat runs (`NumRepeats`) and compare distributions,
  not single runs.

## 9. Troubleshooting

| Symptom | Cause |
|---|---|
| `Unrecognized function 'CentralCamera'` | RVC3 not on the path — run `START_HERE` |
| `Unrecognized function 'RunGA'` | code folders not on the path — run `START_HERE` |
| `Log file not found: .../Results/Logs/...` | no runs on this machine yet — pull them from the simulation machine |
| Wrong `resUncertainty` picked up | `OptiTrackConfig/` is on the path at the same time; re-run `START_HERE` to reset |
| Run is extremely slow | `spacing` too fine — it drives the number of evaluation points |
