# UGV grid-change runbook

Use this after you stop the in-progress batch on **USER-PC**. The OneDrive
copy of the code has been updated; USER-PC needs to pick up those changes
before you do anything else.

## TL;DR

The legacy UGV runs used an isotropic 0.25 m grid (33 × 33 × 3 = 3267 points)
because `batchRunGA` clamped the UGV spacing to `UGV_MaxHeight/2`. UGV runs
were ~7× slower than UAV runs because of this. The code now uses
**anisotropic UGV spacing**: x-y honours the swept `Spacing`, z is fixed
at `UGV_ZSpacing` (default 0.25 m → 3 layers on a 0.5 m slab). Existing
UGV runs are no longer comparable to anything we'll produce going forward,
so they get archived and the master log gets pruned.

## What changed in the code

| File | Change |
|---|---|
| `Geometry/generateTargetSpace.m` | Accepts scalar OR `1×3` vector spacing. Backward-compatible (UAV side unaffected). |
| `batchRunGA.m` | New `UGV_ZSpacing` parameter (default 0.25 m). UGV branch builds `[spacing, spacing, UGV_ZSpacing]` instead of clamping isotropic. Spacing logged in `runLog.Spacing` is still the scalar x-y so old log entries keep their shape. |
| `runCameraOptimiser.m` | Same anisotropic logic for single-run usage. `UGV_zSpacing = 0.25` defined in the user-input block. |
| `spacingSensitivity_UGV.m` | Spacings swept are now `[0.25, 0.4, 0.5, 0.75, 1.0]` (x-y), with z fixed at 0.25 m. |
| `Sensitivity/spacingSensitivityCore.m` | New optional `opts.zSpacing` field. When set, anisotropic grids are built per sweep point. |
| `Analysis/archiveUGVRuns.m` | New utility to move UGV-flagged runs to `Results_UGV_fine_grid_archive/` and prune the master log. Dry-run by default. |

## Step-by-step

### 1. Stop the running batch cleanly on USER-PC

In the MATLAB Command Window on USER-PC:

```matlab
% Ctrl+C — interrupts whatever run is currently executing.
% Wait for MATLAB to return to the prompt, then:
close all;
clear;
```

The batch was at 61 / 160 (last completed UGV run logged at 02:06 +
177 min, so the next print is still mid-run). Ctrl+C will abort that
run; the 71-or-so runs already saved up to that point stay on disk.

### 2. Sync the OneDrive code changes to USER-PC

Whatever sync mechanism you normally use between OneDrive and USER-PC's
local copy (git pull, OneDrive sync, manual copy of the .m files) —
trigger it now. You need these files updated on USER-PC:

```
Geometry/generateTargetSpace.m
batchRunGA.m
runCameraOptimiser.m
spacingSensitivity_UGV.m
Sensitivity/spacingSensitivityCore.m
Analysis/archiveUGVRuns.m         (new file)
```

Quick check inside MATLAB on USER-PC:

```matlab
type batchRunGA.m | head -n 30        % should mention UGV_ZSpacing
which archiveUGVRuns                   % should resolve
```

### 3. Archive the existing UGV runs

```matlab
cd 'C:\Users\USER-PC\GGA-for-MoCap-Camera-Placement\Real-coded GGA\Optimal Camera Placement Genetic Algorithm'
addProjectPaths();

% Dry run — see what would move
plan = archiveUGVRuns();

% If the plan looks right, commit:
plan = archiveUGVRuns('DryRun', false);
```

The dry run prints every UGV entry with its log index, cost, timestamp,
and source path. The commit step:

1. Writes a timestamped backup of `Results/Logs/GGA_RunsLog.mat`
2. Moves each UGV `.mat` and `.txt` from `Results/<N>Cams/` to
   `Results_UGV_fine_grid_archive/<N>Cams/`
3. Saves a pruned log over the original

The backup name will be something like
`Results/Logs/GGA_RunsLog_pre_UGVarchive_20260511_120000.mat` — keep
that around until you're sure the new sweep + batch work.

### 4. Run the UGV grid-sensitivity sweep

```matlab
sweep = spacingSensitivity_UGV();
```

This:

- Loads the GA-best CF3 7-cam **UGV** chromosome from the (now pruned)
  log. If you archived every UGV CF3 run, it will fall back to the
  UAV-best chromosome with a warning — that's still a valid stand-in
  for evaluating how cost varies with x-y resolution.
- Builds anisotropic grids `[xy, xy, 0.25]` for `xy ∈ {0.25, 0.4, 0.5, 0.75, 1.0}` m
- Evaluates CF1, CF2, CF3 for the GA chromosome **and** the OptiTrack
  ad-hoc chromosome at each spacing
- Prints a recommendation table (coarsest spacing whose deviation from
  the 0.25 m reference stays under 2%)
- Saves `Results/Sensitivity/UGV/spacing_sweep_UGV_<timestamp>.mat`
- Writes the three comparison figures into `figures/Sensitivity/UGV/`

Look at the printed table and the deviation plot — pick the **largest**
x-y spacing whose CF3 deviation from the 0.25 m reference is under your
tolerance. Most likely candidates: 0.5 m or 0.75 m. 1.0 m matches UAV
exactly, but UGV is more sensitive to in-plane resolution because of how
the floor-slab points cluster, so verify before committing.

### 5. Relaunch the batch with the chosen spacing

If you decide on 0.5 m x-y for UGV (just an example):

```matlab
% Resume from the existing batch log so the UAV runs aren't redone:
batchRunGA('ResumeLog', 'Results/Logs/BatchLog_<the-timestamp-from-before>.mat', ...
           'Spacings', 0.5, ...        % NB: now the UGV x-y spacing
           'UGV_ZSpacing', 0.25, ...   % 3 layers
           'TargetTypes', 2);          % UGV only
```

Or start a fresh batch:

```matlab
batchRunGA('Spacings', 0.5, 'UGV_ZSpacing', 0.25, 'TargetTypes', 2);
```

**Expected speedup at 0.5 m UGV x-y:** 17 × 17 × 3 = 867 points
(vs the legacy 3267). Roughly a **3-4× faster** UGV run, plus less
random GA noise from over-sampled floor points. At 1.0 m x-y it's
9 × 9 × 3 = 243 points — basically as fast as UAV.

## Sanity checks to do before kicking off a 100-hour batch

1. `which generateTargetSpace` on USER-PC resolves to the updated file
   (i.e. it accepts vector spacing).
2. Run **one** UGV sanity invocation:
   ```matlab
   batchRunGA('CameraRange', 7, 'CostFunctions', 3, 'TargetTypes', 2, ...
              'GridModes', 1, 'Spacings', 0.5, 'UGV_ZSpacing', 0.25, ...
              'NumRepeats', 1, 'SkipWarmStart', true, 'MaxGenerations', 30, ...
              'DryRun', false);
   ```
   It should print `UGV slab: z in [0, 0.5] m, z-spacing fixed at 0.25 m (3 layers)`
   in the pre-flight summary and the run text-summary should show
   `Number of Target Points: 867` (for spacing 0.5).
3. After the sanity run, peek at the new log entry —
   `runLog(end).Spacing` should equal **0.5** (scalar, not a vector),
   `TargetType == 2`, and `NumTargetPoints == 867`.

If any of those fail, do **not** start the long batch — debug first
(or message me with the error).

## Why this is defensible for the thesis

You'll be able to write:

> *Grid-spacing sensitivity for the UGV target space was assessed by
> evaluating the GA-best CF3 7-camera configuration over in-plane
> spacings of 0.25 – 1.0 m, with the slab-thickness z-spacing held
> fixed at 0.25 m (3 layers). The combined-cost deviation from the
> 0.25 m reference grid remained below X % at Y m, motivating
> Y m as the production in-plane spacing for the full UGV batch.*

That's the standard discretisation-convergence argument and will hold
up to examiner pushback far better than the previous hardcoded 0.25 m.
