# Results folders

Every GA run is saved by `saveResults` into the folder for its `runConfig`
preset and logged in `Logs/GGA_RunsLog.mat` with a `Preset` field.
`loadGARuns` returns only `optitrack_lab` runs unless you ask for another
preset (`loadGARuns('Preset', 'openmv_2m')`), so the sweeps never mix in an
analysis. Each preset also has its own CF3 normalisation table
(`normTable_<preset>.mat`) and seed archive (`SeedArchive/<preset>/`).

All target grids are the UAV uniform grid at 1 m x-y spacing unless stated.
"Flight space" means the measured flyable area (2026-11-17): 5 m in x
(between the two walls holding the installed OptiTrack cameras) by 6 m in y
(along those walls), centred on the origin.

## Camera-placement sweeps

| Folder | Preset | Hardware | Mounting | Capture volume / targets | Cameras | Purpose |
|---|---|---|---|---|---|---|
| `6Cams/`, `7Cams/` | `optitrack_lab` | OptiTrack (narrow/wide lens) | walls, ceiling and tripods (0.3-2.5 m) | full lab, 8 x 8 x 4 m; UAV and UGV targets, uniform and centre-weighted grids | 6, 7 | Main GA results for the thesis (GA vs the manually posed rig). Norm table `normTable.mat`. |
| `openmv_4m/` | `openmv_4m` (logged as `openmv_lab`) | OpenMV RT1062 | floor tripods 0.75-1.75 m, 0.5-1 m band just outside the footprint | 8 x 8 x 4 m (405 points) | 6-10 | First OpenMV sweep (2026-11-16). Coverage limited at the footprint edge. |
| `openmv_2m/` | `openmv_2m` | OpenMV RT1062 | floor tripods 0.75-1.75 m, 1.5 m band on every side | flight space, 0.2-2 m (126 points) | 6-12 | OpenMV in the real flight space. |
| `optitrack_tripod/` | `optitrack_tripod` | OptiTrack | as `openmv_2m` | as `openmv_2m` | 6, 7 | Like-for-like OptiTrack reference for `openmv_2m`. Not comparable with the full-room lab runs. |
| `openmv_cab_2m/` | `openmv_cab_2m` | OpenMV RT1062 | floor tripods 0.5-2.0 m in front of the cabinets (x sides); floor or table tripods 0.5-2.75 m at the y ends; tripods on the 0.9 m cabinets, heads 1.4-2.9 m | flight space, 0.2-2 m (126 points) | 6-12 | Effect of raised (cabinet/table) mounting. |
| `openmv_floor_2m/` | `openmv_floor_2m` | OpenMV RT1062 | floor tripods 0.5-2.0 m, 1.5 m band on every side | flight space, 0.2-2 m (126 points) | 6-12 | Control for `openmv_cab_2m`: same tripod range, no cabinets or tables. |
| `openmv_cab_4m/` | `openmv_cab_4m` | OpenMV RT1062 | as `openmv_cab_2m` | flight space, 0.2-4 m (5 layers, 210 points) | 6-12 | Supervisor's request: full 4 m height with cabinet mounting. |
| `openmv_cab_2m_fp/`, `openmv_floor_2m_fp/`, `openmv_cab_4m_fp/` | `*_fp` | OpenMV RT1062 | as the preset without `_fp`, plus the tripod footprint: legs 0.3 m from the centre kept clear of walls, cabinet fronts and cabinet-top edges (cabinet tripods x 4.3-4.45 m, y -2.45 to 2.45 m); any two tripods at least 0.6 m apart; every camera level or tilted up, never down (`Mount.MinElevation = 0`) | as without `_fp` | 6-12 | Buildable rigs. The OpenMV screws into the ball head through its back, so the head's 90-degree slot stops it at level: it can look level or up only. The earlier OpenMV sweeps allowed any tilt (the cab_2m cabinet cameras look 25-30 degrees down), put two tripods 0.22 m apart and put tripods on the cabinets' back edge; repairing those rigs afterwards (`legaliseRig`) roughly doubles their occlusion even before the tilt limit, so these sweeps re-optimise with all three constraints. |

Cabinets (measured 2026-11-20): 0.9 m tall, 0.75 m deep, 5.5 m long, on the
two OptiTrack walls behind the 1.5 m floor band (|x| 4-4.75 m), assumed
centred along those walls (y -2.75 to 2.75 m). The y-end tables are assumed
to be about 0.75 m tall.

Comparisons that isolate one change:
- cabinets vs floor: `openmv_cab_2m` vs `openmv_floor_2m`
- tripod range 0.75-1.75 vs 0.5-2.0 m: `openmv_floor_2m` vs `openmv_2m`
- 2 m vs 4 m volume (cabinet mounting): `openmv_cab_2m` vs `openmv_cab_4m`
- OpenMV vs OptiTrack: `openmv_2m` vs `optitrack_tripod`

Compare across presets on raw resolution uncertainty [mm], dynamic occlusion
[deg] and coverage (`costBreakdown`, `perTargetCoverage`,
`Plotting/plotFlightSpaceSweep`), not on CF3 cost: each preset and camera
count has its own normalisation.

## Other folders

| Folder | Contents |
|---|---|
| `Logs/` | `GGA_RunsLog.mat` (every run, all presets) and one `BatchLog_<date>.mat` per batch (schedule, settings, status; used to resume). |
| `SeedArchive/<preset>/` | Warm-start pool per condition: each run's best and the top 10% of its final population. |
| `Analysis/` | Text reports (`compareGenerations` 100 vs 150 generations, re-run summaries). |
| `Tuning/` | GA tuning studies: `tuneMutation`, `ABParameterisation` (rotation-vector vs Euler), `ABWarmStart` (warm-start seeding A/B study), `MacPilot`. |
| `Sensitivity/` | Mounting, orientation and grid-spacing sensitivity of the lab rigs. |
| `Sweep_PopVsGen/` | Population size vs generations sweep (July 2026). |
| `Feasibility/` | Configuration-independent coverage feasibility map. |
| `normTable*.mat` | CF3 utopia/nadir tables; `normTable.mat` is the lab table, `_backup_<date>` the superseded ones. |
| `Handover/` | Not tracked: local hand-over material. |
