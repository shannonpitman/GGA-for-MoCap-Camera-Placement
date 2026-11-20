function cfg = runConfig(preset, varargin)
%RUNCONFIG  Every setting that changes from run to run, in one place.
%
%   cfg = runConfig()                          % 'optitrack_lab'
%   cfg = runConfig('lowcost_tripod')
%   cfg = runConfig('optitrack_lab', 'Volume', [-3 3; -3 3; 0 3])
%
%   Presets
%     'optitrack_lab'   Lab OptiTrack rig: 8 x 8 x 4 m capture volume,
%                       cameras on walls, ceiling or tripods.
%     'lowcost_tripod'  Low-cost cameras on tripods only.
%     'openmv_2m'       OpenMV RT1062 cameras (single lens, 8 m IR range) on
%                       tripods only (head height 0.75-1.75 m) around the
%                       measured flight space: 5 m (x, between the walls
%                       holding the installed OptiTrack cameras) by 6 m (y,
%                       along those walls), targets 0.2-2 m high (the
%                       marker board's realistic sweep), 1.5 m tripod band
%                       on every side (x +-4, y +-4.5); UAV uniform grid
%                       6 x 7 x 3 = 126 points, 6-12 cameras.
%     'optitrack_tripod' OptiTrack cameras in the same flight space on the
%                       same tripods, 6-7 cameras: the like-for-like
%                       reference for openmv_2m. Kept apart from the
%                       full-room lab runs (Results/optitrack_tripod/).
%     'openmv_cab_2m'   OpenMV in the same flight space (0.2-2 m) with three
%                       kinds of mounting (2026-11-20): floor tripods
%                       (heads 0.5-2.0 m) in the band in front of the
%                       cabinets on the two x sides; tripods on the floor or
%                       a table (0.5-2.75 m, table ~0.75 m assumed) at the
%                       two y ends; tripods on the 0.9 m cabinets (0.75 m
%                       deep, 5.5 m long, centred on those walls), heads
%                       1.4-2.9 m. 6-12 cameras.
%     'openmv_cab_4m'   As openmv_cab_2m for targets 0.2-4 m (5 layers, 210
%                       points).
%     'openmv_floor_2m' Control for openmv_cab_2m: floor tripods only, heads
%                       0.5-2.0 m, in the 1.5 m band on every side.
%     'openmv_4m'       OpenMV on tripods in the full 8 x 8 x 4 m lab volume,
%                       6-10 cameras (first sweep, 2026-11-16). 'openmv_lab'
%                       is its old name and still maps to it.
%
%   Any field below can be overridden with a name-value pair. Unknown names
%   are an error, so a typo cannot silently fall back to a default.
%   buildRunSpecs(cfg, ...) turns a cfg into the specs/problem/params that
%   RunGA needs; batchRunGA and runCameraOptimiser both take their defaults
%   from here.

    if nargin < 1 || isempty(preset)
        preset = 'optitrack_lab';
    end

    %% Defaults shared by every preset
    cfg.Preset = preset;

    % Hardware profile (see setupHardwareSpecs)
    cfg.Hardware = 'optitrack';

    % Experimental design (batchRunGA sweeps these)
    cfg.CameraRange   = [6 7];
    cfg.CostFunctions = [1 2 3];
    cfg.TargetTypes   = [1 2];        % 1 = UAV (full volume), 2 = UGV (floor slab)
    cfg.GridModes     = [1 2];        % 1 = uniform, 2 = centre-weighted
    cfg.Spacings      = 1.0;          % x-y grid spacing [m]
    cfg.NumRepeats    = 5;
    cfg.SkipWarmStart = false;

    % Workspace
    cfg.Volume        = [-4 4; -4 4; 0 4];   % capture volume [m]
    cfg.ZSpacing      = [];                  % UAV grid z spacing [m]; [] = x-y spacing
    cfg.UGV_MaxHeight = 0.5;                 % UGV slab height [m]
    cfg.UGV_ZSpacing  = 0.25;                % UGV slab z spacing [m]

    % Camera search box [x y z rx ry rz]: position [m] and rotation vector
    % (exponential map, |r| <= pi, see genesToRotm). Positions are further
    % restricted to the mountable regions in cfg.Mount.
    cfg.CamLowerBounds = [-5 -4.5 0   -pi -pi -pi];
    cfg.CamUpperBounds = [ 5  4.5 4.8  pi  pi  pi];

    % Mountable camera positions (see mountRegions):
    %   Model 'walls_ceiling_tripod' | 'walls_ceiling' | 'tripod' | 'box'
    %   TripodHeight [zMin zMax] tripod head height range [m]
    %   Regions: custom region struct array; overrides Model when non-empty
    cfg.Mount = struct('Model', 'walls_ceiling_tripod', ...
                       'TripodHeight', [0.3 2.5], 'Regions', []);

    % Cost function
    cfg.Weights = [0.5 0.5];                 % [resolution, occlusion]

    % GA
    % 150 since 2026-11: the warm-start A/B study (Results/Tuning/ABWarmStart)
    % showed runs still improving at generation 100. Every run also records
    % its generation-SnapshotGeneration best, so 100 vs 150 is measured from
    % the same run (the GA does not depend on MaxGenerations).
    cfg.MaxGenerations  = 150;
    cfg.SnapshotGeneration = 100;
    cfg.PopulationSize  = [];                % [] = numCams * genesPerCam * 10
    cfg.PopulationScale = 10;
    cfg.CrossoverFraction = 1;
    % Mutation, tuned with ParameterTesting/tuneMutation (2026-10-02, simulation
    % machine, Results/Tuning/tuneMutation_20261002_180137): best validated
    % setting "tuned #1" (mu 0.294, 0.165 m, 3.3 deg; median J 0.576, IQR 0.012
    % vs 0.629 / 0.056 for the previous 0.5 / 0.1 / 0.1). mu rounded to 0.3.
    cfg.MutationRate     = 0.3;              % per-gene probability
    cfg.MutationSigmaPos = 0.165;            % position step s.d. [m]
    cfg.MutationSigmaRot = deg2rad(3.3);     % rotation-vector step s.d. [rad]
    cfg.TournamentSize  = 3;

    % Warm start (see WarmStart/selectWarmSeeds). Warm runs inject seeds from
    % Results/SeedArchive/<preset>/<condition>.mat, which every saved run
    % (cold and warm) adds to.
    %   WarmSeedMode      'top_n' (best distinct archived solutions) or
    %                     'best_mutant' (best + mutated copies, the old scheme)
    %   WarmSeedCount     seeds injected; [] = WarmSeedFraction * population
    %   WarmSeedMinDist   diversity threshold for 'top_n' [m], 0 = off
    %   ArchiveFraction   share of each run's final population archived
    cfg.WarmSeedMode     = 'top_n';
    cfg.WarmSeedCount    = 10;        % A/B study 2026-11: top_n x 10 best
    cfg.WarmSeedFraction = 0.05;
    cfg.WarmSeedMinDist  = 0.25;
    cfg.ArchiveFraction  = 0.10;
    %   FreezeWarmPool    true: a batch's warm runs of one condition all seed
    %                     from the archive as it stood at that condition's
    %                     first warm run (independent repeats). false: each
    %                     warm run also sees the batch's earlier warm runs.
    cfg.FreezeWarmPool   = true;

    %% Preset-specific values
    switch lower(preset)
        case 'optitrack_lab'
            % defaults above

        case 'lowcost_tripod'
            cfg.Hardware = 'lowcost';
            cfg.Mount.Model = 'tripod';

        case {'openmv_2m', 'optitrack_tripod'}
            % Measured flight space (2026-11-17): 5 m x 6 m, board swept
            % from 0.2 m to 2 m; 1.5 m of floor for tripods on every side
            % before the walls/cabinets.
            cfg.Volume         = [-2.5 2.5; -3 3; 0.2 2.0];
            cfg.ZSpacing       = 0.9;               % layers at 0.2, 1.1, 2.0 m
            cfg.CamLowerBounds = [-4 -4.5 0   -pi -pi -pi];
            cfg.CamUpperBounds = [ 4  4.5 4.8  pi  pi  pi];
            cfg.Mount.Model    = 'tripod';
            cfg.Mount.TripodHeight = [0.75 1.75];  % user's tripods; 2 m possible but unstable [m]
            cfg.TargetTypes = 1;          % UAV
            cfg.GridModes   = 1;          % uniform
            if strcmpi(preset, 'openmv_2m')
                cfg.Hardware    = 'openmv_rt1062';
                cfg.CameraRange = 6:12;   % cameras available
            else
                cfg.CameraRange = [6 7];  % as the lab rig
            end

        case {'openmv_cab_2m', 'openmv_cab_4m', 'openmv_floor_2m'}
            cfg.Hardware       = 'openmv_rt1062';
            cfg.CameraRange    = 6:12;
            cfg.TargetTypes    = 1;
            cfg.GridModes      = 1;
            cfg.Volume         = [-2.5 2.5; -3 3; 0.2 2.0];
            cfg.ZSpacing       = 0.9;               % 0.2, 1.1, 2.0 m
            cfg.Mount.TripodHeight = [0.5 2.0];     % tripods remeasured 2026-11-20
            if strcmpi(preset, 'openmv_floor_2m')
                cfg.CamLowerBounds = [-4 -4.5 0   -pi -pi -pi];
                cfg.CamUpperBounds = [ 4  4.5 4.8  pi  pi  pi];
                cfg.Mount.Model    = 'tripod';
            else
                % room to the back of the cabinets (x) and the 1.5 m end bands (y)
                cfg.CamLowerBounds = [-4.75 -4.5 0   -pi -pi -pi];
                cfg.CamUpperBounds = [ 4.75  4.5 4.8  pi  pi  pi];
                cfg.Mount.Model    = 'tripod_cabinets';
                cfg.Mount.EndHeight = [0.5 2.75];   % floor tripod, or tripod on a ~0.75 m table
                cfg.Mount.Cabinet  = struct('XDepth', [4 4.75], 'YSpan', [-2.75 2.75], 'Height', 0.9);
                if strcmpi(preset, 'openmv_cab_4m')
                    cfg.Volume(3, :) = [0.2 4.0];
                    cfg.ZSpacing     = 0.95;        % 0.2, 1.15, 2.1, 3.05, 4.0 m
                end
            end

        case {'openmv_4m', 'openmv_lab'}
            cfg.Preset      = 'openmv_4m';         % results, archive and table under this name
            cfg.Hardware    = 'openmv_rt1062';
            cfg.Mount.Model = 'tripod';
            cfg.Mount.TripodHeight = [0.75 1.75];
            cfg.TargetTypes = 1;
            cfg.GridModes   = 1;
            cfg.CameraRange = 6:10;

        otherwise
            error('runConfig:UnknownPreset', ...
                'Unknown preset "%s". Use ''optitrack_lab'', ''openmv_2m'', ''optitrack_tripod'', ''openmv_cab_2m'', ''openmv_cab_4m'', ''openmv_floor_2m'', ''openmv_4m'' or ''lowcost_tripod''.', preset);
    end

    %% Name-value overrides
    if mod(numel(varargin), 2) ~= 0
        error('runConfig:BadOverrides', 'Overrides must be name-value pairs.');
    end
    for k = 1:2:numel(varargin)
        name = varargin{k};
        if ~isfield(cfg, name)
            error('runConfig:UnknownField', 'runConfig has no field "%s".', name);
        end
        cfg.(name) = varargin{k+1};
    end
end
