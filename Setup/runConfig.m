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
    cfg.UGV_MaxHeight = 0.5;                 % UGV slab height [m]
    cfg.UGV_ZSpacing  = 0.25;                % UGV slab z spacing [m]

    % Camera search box [x y z alpha beta gamma]. Positions are further
    % restricted to the mountable regions in cfg.Mount.
    cfg.CamLowerBounds = [-5 -4.5 0   -pi -pi/2 -pi];
    cfg.CamUpperBounds = [ 5  4.5 4.8  pi  pi/2  pi];

    % Mountable regions (filled in by the preset below; see projectToMount)
    cfg.Mount = struct('Regions', {{}});

    % Cost function
    cfg.Weights = [0.5 0.5];                 % [resolution, occlusion]

    % GA
    cfg.MaxGenerations  = 100;
    cfg.PopulationSize  = [];                % [] = numCams * genesPerCam * 10
    cfg.PopulationScale = 10;
    cfg.CrossoverFraction = 1;
    cfg.MutationRate    = 0.5;               % per-gene probability
    cfg.MutationSigma   = 0.1;
    cfg.TournamentSize  = 3;

    %% Preset-specific values
    switch lower(preset)
        case 'optitrack_lab'
            % defaults above

        case 'lowcost_tripod'
            cfg.Hardware = 'lowcost';

        otherwise
            error('runConfig:UnknownPreset', ...
                'Unknown preset "%s". Use ''optitrack_lab'' or ''lowcost_tripod''.', preset);
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
