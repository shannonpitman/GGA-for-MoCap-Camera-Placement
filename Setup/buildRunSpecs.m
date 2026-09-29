function [specs, problem, params] = buildRunSpecs(cfg, numCams, costFunctionType, targetType, gridMode, spacing, varargin)
%BUILDRUNSPECS  specs, problem and GA params for one run, from a runConfig.
%
%   [specs, problem, params] = buildRunSpecs(cfg, numCams, cf, tt, gm, spacing)
%
%   The single place where a run's target space, hardware, cost parameters,
%   search bounds and GA settings are assembled. batchRunGA,
%   runCameraOptimiser, the normalisation builder and the sensitivity tools
%   all go through here, so every cost is evaluated on an identical setup.
%
%   cfg              - struct from runConfig
%   costFunctionType - 1 = resolution, 2 = occlusion, 3 = combined
%   targetType       - 1 = UAV (full volume), 2 = UGV (floor slab)
%   gridMode         - 1 = uniform, 2 = centre-weighted
%   spacing          - x-y grid spacing [m]
%
%   Name-value: 'UseNormTable' (default true). The normalisation builders
%   pass false so they never read the table they are writing.

    p = inputParser;
    addParameter(p, 'UseNormTable', true, @islogical);
    parse(p, varargin{:});

    %% Target volume and grid spacing
    volume = cfg.Volume;
    if targetType == 2
        % UGV floor slab: x-y honour the swept spacing, z is fixed so the
        % slab keeps its layers when x-y is coarsened.
        volume(3, :) = [0, cfg.UGV_MaxHeight];
        zSpacing = min(cfg.UGV_ZSpacing, cfg.UGV_MaxHeight);
        targetSpacing = [spacing, spacing, zSpacing];
    else
        targetSpacing = spacing;
    end

    %% Specs
    specs = setupHardwareSpecs(numCams, cfg.Hardware);
    specs.RunConfig = cfg;
    specs.Parameterisation = 'rotvec';   % orientation genes = rotation vector
    specs.warmStart = false;
    specs.warmChromosomes = [];

    specs.WeightUncertainty = cfg.Weights(1);
    specs.WeightOcclusion   = cfg.Weights(2);

    specs.TargetType = targetType;
    specs.TargetMode = gridMode;
    specs.Target     = generateTargetSpace(volume, gridMode, targetSpacing);
    specs.NumPoints  = size(specs.Target, 1);
    specs.spacing    = spacing;
    if targetType == 2
        specs.spacingZ = zSpacing;
    end

    specs.SectionCentres = generateSectionCentres(numCams, volume);
    specs.MountRegions   = mountRegions(cfg);
    specs.UseNormTable = p.Results.UseNormTable;
    specs = setupCostParams(specs);

    %% Problem and GA parameters
    problem = setupProblem(numCams, costFunctionType, cfg.CamUpperBounds, cfg.CamLowerBounds);

    if isempty(cfg.PopulationSize)
        popSize = numCams * numel(cfg.CamLowerBounds) * cfg.PopulationScale;
    else
        popSize = cfg.PopulationSize;
    end
    params = setupGAparams(cfg.MaxGenerations, popSize, cfg);
end
