function sweep = spacingSensitivity_UAV()
%SPACINGSENSITIVITY_UAV  Grid-spacing sensitivity sweep for UAV target space.
%
%   Compares the GA-best CF3 7-cam configuration against the OptiTrack
%   ad-hoc rig over the full UAV flight envelope, sweeping grid spacings
%   from 0.1 m (microUAV marker spacing) to 2.0 m.
%
%   USAGE:
%       sweep = spacingSensitivity_UAV()
%
%   See also: spacingSensitivity_UGV, spacingSensitivityCore,
%             recommendSpacing, plotSpacingSensitivity.

    addProjectPaths();

    %% UAV-specific inputs
    opts.modeTag      = 'UAV';
    opts.targetType   = 1;                       % UAV: full flight volume
    opts.volume       = [-4 4; -4 4; 0 4];       % MS.G flight envelope
    opts.spacings     = [0.1, 0.15, 0.25, 0.4, 0.6, 1.0, 1.5, 2.0];
    opts.targetMode   = 1;                       % uniform grid
    opts.numCams      = 7;
    opts.weightUnc    = 0.5;
    opts.weightOcc    = 0.5;
    opts.tolerancePct = 2;                       % flag deviations > 2% as too coarse

    %% Configurations
    [gaChrom, bestRun] = loadBestCF3Config(opts.numCams, opts.targetType, opts.targetMode);
    opts.configs(1).name       = sprintf('GA-best CF3 (logged cost=%.5f)', bestRun.BestCost);
    opts.configs(1).chromosome = gaChrom;

    opts.configs(2).name       = 'OptiTrack ad-hoc';
    opts.configs(2).chromosome = buildOptiTrackChromosome();

    %% Run
    sweep = spacingSensitivityCore(opts);
end
