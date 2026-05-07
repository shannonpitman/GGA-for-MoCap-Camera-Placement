function sweep = spacingSensitivity_UGV()
%SPACINGSENSITIVITY_UGV  Grid-spacing sensitivity sweep for UGV target space.
%
%   Same comparison as spacingSensitivity_UAV but with the floor-slab
%   target space used for UGV optimisation: the workspace volume's z range
%   is capped at UGV_maxHeight (default 0.5 m), and the spacing sweep
%   tightens its upper bound to UGV_maxHeight/2 = 0.25 m so that at least
%   two layers of target points fit in the slab.
%
%   USAGE:
%       sweep = spacingSensitivity_UGV()
%
%   NOTES:
%     * The GA-best CF3 7-cam config loaded here is filtered to runs with
%       TargetType=2 (UGV). If you have not yet run a UGV batch, this
%       script will error — fall back to runCameraOptimiser with
%       targetType=2 first, or temporarily edit the call to
%       loadBestCF3Config to relax the filter.
%     * Spacings below the marker spacing of 0.05 m (typical UGV) are
%       not added by default; extend opts.spacings if the application
%       has finer markers.
%
%   See also: spacingSensitivity_UAV, spacingSensitivityCore.

    addProjectPaths();

    %% UGV-specific inputs
    UGV_maxHeight = 0.5;                         % m, floor slab height
    opts.modeTag      = 'UGV';
    opts.targetType   = 2;                       % UGV: floor slab
    opts.volume       = [-4 4; -4 4; 0 UGV_maxHeight];
    % Spacing sweep capped at UGV_maxHeight/2 to ensure >=2 z-layers
    opts.spacings     = [0.05, 0.075, 0.1, 0.15, 0.2, UGV_maxHeight/2];
    opts.targetMode   = 1;                       % uniform grid
    opts.numCams      = 7;
    opts.weightUnc    = 0.5;
    opts.weightOcc    = 0.5;
    opts.tolerancePct = 2;

    %% Configurations
    try
        [gaChrom, bestRun] = loadBestCF3Config(opts.numCams, opts.targetType, opts.targetMode);
        opts.configs(1).name       = sprintf('GA-best CF3 UGV (logged cost=%.5f)', bestRun.BestCost);
        opts.configs(1).chromosome = gaChrom;
    catch ME
        if strcmp(ME.identifier, 'loadBestCF3Config:NoMatch')
            warning(['No UGV-tagged CF3 runs found in master log. ' ...
                     'Falling back to the UAV-best chromosome — re-run a ' ...
                     'UGV CF3 batch and re-evaluate this sweep when available.']);
            [gaChrom, bestRun] = loadBestCF3Config(opts.numCams, 1, opts.targetMode);
            opts.configs(1).name = sprintf('GA-best CF3 UAV-fallback (cost=%.5f)', ...
                                           bestRun.BestCost);
            opts.configs(1).chromosome = gaChrom;
        else
            rethrow(ME);
        end
    end

    opts.configs(2).name       = 'OptiTrack ad-hoc';
    opts.configs(2).chromosome = buildOptiTrackChromosome();

    %% Run
    sweep = spacingSensitivityCore(opts);
end
