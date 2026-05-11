function sweep = spacingSensitivity_UGV(varargin)
%SPACINGSENSITIVITY_UGV  Grid-spacing sensitivity sweep for UGV target space.
%
%   Sweeps the in-plane (x-y) grid spacing on the floor slab while holding
%   the z (slab-thickness) spacing fixed. This decouples the slab-layer
%   count from the x-y resolution so the GA-cost surface can be evaluated
%   at coarse in-plane grids comparable to the UAV runs, without losing
%   the 3-layer sampling of the 0.5 m slab.
%
%   Defaults reproduce the original UGV slab (UGV_maxHeight = 0.5 m, z step
%   0.25 m -> z = [0, 0.25, 0.5], 3 layers) and sweep x-y from 0.25 m
%   (matching the legacy isotropic grid) up to 1.0 m (matching UAV).
%
%   USAGE:
%       sweep = spacingSensitivity_UGV()
%       sweep = spacingSensitivity_UGV('TolerancePct', 5)
%       sweep = spacingSensitivity_UGV('Spacings', [0.25 0.5 0.75 1.0], ...
%                                      'TolerancePct', 10)
%
%   NOTES:
%     * Loads the GA-best CF3 7-cam UGV chromosome from GGA_RunsLog.mat.
%       Falls back to the UAV-best chromosome with a warning if no UGV
%       CF3 runs are logged yet.
%     * Pass the recommended x-y spacing back into batchRunGA via the
%       Spacings parameter; UGV_ZSpacing remains 0.25 m unless you have
%       a separate reason to vary z.
%
%   See also: spacingSensitivity_UAV, spacingSensitivityCore,
%             recommendSpacing, plotSpacingSensitivity.

    addProjectPaths();

    %% Optional overrides
    p = inputParser;
    addParameter(p, 'Spacings',     [0.25, 0.4, 0.5, 0.75, 1.0], @isnumeric);
    addParameter(p, 'ZSpacing',     0.25,                        @isnumeric);
    addParameter(p, 'TolerancePct', 2,                           @isnumeric);
    addParameter(p, 'UGVmaxHeight', 0.5,                         @isnumeric);
    parse(p, varargin{:});
    args = p.Results;

    %% UGV-specific inputs
    opts.modeTag      = 'UGV';
    opts.targetType   = 2;                        % UGV: floor slab
    opts.volume       = [-4 4; -4 4; 0 args.UGVmaxHeight];
    % x-y spacings to test. 0.25 m anchors to the legacy isotropic grid
    % so the recommendation pass can quote deviation relative to it.
    opts.spacings     = args.Spacings;
    opts.zSpacing     = args.ZSpacing;            % fixed z step (3 layers on 0.5 m slab)
    opts.targetMode   = 1;                        % uniform grid
    opts.numCams      = 7;
    opts.weightUnc    = 0.5;
    opts.weightOcc    = 0.5;
    opts.tolerancePct = args.TolerancePct;        % deviation threshold (%)

    %% Configurations
    try
        [gaChrom, bestRun] = loadBestCF3Config(opts.numCams, opts.targetType, opts.targetMode);
        opts.configs(1).name       = sprintf('GA-best CF3 UGV (logged cost=%.5f)', bestRun.BestCost);
        opts.configs(1).chromosome = gaChrom;
    catch ME
        if strcmp(ME.identifier, 'loadBestCF3Config:NoMatch')
            warning(['No UGV-tagged CF3 runs found in master log. ' ...
                     'Falling back to the UAV-best chromosome - re-run a ' ...
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
