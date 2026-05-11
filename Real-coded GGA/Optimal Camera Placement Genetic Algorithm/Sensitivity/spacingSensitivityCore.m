function sweep = spacingSensitivityCore(opts)
%SPACINGSENSITIVITYCORE  Run the spacing sensitivity sweep (UAV or UGV).
%
%   sweep = spacingSensitivityCore(opts) runs the grid-spacing sensitivity
%   sweep for the configurations and parameters supplied in opts.
%
%   opts fields (all required unless noted):
%       modeTag        'UAV' or 'UGV' — used in filenames + figure dirs
%       targetType     1 = UAV (full volume), 2 = UGV (floor slab)
%       volume         3x2 [xLo xHi; yLo yHi; zLo zHi] (m)
%       spacings       row vector of grid spacings to sweep (m). When
%                      opts.zSpacing is set these are interpreted as x-y
%                      only; otherwise they are applied isotropically.
%       targetMode     1 = uniform grid, 2 = inverse-CDF concentrated
%       numCams        camera count (must match the configurations supplied)
%       weightUnc      CF3 weight on resolution uncertainty
%       weightOcc      CF3 weight on dynamic occlusion
%       configs        struct array with fields .name and .chromosome
%       zSpacing       (optional) scalar z-axis spacing held fixed across
%                      the sweep. When supplied, generateTargetSpace is
%                      called with [s s zSpacing] so the in-plane grid is
%                      swept independently of slab-layer count.
%       tolerancePct   (optional, [] to skip) deviation tolerance (%) for
%                      recommendSpacing pass + plot annotations
%
%   Returns the sweep struct that's also written to
%   Results/Sensitivity/<modeTag>/spacing_sweep_<modeTag>_<timestamp>.mat.
%
%   Side effects:
%     - Saves three figures to figures/Sensitivity/<modeTag>/
%     - If tolerancePct is supplied, prints recommendation table.

    %% Validate inputs
    required = {'modeTag','targetType','volume','spacings','targetMode', ...
                'numCams','weightUnc','weightOcc','configs'};
    for k = 1:numel(required)
        assert(isfield(opts, required{k}), ...
            'spacingSensitivityCore:MissingField', ...
            'opts.%s is required.', required{k});
    end
    if ~isfield(opts, 'tolerancePct'), opts.tolerancePct = []; end
    if ~isfield(opts, 'zSpacing'),     opts.zSpacing     = []; end

    addProjectPaths();
    projectRoot = fileparts(fileparts(mfilename('fullpath')));   % .../<root>

    %% Hardware spec (used for ALL configurations to isolate spacing effect)
    specs = setupHardwareSpecs(opts.numCams);
    specs.WeightUncertainty = opts.weightUnc;
    specs.WeightOcclusion   = opts.weightOcc;
    specs.TargetType        = opts.targetType;
    specs.TargetMode        = opts.targetMode;

    configs = opts.configs;
    spacings = opts.spacings;
    nS = numel(spacings);
    nC = numel(configs);
    cfNames  = {'CF1_resUncert', 'CF2_dynOccl', 'CF3_combined'};
    cfLabels = {'Resolution Uncertainty (raw mean)', ...
                'Dynamic Occlusion (raw mean)', ...
                'Combined CF3 (normalised)'};
    nF = numel(cfNames);

    blank = struct('config','', 'cf','', 'spacing',NaN, ...
                   'cost',NaN, 'numPoints',NaN, 'evalTime',NaN);
    results = repmat(blank, nC*nS*nF, 1);
    idx = 0;

    %% Banner
    fprintf('\n======================================================\n');
    fprintf('  Spacing sensitivity sweep — %s, %d cams\n', ...
        upper(opts.modeTag), opts.numCams);
    fprintf('======================================================\n');
    fprintf('Volume   : x=[%g %g], y=[%g %g], z=[%g %g] m\n', ...
        opts.volume(1,1), opts.volume(1,2), opts.volume(2,1), opts.volume(2,2), ...
        opts.volume(3,1), opts.volume(3,2));
    fprintf('Spacings : '); fprintf('%.3g ', spacings);
    if isempty(opts.zSpacing)
        fprintf('m (isotropic)\n');
    else
        fprintf('m (x-y); z held fixed at %.3g m\n', opts.zSpacing);
    end
    fprintf('Configs  :\n');
    for c = 1:nC, fprintf('  [%d] %s\n', c, configs(c).name); end
    fprintf('Hardware : %dx%d res, focal %.1f/%.1f mm, range %.1f/%.1f m\n', ...
        specs.Resolution(1), specs.Resolution(2), ...
        specs.Focal*1e3, specs.FocalWide*1e3, specs.Range, specs.RangeWide);
    fprintf('------------------------------------------------------\n');

    sweepStart = tic;

    %% Sweep
    for c = 1:nC
        chrom = configs(c).chromosome;
        fprintf('\n--- Config %d/%d: %s ---\n', c, nC, configs(c).name);

        for s = 1:nS
            spacing = spacings(s);

            if isempty(opts.zSpacing)
                gridSpacing = spacing;                            % isotropic
            else
                gridSpacing = [spacing, spacing, opts.zSpacing];  % anisotropic
            end

            specs.Target    = generateTargetSpace(opts.volume, opts.targetMode, gridSpacing);
            specs.NumPoints = size(specs.Target, 1);
            specs.spacing   = spacing;                % scalar x-y for logging
            if ~isempty(opts.zSpacing)
                specs.spacingZ = opts.zSpacing;
            end
            specs           = setupCostParams(specs);

            [cameras, CamCenters] = setupCameras(chrom, opts.numCams, ...
                specs.Resolution, specs.Focal, specs.FocalWide, specs.PrincipalPoint);

            fprintf('  spacing=%.3g m  | N=%6d pts  | ', spacing, specs.NumPoints);

            tStart = tic;
            costUnc = resUncertainty(specs, cameras, CamCenters);
            tUnc = toc(tStart);

            tStart = tic;
            costOcc = dynamicOcclusion(specs, cameras, CamCenters);
            tOcc = toc(tStart);

            cost3 = opts.weightUnc * (costUnc / specs.PreComputed.uncertNorm) + ...
                    opts.weightOcc * (costOcc / specs.PreComputed.occlNorm);
            t3    = tUnc + tOcc;

            triplet  = {costUnc, costOcc, cost3};
            tTriplet = {tUnc,    tOcc,    t3};
            for f = 1:nF
                idx = idx + 1;
                results(idx).config    = configs(c).name;
                results(idx).cf        = cfNames{f};
                results(idx).spacing   = spacing;
                results(idx).cost      = triplet{f};
                results(idx).numPoints = specs.NumPoints;
                results(idx).evalTime  = tTriplet{f};
            end

            fprintf('CF1=%.5f (%.1fs) | CF2=%.5f (%.1fs) | CF3=%.5f\n', ...
                costUnc, tUnc, costOcc, tOcc, cost3);
        end
    end

    sweepElapsed = toc(sweepStart);
    fprintf('\nSweep complete in %.1f s (%.2f min)\n', sweepElapsed, sweepElapsed/60);

    %% Pack sweep struct
    sweep = struct();
    sweep.modeTag        = opts.modeTag;
    sweep.timestamp      = datetime('now');
    sweep.spacings       = spacings;
    sweep.volume         = opts.volume;
    sweep.targetMode     = opts.targetMode;
    sweep.targetType     = opts.targetType;
    sweep.numCams        = opts.numCams;
    sweep.zSpacing       = opts.zSpacing;       % [] for isotropic sweeps
    sweep.configs        = configs;
    sweep.cfNames        = cfNames;
    sweep.cfLabels       = cfLabels;
    sweep.results        = results;
    sweep.weightUnc      = opts.weightUnc;
    sweep.weightOcc      = opts.weightOcc;
    sweep.uncertNorm     = specs.PreComputed.uncertNorm;
    sweep.occlNorm       = specs.PreComputed.occlNorm;
    sweep.totalSweepTime = sweepElapsed;

    %% Save .mat
    sweepDir = fullfile(projectRoot, 'Results', 'Sensitivity', opts.modeTag);
    if ~isfolder(sweepDir), mkdir(sweepDir); end
    ts = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    matFile = fullfile(sweepDir, sprintf('spacing_sweep_%s_%s.mat', opts.modeTag, ts));
    save(matFile, 'sweep');
    fprintf('Sweep saved to: %s\n', matFile);

    %% Optional tolerance recommendations
    recs = [];
    if ~isempty(opts.tolerancePct)
        recs = recommendSpacing(sweep, opts.tolerancePct);
    end

    %% Plot
    plotSpacingSensitivity(sweep, projectRoot, ts, opts.tolerancePct, recs);
end
