function out = evaluateOptiTrackCost(varargin)
% EVALUATEOPTITRACKCOST  Evaluate the OptiTrack 7-cam config under given conditions.
%
%   out = evaluateOptiTrackCost('Name', Value, ...)
%
%   Returns a struct with fields:
%     CF1   - raw resolution-uncertainty cost
%     CF2   - raw dynamic-occlusion cost
%     CF3   - weight-combined CF3 cost, built via cf3Terms: each component
%             min-max scaled as (raw - utopia)/norm then weighted, at the
%             same target space (utopia = 0 when no calibrated normTable)
%     WeightedUnc - the J_uncertainty term as it appears inside CF3, i.e.
%                   WeightedUnc + WeightedOcc == CF3
%     WeightedOcc - the J_occlusion term
%     UncertNorm, OcclNorm - the normalisation constants used above
%     numCams, targetType, gridMode, spacing
%
%   The function reuses the same pipeline as runCameraOptimiser /
%   spacingSensitivityCore: same hardware specs, same target-space
%   builder, same cost-param precomputations. The output is therefore
%   directly comparable to a GA's BestCost at the same conditions and
%   safe to overlay on box plots / convergence curves as an "ad-hoc
%   baseline".
%
%   The result is cached per (TargetType, GridMode, Spacing, weights)
%   tuple in a function-scope persistent map, so repeated calls in a
%   plotting loop are cheap.
%
%   Name-Value Parameters:
%     'TargetType'   - 1 (UAV) or 2 (UGV). Default: 1
%     'GridMode'     - 1 (Uniform) or 2 (Normal). Default: 1
%     'Spacing'      - x-y grid spacing in m. Default: 1.0
%     'UGVmaxHeight' - max height of UGV slab. Default: 0.5
%     'UGVzSpacing'  - z step on UGV slab. Default: 0.25
%     'WeightUnc'    - CF3 uncertainty weight. Default: 0.5
%     'WeightOcc'    - CF3 occlusion weight. Default: 0.5
%     'Verbose'      - true to print eval line (default: false)

    p = inputParser;
    addParameter(p, 'TargetType',   1,    @isnumeric);
    addParameter(p, 'GridMode',     1,    @isnumeric);
    addParameter(p, 'Spacing',      1.0,  @isnumeric);
    addParameter(p, 'UGVmaxHeight', 0.5,  @isnumeric);
    addParameter(p, 'UGVzSpacing',  0.25, @isnumeric);
    addParameter(p, 'WeightUnc',    0.5,  @isnumeric);
    addParameter(p, 'WeightOcc',    0.5,  @isnumeric);
    addParameter(p, 'Verbose',      false,@islogical);
    addParameter(p, 'Preset',       'optitrack_lab', @ischar);   % see runConfig
    parse(p, varargin{:});
    opts = p.Results;

    %% Cache lookup ---------------------------------------------------
    persistent cache
    if isempty(cache)
        cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    end
    cacheKey = sprintf('%s_TT%d_GM%d_sp%.4f_wu%.3f_wo%.3f', opts.Preset, ...
        opts.TargetType, opts.GridMode, opts.Spacing, ...
        opts.WeightUnc, opts.WeightOcc);
    if isKey(cache, cacheKey)
        out = cache(cacheKey);
        return;
    end

    %% Project path ---------------------------------------------------
    addProjectPaths();

    numCams = 7;

    %% Specs: same builder as the GA runs
    run = runConfig(opts.Preset, 'UGV_MaxHeight', opts.UGVmaxHeight, ...
        'UGV_ZSpacing', opts.UGVzSpacing, 'Weights', [opts.WeightUnc, opts.WeightOcc]);
    specs = buildRunSpecs(run, numCams, 3, opts.TargetType, opts.GridMode, opts.Spacing);

    %% Chromosome + cameras ------------------------------------------
    chrom = buildOptiTrackChromosome();
    [cameras, CamCenters] = setupCameras(chrom, numCams, ...
        specs.Resolution, specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);

    %% Evaluate -------------------------------------------------------
    tStart = tic;
    costUnc = resUncertainty(specs, cameras, CamCenters);
    costOcc = dynamicOcclusion(specs, cameras, CamCenters);
    [cost3, weightedUnc, weightedOcc] = cf3Terms(costUnc, costOcc, specs);
    tEval = toc(tStart);

    out.CF1        = costUnc;
    out.CF2        = costOcc;
    out.CF3        = cost3;
    out.WeightedUnc = weightedUnc;
    out.WeightedOcc = weightedOcc;
    out.UncertNorm  = specs.PreComputed.uncertNorm;
    out.OcclNorm    = specs.PreComputed.occlNorm;
    out.numCams    = numCams;
    out.targetType = opts.TargetType;
    out.gridMode   = opts.GridMode;
    out.spacing    = opts.Spacing;
    out.numPoints  = specs.NumPoints;
    out.evalTime   = tEval;

    if opts.Verbose
        fprintf(['evaluateOptiTrackCost: TT=%d GM=%d sp=%.2f m  -> ' ...
                 'CF1=%.5f CF2=%.5f CF3=%.5f  (N=%d pts, %.2fs)\n'], ...
                 opts.TargetType, opts.GridMode, opts.Spacing, ...
                 costUnc, costOcc, cost3, specs.NumPoints, tEval);
    end

    cache(cacheKey) = out;
end
