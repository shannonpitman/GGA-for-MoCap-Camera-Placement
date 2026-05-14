function saveOptiTrackAsRun(varargin)
%SAVEOPTITRACKASRUN  Persist the OptiTrack ad-hoc 7-cam config to Results/.
%
%   saveOptiTrackAsRun() builds the OptiTrack chromosome via
%   buildOptiTrackChromosome, evaluates it at the standard 0.5 m UAV
%   uniform grid, and writes a per-run .mat plus a companion .txt file
%   into Results/7Cams/, using the saveResults format so that
%   analyseConfiguration('RunFile', ...) can load it directly.
%
%   Both filenames and the .txt body explicitly mark this as the
%   "OptiTrackAdHoc" configuration so it isn't mistaken for a GA result.
%   The master GGA_RunsLog is intentionally NOT touched — this entry is
%   a measured baseline, not an optimised run.
%
%   USAGE:
%       saveOptiTrackAsRun()                       % defaults below
%       saveOptiTrackAsRun('Spacing', 0.25)        % finer evaluation grid
%       saveOptiTrackAsRun('TargetType', 2)        % UGV slab (z<=0.5 m)
%
%   PARAMETERS (Name-Value):
%       Spacing      grid spacing for the cost evaluation (default 0.5 m)
%       Volume       3x2 workspace volume (default UAV envelope)
%       TargetType   1 = UAV (full volume), 2 = UGV (floor slab)
%       TargetMode   1 = uniform grid, 2 = inverse-CDF concentrated
%       WeightUnc    CF3 weight on resolution uncertainty (default 0.5)
%       WeightOcc    CF3 weight on dynamic occlusion       (default 0.5)
%
%   Run once per change to the OptiTrack rig measurements, then re-use the
%   resulting .mat from analyseConfiguration / sanity scripts.

    addProjectPaths();
    projectRoot = fileparts(mfilename('fullpath'));

    %% Parse inputs
    p = inputParser;
    addParameter(p, 'Spacing',    0.5,                @(x) isnumeric(x) && x > 0);
    addParameter(p, 'Volume',     [-4 4; -4 4; 0 4],  @isnumeric);
    addParameter(p, 'TargetType', 1,                  @(x) any(x == [1 2]));
    addParameter(p, 'TargetMode', 1,                  @(x) any(x == [1 2]));
    addParameter(p, 'WeightUnc',  0.5,                @isnumeric);
    addParameter(p, 'WeightOcc',  0.5,                @isnumeric);
    parse(p, varargin{:});
    opt = p.Results;

    if opt.TargetType == 2
        opt.Volume(3,:) = [0, 0.5];           % UGV slab
        if opt.Spacing > 0.25
            warning('saveOptiTrackAsRun:UGVspacing', ...
                'Spacing %.2f m is too coarse for UGV slab; using 0.25 m.', opt.Spacing);
            opt.Spacing = 0.25;
        end
    end

    %% Build chromosome + specs
    numCams = 7;
    chrom = buildOptiTrackChromosome();

    specs = setupHardwareSpecs(numCams);
    specs.WeightUncertainty = opt.WeightUnc;
    specs.WeightOcclusion   = opt.WeightOcc;
    specs.TargetType        = opt.TargetType;
    specs.TargetMode        = opt.TargetMode;
    specs.Target            = generateTargetSpace(opt.Volume, opt.TargetMode, opt.Spacing);
    specs.NumPoints         = size(specs.Target, 1);
    specs.spacing           = opt.Spacing;
    specs                   = setupCostParams(specs);
    specs.warmStart         = false;
    specs.warmChromosomes   = [];

    %% Evaluate cost components
    [cameras, CamCenters] = setupCameras(chrom, numCams, ...
        specs.Resolution, specs.Focal, specs.FocalWide, specs.PrincipalPoint, specs.PixelSize);

    fprintf('\nEvaluating OptiTrack ad-hoc rig (N=%d points)...\n', specs.NumPoints);
    tStart = tic;
    costUnc = resUncertainty(specs, cameras, CamCenters);
    costOcc = dynamicOcclusion(specs, cameras, CamCenters);
    cost3   = opt.WeightUnc * (costUnc / specs.PreComputed.uncertNorm) + ...
              opt.WeightOcc * (costOcc / specs.PreComputed.occlNorm);
    elapsed = toc(tStart);

    fprintf('  raw uncertainty = %.5f\n', costUnc);
    fprintf('  raw occlusion   = %.5f\n', costOcc);
    fprintf('  combined CF3    = %.5f   (eval %.1f s)\n', cost3, elapsed);

    %% Coverage stats (inlined — no plotting)
    coverageStats = computeCoverageStats(specs, cameras, CamCenters);

    %% Build saveData struct (mirrors saveResults.m schema)
    sd = struct();
    sd.BestSolution.Chromosome = chrom;
    sd.BestSolution.Cost       = cost3;
    sd.BestCost                = cost3;
    sd.CameraConfiguration     = reshape(chrom, 6, numCams)';
    sd.Specifications          = specs;
    sd.GAParams                = struct();          % no GA was used
    sd.ConvergenceHistory      = [];                % no convergence
    sd.AvgCostHistory          = [];
    sd.TopTenAvgCostHistory    = [];
    sd.ElapsedTime             = elapsed;
    sd.Timestamp               = datetime('now');
    sd.CoverageStats           = coverageStats;

    for i = 1:numCams
        i0 = (i-1)*6;
        sd.Cameras(i).Position           = chrom(i0+1:i0+3);
        sd.Cameras(i).Orientation        = chrom(i0+4:i0+6);
        sd.Cameras(i).OrientationDegrees = rad2deg(chrom(i0+4:i0+6));
    end

    % Tag this as the ad-hoc reference (not a GA run)
    sd.IsOptiTrackAdHoc        = true;
    sd.Source                  = 'OptiTrack Motive measured rig (Y-up -> Z-up transformed)';
    sd.RawCosts.Uncertainty    = costUnc;
    sd.RawCosts.Occlusion      = costOcc;
    sd.RawCosts.CombinedCF3    = cost3;

    %% Write .mat
    saveData = sd; %#ok<NASGU>  (kept for analyseConfiguration compatibility)
    runDir = fullfile(projectRoot, 'Results', sprintf('%dCams', numCams));
    if ~isfolder(runDir), mkdir(runDir); end

    ts = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    matFilename = sprintf('%dCams_OptiTrackAdHoc_%s.mat', numCams, ts);
    matPath = fullfile(runDir, matFilename);
    save(matPath, 'saveData');
    fprintf('\nSaved chromosome+stats to: %s\n', matPath);

    %% Write .txt summary (clearly labelled)
    txtFilename = sprintf('%dCams_OptiTrackAdHoc_%s.txt', numCams, ts);
    txtPath = fullfile(runDir, txtFilename);
    fid = fopen(txtPath, 'w');
    targetTypeNames = {'UAV (Full Volume)', 'UGV (Floor Slab)'};
    gridModeNames   = {'Uniform Grid', 'Normal Grid'};

    fprintf(fid, 'OPTITRACK AD-HOC CAMERA CONFIGURATION\n');
    fprintf(fid, '=====================================\n');
    fprintf(fid, '** This is NOT a GA-optimised result. **\n');
    fprintf(fid, 'It is the as-measured OptiTrack Motive rig converted into\n');
    fprintf(fid, 'the GA chromosome representation, evaluated through the\n');
    fprintf(fid, 'GA cost function for direct comparison against optimised runs.\n\n');
    fprintf(fid, 'Timestamp: %s\n', char(sd.Timestamp));
    fprintf(fid, 'Number of Cameras: %d\n', numCams);
    fprintf(fid, 'Source: %s\n', sd.Source);
    fprintf(fid, 'Target Type: %s\n', targetTypeNames{opt.TargetType});
    fprintf(fid, 'Grid Mode:   %s\n', gridModeNames{opt.TargetMode});
    fprintf(fid, 'Grid Spacing: %.3f m   (N = %d target points)\n', ...
        opt.Spacing, specs.NumPoints);
    fprintf(fid, 'Hardware: %dx%d res, focal %.1f/%.1f mm (narrow/wide), range %.1f/%.1f m\n', ...
        specs.Resolution(1), specs.Resolution(2), ...
        specs.Focal*1e3, specs.FocalWide*1e3, specs.Range, specs.RangeWide);
    fprintf(fid, '\nEvaluated Costs (alternating narrow/wide focal per setupCameras):\n');
    fprintf(fid, '  Resolution Uncertainty (raw mean): %.6f\n', costUnc);
    fprintf(fid, '  Dynamic Occlusion       (raw mean): %.6f\n', costOcc);
    fprintf(fid, '  Combined CF3 (normalised):          %.6f\n', cost3);
    fprintf(fid, '  Cost Function Weights: %.2f resolution + %.2f occlusion\n', ...
        opt.WeightUnc, opt.WeightOcc);
    fprintf(fid, '  Evaluation time: %.2f s\n', elapsed);

    fprintf(fid, '\nCamera Configurations (after Y-up -> Z-up transform):\n');
    fprintf(fid, '----------------------------------------------------\n');
    for i = 1:numCams
        fprintf(fid, '\nCamera %d:\n', i);
        fprintf(fid, '  Position (m):       [%.3f, %.3f, %.3f]\n', sd.Cameras(i).Position);
        fprintf(fid, '  Orientation (rad):  [%.3f, %.3f, %.3f]\n', sd.Cameras(i).Orientation);
        fprintf(fid, '  Orientation (deg):  [%.1f, %.1f, %.1f]\n', sd.Cameras(i).OrientationDegrees);
    end

    fprintf(fid, '\nCamera Coverage Statistics:\n');
    fprintf(fid, '---------------------------\n');
    fprintf(fid, 'Total target points: %d\n', coverageStats.numPoints);
    fprintf(fid, 'Points with 0 cameras: %d (%.1f%%)\n', ...
        coverageStats.zeroCameras, coverageStats.zeroCamerasPercent);
    fprintf(fid, 'Points with 1 camera:  %d (%.1f%%)\n', ...
        coverageStats.oneCamera, coverageStats.oneCameraPercent);
    fprintf(fid, 'Points with 2+ cameras: %d (%.1f%%)\n', ...
        coverageStats.twoPlusCameras, coverageStats.twoPlusCamerasPercent);
    fprintf(fid, 'Average coverage: %.2f cameras per point\n', coverageStats.avgCoverage);
    fprintf(fid, 'Max coverage: %d | Min: %d | Median: %.2f\n', ...
        coverageStats.maxCoverage, coverageStats.minCoverage, coverageStats.medianCoverage);

    fprintf(fid, '\nNote: The master log (GGA_RunsLog.mat) is intentionally NOT modified\n');
    fprintf(fid, 'by this script; this file lives only as an annotated reference.\n');
    fclose(fid);

    fprintf('Wrote summary to:        %s\n', txtPath);
    fprintf('\nLoad it later via: analyseConfiguration(''RunFile'', ''%s'')\n', matFilename);
end


function stats = computeCoverageStats(specs, cameras, CamCenters)
%COMPUTECOVERAGESTATS  Coverage statistics without the plotting overhead.
    numCams      = specs.Cams;
    resolution   = specs.Resolution;
    target       = specs.Target;
    maxRange     = specs.PreComputed.maxCameraRange;
    maxRangeWide = specs.PreComputed.maxCameraRangeWide;
    focalWide    = specs.FocalWide;

    nP = size(target, 1);
    cov = zeros(nP, 1);
    parfor q = 1:nP
        [vis, ~] = findVisibleCameras(target(q,:), cameras, CamCenters, ...
            numCams, resolution, maxRange, maxRangeWide, focalWide);
        cov(q) = numel(vis);
    end

    stats.numPoints           = nP;
    stats.zeroCameras         = sum(cov == 0);
    stats.oneCamera           = sum(cov == 1);
    stats.twoPlusCameras      = sum(cov >= 2);
    stats.zeroCamerasPercent  = 100 * stats.zeroCameras / nP;
    stats.oneCameraPercent    = 100 * stats.oneCamera / nP;
    stats.twoPlusCamerasPercent = 100 * stats.twoPlusCameras / nP;
    stats.avgCoverage         = mean(cov);
    stats.maxCoverage         = max(cov);
    stats.minCoverage         = min(cov);
    stats.medianCoverage      = median(cov);
end
