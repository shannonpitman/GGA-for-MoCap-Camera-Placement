function stats = sanityCheckCoverage(varargin)
% SANITYCHECKCOVERAGE  Per-camera and per-point coverage diagnostic.
%
% After applying the lens-range / uv-indexing / isInsidePlanes fixes, run
% this on a saved GA result to confirm the corrected coverage numbers
% match what you see when the camera FOVs are plotted.
%
% USAGE:
%   sanityCheckCoverage()                              % Best 7C CF3 UAV Uniform
%   sanityCheckCoverage('NumCameras', 8)               % Best 8C CF3 UAV Uniform
%   sanityCheckCoverage('TargetType', 2)               % Best UGV result
%   sanityCheckCoverage('GridMode', 2)                 % Best Normal grid result
%   sanityCheckCoverage('RunFile', '7Cams_Run_X.mat')  % Specific result file
%   sanityCheckCoverage('Plot', false)                 % Skip figure
%
% Returns a struct of summary stats so you can compare across runs.

    addProjectPaths();

    %% Parse inputs
    defaultLog = fullfile(fileparts(mfilename('fullpath')), ...
                          'Results', 'Logs', 'GGA_RunsLog.mat');

    p = inputParser;
    addParameter(p, 'LogFile',    defaultLog, @ischar);
    addParameter(p, 'NumCameras', 7,          @isnumeric);
    addParameter(p, 'TargetType', 1,          @isnumeric);
    addParameter(p, 'GridMode',   1,          @isnumeric);
    addParameter(p, 'RunFile',    '',         @ischar);
    addParameter(p, 'Plot',       true,       @islogical);
    parse(p, varargin{:});
    opts = p.Results;

    %% Load saved data
    if ~isempty(opts.RunFile)
        runFile = resolveRunPath(opts.RunFile);
        tmp = load(runFile, 'saveData');
        sd  = tmp.saveData;
        fprintf('Loaded: %s\n', runFile);
    else
        load(opts.LogFile, 'runLog');
        if ~isfield(runLog, 'TargetType'), [runLog.TargetType] = deal(NaN); end
        if ~isfield(runLog, 'GridMode'),   [runLog.GridMode]   = deal(NaN); end

        mask = ([runLog.NumCameras]       == opts.NumCameras) & ...
               ([runLog.CostFunctionType] == 3) & ...
               ([runLog.TargetType]       == opts.TargetType) & ...
               ([runLog.GridMode]         == opts.GridMode);
        filtered = runLog(mask);
        if isempty(filtered)
            error('No matching CF3 runs found for %dC TT%d GM%d.', ...
                opts.NumCameras, opts.TargetType, opts.GridMode);
        end
        [~, bestIdx] = min([filtered.BestCost]);
        bestRun  = filtered(bestIdx);
        bestPath = resolveRunPath(bestRun.RunFilename, bestRun.NumCameras);
        tmp = load(bestPath, 'saveData');
        sd  = tmp.saveData;
        fprintf('Loaded best CF3: %s\n', bestPath);
    end

    %% Build cameras with the SAME specs the cost function uses
    specs       = sd.Specifications;
    specs       = backfillLegacySpecs(specs);   % top-up legacy saves
    chrom       = sd.BestSolution.Chromosome;
    numCams     = specs.Cams;
    targetSpace = specs.Target;
    numPoints   = size(targetSpace, 1);

    [cameras, camCenters] = setupCameras(chrom, numCams, specs.Resolution, ...
        specs.Focal, specs.FocalWide, specs.PrincipalPoint);

    maxRange     = specs.PreComputed.maxCameraRange;
    maxRangeWide = specs.PreComputed.maxCameraRangeWide;
    focalWide    = specs.FocalWide;

    %% --- Per-camera summary ---
    fprintf('\n');
    fprintf('==========================================================\n');
    fprintf('  PER-CAMERA SUMMARY\n');
    fprintf('==========================================================\n');
    fprintf('  %-4s %-7s %-10s %-10s %-12s\n', ...
        'Cam', 'Lens', 'f [mm]', 'Range [m]', '#pts seen');

    pointsSeenByCam = zeros(numCams, 1);
    for i = 1:numCams
        if cameras{i}.f == focalWide
            lensStr  = 'wide';
            rng_i    = maxRangeWide;
        else
            lensStr  = 'narrow';
            rng_i    = maxRange;
        end

        cnt = 0;
        for q = 1:numPoints
            point = targetSpace(q,:);
            uv = cameras{i}.project(point);
            if uv(1) >= 1 && uv(1) <= specs.Resolution(1) && ...
               uv(2) >= 1 && uv(2) <= specs.Resolution(2)
                d = norm(point - camCenters(:,i)');
                if d > 0 && d <= rng_i
                    cnt = cnt + 1;
                end
            end
        end
        pointsSeenByCam(i) = cnt;
        fprintf('  %-4d %-7s %-10.2f %-10.2f %-12d\n', ...
            i, lensStr, cameras{i}.f * 1000, rng_i, cnt);
    end

    %% --- Per-point coverage via the (now corrected) findVisibleCameras ---
    coverage = zeros(numPoints, 1);
    for q = 1:numPoints
        [vis, ~] = findVisibleCameras(targetSpace(q,:), cameras, camCenters, ...
            numCams, specs.Resolution, maxRange, maxRangeWide, focalWide);
        coverage(q) = numel(vis);
    end

    stats.minCov     = min(coverage);
    stats.maxCov     = max(coverage);
    stats.avgCov     = mean(coverage);
    stats.medianCov  = median(coverage);
    stats.zeroPts    = sum(coverage == 0);
    stats.onePt      = sum(coverage == 1);
    stats.twoPlusPts = sum(coverage >= 2);
    stats.numPoints  = numPoints;
    stats.coverage   = coverage;

    %% --- Coverage summary ---
    fprintf('\n');
    fprintf('==========================================================\n');
    fprintf('  COVERAGE SUMMARY (post-fix)\n');
    fprintf('==========================================================\n');
    fprintf('  Target points: %d\n', numPoints);
    fprintf('  Min cameras per point:    %d\n', stats.minCov);
    fprintf('  Max cameras per point:    %d\n', stats.maxCov);
    fprintf('  Avg cameras per point:    %.2f\n', stats.avgCov);
    fprintf('  Median cameras per point: %.1f\n', stats.medianCov);
    fprintf('  Points seen by 0 cams:  %d (%.1f%%)\n', ...
        stats.zeroPts, 100*stats.zeroPts/numPoints);
    fprintf('  Points seen by 1 cam:   %d (%.1f%%)\n', ...
        stats.onePt,   100*stats.onePt/numPoints);
    fprintf('  Points seen by 2+ cams: %d (%.1f%%)\n', ...
        stats.twoPlusPts, 100*stats.twoPlusPts/numPoints);
    fprintf('==========================================================\n\n');

    %% --- Optional figure: FOV frustums + target points coloured by coverage ---
    if opts.Plot
        figure('Name', 'Sanity check: FOV frustums + target coverage', ...
            'Position', [100 100 1000 700]);
        hold on;

        for i = 1:numCams
            if cameras{i}.f == focalWide
                rng_i = maxRangeWide;
                col_i = [0.85 0.40 0.10];
            else
                rng_i = maxRange;
                col_i = [0.10 0.45 0.75];
            end
            plotCameraFOV(cameras{i}, camCenters(:,i), rng_i, ...
                'Color', col_i, 'Label', sprintf('cam%d', i));
        end

        % Scatter target points coloured by coverage; zero-coverage in red
        cov = coverage;
        scatter3(targetSpace(:,1), targetSpace(:,2), targetSpace(:,3), ...
            18, cov, 'filled', 'MarkerFaceAlpha', 0.65);
        % Highlight zero-coverage points
        if any(cov == 0)
            zPts = targetSpace(cov == 0, :);
            scatter3(zPts(:,1), zPts(:,2), zPts(:,3), 60, ...
                'MarkerEdgeColor', [0.85 0.10 0.10], ...
                'MarkerFaceColor', [0.95 0.20 0.20], 'LineWidth', 1.0);
        end

        cb = colorbar; cb.Label.String = 'Cameras seeing point';
        colormap(parula);
        clim([0 max(numCams, 1)]);

        axis equal; grid on; view(45, 25);
        xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
        title(sprintf('%dC sanity check  |  min cov = %d, 0-cov pts = %d', ...
            numCams, stats.minCov, stats.zeroPts));
        hold off;
    end
end
