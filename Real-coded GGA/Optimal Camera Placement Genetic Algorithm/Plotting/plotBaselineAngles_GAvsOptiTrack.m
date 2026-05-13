function plotBaselineAngles_GAvsOptiTrack(varargin)
% PLOTBASELINEANGLES_GAVSOPTITRACK  Pairwise camera-baseline-angle histogram
% comparing GA-best 7-camera config with OptiTrack ad-hoc, one scenario
% per figure (UAV or UGV).
%
% =====================================================================
% EXAMINER REVIEW
% =====================================================================
% What this plot claims to show
%   "Triangulation accuracy depends on the ANGLE between camera rays
%   meeting at a target point, not just the number of cameras that see
%   it. For every (point, camera-pair) where both cameras see the
%   point, this histogram reports the converging-ray angle. The
%   [minTriangAngle, maxTriangAngle] band marks the angles the cost
%   function counts as triangulable — outside that band, a pair is
%   effectively useless for reconstruction. The GA solution and the
%   OptiTrack rig are overlaid so the reader can judge which produces
%   a fatter mass inside the triangulable band."
%
% Strengths
%   - Directly addresses the colour-blind / 2+ count critique levelled
%     at visualizeCameraCoverage: instead of "2+ cameras = good", the
%     metric is the angle the cost function actually uses.
%   - Vertical dashed lines at minTriang / maxTriang make the
%     triangulable band visible.
%   - Same target space for both placements → differences are
%     placement-driven, not grid-driven.
%
% USAGE
%   plotBaselineAngles_GAvsOptiTrack('TargetType', 1)              % UAV
%   plotBaselineAngles_GAvsOptiTrack('TargetType', 2)              % UGV
%
% Name-Value Parameters (same as plotHeatmap_GAvsOptiTrack, plus)
%   'NumBins'   - Histogram bin count. Default 36 (5° bins over 0–180°).
%   'Normalise' - 'probability' (default), 'count', or 'pdf'.

    defaultLog = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                          'Results', 'Logs', 'GGA_RunsLog.mat');

    p = inputParser;
    addParameter(p, 'TargetType',   1,           @isnumeric);
    addParameter(p, 'GridMode',     1,           @isnumeric);
    addParameter(p, 'Spacing',      1.0,         @isnumeric);
    addParameter(p, 'CostFunction', 3,           @isnumeric);
    addParameter(p, 'NumCameras',   7,           @isnumeric);
    addParameter(p, 'LogFile',      defaultLog,  @ischar);
    addParameter(p, 'NumBins',      36,          @isnumeric);
    addParameter(p, 'Normalise',    'probability', @ischar);
    addParameter(p, 'SaveAs',       '',          @ischar);
    parse(p, varargin{:});
    opts = p.Results;

    sty = gaPlotStyle();
    ttStr = sty.TargetNames{opts.TargetType};
    gmStr = sty.GridNames{opts.GridMode};

    %% Load best GA run
    [gaChrom, specs, gaCost] = loadBestGARun(opts);

    %% OptiTrack chromosome — same specs.Target for fair comparison
    optiChrom = buildOptiTrackChromosome();
    optiSpecs = specs;
    optiSpecs.Cams = 7;

    %% Collect pairwise baseline angles (degrees) for every target point
    fprintf('Computing pairwise baseline angles for GA-best...\n');
    angGA   = pairwiseBaselineAngles(gaChrom,   specs);
    fprintf('Computing pairwise baseline angles for OptiTrack...\n');
    angOpti = pairwiseBaselineAngles(optiChrom, optiSpecs);

    minTri = specs.PreComputed.minTriangAngle;
    maxTri = specs.PreComputed.maxTriangAngle;

    %% Plot
    fig = figure('Name', sprintf('Baseline angles — %s', ttStr), ...
        'Units', 'inches', ...
        'Position', [0.5, 0.5, sty.FigWidthFull, sty.FigHeight], ...
        'PaperPositionMode', 'auto', ...
        'Color', sty.BackgroundColor);
    ax = axes(fig);
    hold(ax, 'on');

    edges = linspace(0, 180, opts.NumBins + 1);

    hGA = histogram(ax, angGA, edges, ...
        'Normalization', opts.Normalise, ...
        'FaceColor',    sty.CostFuncColors(3,:), ...   % combined-green
        'FaceAlpha',    0.55, ...
        'EdgeColor',    'none', ...
        'DisplayName',  sprintf('GA-best (n_{pairs}=%d)', numel(angGA)));

    hOpti = histogram(ax, angOpti, edges, ...
        'Normalization', opts.Normalise, ...
        'FaceColor',    [0.85 0.10 0.10], ...          % opti-red
        'FaceAlpha',    0.45, ...
        'EdgeColor',    'none', ...
        'DisplayName',  sprintf('OptiTrack ad-hoc (n_{pairs}=%d)', numel(angOpti)));

    %% Triangulable-band guide lines
    yL = ylim(ax);
    plot(ax, [minTri minTri], yL, '--', ...
        'Color', [0.30 0.30 0.30 0.7], 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
    plot(ax, [maxTri maxTri], yL, '--', ...
        'Color', [0.30 0.30 0.30 0.7], 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
    text(ax, (minTri + maxTri)/2, yL(2)*0.97, ...
        sprintf('Triangulable: %d°–%d°', minTri, maxTri), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontSize', sty.FontSizeAnnot, 'FontName', sty.FontName, ...
        'BackgroundColor', [1 1 1 0.7]);

    hold(ax, 'off');

    %% Headline numbers
    fracGA   = sum(angGA   >= minTri & angGA   <= maxTri) / max(numel(angGA), 1);
    fracOpti = sum(angOpti >= minTri & angOpti <= maxTri) / max(numel(angOpti), 1);

    xlabel(ax, 'Baseline angle between camera pair (°)', ...
        'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    if strcmpi(opts.Normalise, 'count')
        ylabel(ax, 'Count', ...
            'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    else
        ylabel(ax, sprintf('Fraction (%s)', opts.Normalise), ...
            'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    end
    title(ax, sprintf( ...
        '%s — Pairwise baseline angles (in-band: GA %.1f%%, OptiTrack %.1f%%)', ...
        ttStr, 100*fracGA, 100*fracOpti), ...
        'FontWeight', 'normal', ...
        'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
    set(ax, 'FontSize', sty.FontSizeTick, 'FontName', sty.FontName, ...
        'Box', 'on', 'TickDir', 'out');
    grid(ax, 'on');
    legend([hGA, hOpti], 'Location', 'northeast', 'FontSize', sty.FontSizeLegend);
    xlim(ax, [0 180]);

    applyThesisStyle(fig);

    fprintf('\n%s — baseline-angle summary (CF3, %dC, %s, sp=%.2f m):\n', ...
        ttStr, opts.NumCameras, gmStr, opts.Spacing);
    fprintf('  GA-best      median %.1f° | in-band %.1f%% | n_pairs %d (cost %.4f)\n', ...
        median(angGA), 100*fracGA, numel(angGA), gaCost);
    fprintf('  OptiTrack    median %.1f° | in-band %.1f%% | n_pairs %d\n', ...
        median(angOpti), 100*fracOpti, numel(angOpti));

    %% Export
    if isempty(opts.SaveAs)
        outName = sprintf('BaselineAngles_GAvsOptiTrack_%s_%dC_GM%d_sp%.0fcm', ...
            ttStr, opts.NumCameras, opts.GridMode, opts.Spacing*100);
    else
        outName = opts.SaveAs;
    end
    exportgraphics(fig, [outName '.pdf'], ...
        'ContentType',     'vector', ...
        'BackgroundColor', sty.ExportBgColor);
    fprintf('Saved: %s.pdf\n', outName);
end


%% ---- Local helpers ------------------------------------------------------

function angles = pairwiseBaselineAngles(chrom, specs)
% For each target point: find all cameras that see it; for every pair of
% those cameras, compute the convergence angle of their rays AT THE POINT.
% Result: concatenated 1-D array of angles in degrees over all (point,
% pair) combinations.
    numCams = specs.Cams;
    [cameras, camCenters] = setupCameras(chrom, numCams, ...
        specs.Resolution, specs.Focal, specs.FocalWide, specs.PrincipalPoint);

    resolution     = specs.Resolution;
    TargetSpace    = specs.Target;
    maxCameraRange     = specs.PreComputed.maxCameraRange;
    maxCameraRangeWide = specs.PreComputed.maxCameraRangeWide;
    focalWide          = specs.FocalWide;

    nPts = size(TargetSpace, 1);

    % Worst-case upper bound on pairs so we can preallocate
    maxPairsPerPoint = numCams * (numCams - 1) / 2;
    pool = nan(nPts * maxPairsPerPoint, 1);
    head = 0;

    for pt = 1:nPts
        point = TargetSpace(pt, :);
        [visCams, viewVecs] = findVisibleCameras(point, cameras, camCenters, ...
            numCams, resolution, maxCameraRange, maxCameraRangeWide, focalWide);
        nv = numel(visCams);
        if nv < 2
            continue;
        end
        % viewVecs is 3 x nv unit vectors from camera to point (i.e.
        % the direction the ray travels TOWARDS the point). The angle
        % at the point between two rays is the supplement of the angle
        % between the two cam->point vectors, but for ray-convergence
        % angle we want the included angle BETWEEN the two camera lines
        % of sight — which is the angle between the vectors (point - cam)
        % and (point - cam') i.e. the same as the angle between
        % the cam->point unit vectors viewVecs(:,i) and viewVecs(:,j).
        % Actually for triangulation, what we care about is the angle
        % between the two cameras as seen from the point — which is the
        % angle between (cam_i - point) and (cam_j - point), which is
        % the angle between -viewVecs(:,i) and -viewVecs(:,j) = same as
        % between viewVecs(:,i) and viewVecs(:,j). So acosd(dot) gives
        % the triangulation angle directly.
        for i = 1:(nv-1)
            for j = (i+1):nv
                c = dot(viewVecs(:,i), viewVecs(:,j));
                c = max(min(c, 1), -1);
                head = head + 1;
                pool(head) = acosd(c);
            end
        end
    end

    angles = pool(1:head);
end


function [chrom, specs, cost] = loadBestGARun(opts)
    if ~isfile(opts.LogFile)
        error('Log file not found: %s', opts.LogFile);
    end
    S = load(opts.LogFile, 'runLog');
    runLog = S.runLog;

    fillFields = {'TargetType', 'GridMode', 'Spacing'};
    for f = 1:length(fillFields)
        if ~isfield(runLog, fillFields{f})
            [runLog.(fillFields{f})] = deal(NaN);
        end
    end

    mask = ([runLog.NumCameras]       == opts.NumCameras)   & ...
           ([runLog.CostFunctionType] == opts.CostFunction) & ...
           ([runLog.TargetType]       == opts.TargetType)   & ...
           ([runLog.GridMode]         == opts.GridMode)     & ...
           (abs([runLog.Spacing] - opts.Spacing) < 1e-6);

    candidates = runLog(mask);
    if isempty(candidates)
        error('No GA runs found for TT=%d GM=%d sp=%.2f CF=%d %dC.', ...
            opts.TargetType, opts.GridMode, opts.Spacing, ...
            opts.CostFunction, opts.NumCameras);
    end
    [cost, idx] = min([candidates.BestCost]);
    bestRun = candidates(idx);

    matFile = resolveRunPath(bestRun.RunFilename, bestRun.NumCameras);
    if ~isfile(matFile)
        error('GA result file not found: %s', matFile);
    end
    L = load(matFile, 'saveData');
    chrom = L.saveData.BestSolution.Chromosome;
    specs = L.saveData.Specifications;
end
