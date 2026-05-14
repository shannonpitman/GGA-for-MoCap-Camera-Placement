function plotHeatmap_GAvsOptiTrack(varargin)
% PLOTHEATMAP_GAVSOPTITRACK  Coverage heat-map: GA best vs OptiTrack ad-hoc.
%
% =====================================================================
% EXAMINER REVIEW
% =====================================================================
% What this plot claims to show
%   "For a single scenario (UAV or UGV) at 7 cameras on a uniform grid,
%   the GA's best configuration achieves a more even spatial coverage
%   of the target volume than the ad-hoc OptiTrack rig measured in the
%   lab. The reader can see WHERE the ad-hoc rig under-covers (zero or
%   one camera per point) and how the GA redistributes cameras to
%   reduce that under-coverage."
%
% Layout (per scenario)
%   2x2 figure:
%     Top-Left : GA-best 3D scatter         | Top-Right : OptiTrack 3D
%     Bot-Left : GA-best XY worst-case      | Bot-Right : OptiTrack XY
%   Within a scenario the four panels share a single colour scale of
%   [0 .. globalMax(visibleCams)] so colours are directly comparable
%   left-to-right. Scenarios are NOT scaled against each other (UGV's
%   floor slab and UAV's full volume are different problems).
%
% Strengths over the older Uniform-vs-Normal coverage heatmap
%   - Uses the exact same target space for both configurations, so
%     visibility differences are due solely to camera placement.
%   - Pins to CF3 / 7-cam / Uniform by default (the conditions the GA
%     was trained against the OptiTrack baseline on).
%
% USAGE
%   plotHeatmap_GAvsOptiTrack('TargetType', 1)               % UAV
%   plotHeatmap_GAvsOptiTrack('TargetType', 2)               % UGV
%   plotHeatmap_GAvsOptiTrack('TargetType', 1, 'Spacing', 0.5)
%
% Name-Value Parameters
%   'TargetType'   - 1 (UAV) | 2 (UGV). Default 1.
%   'GridMode'     - 1 (Uniform) | 2 (Normal). Default 1.
%   'Spacing'      - Grid spacing in m. Default 1.0.
%   'CostFunction' - 1, 2, or 3. Default 3 (combined).
%   'NumCameras'   - GA camera count. Default 7. Must be 7 for a fair
%                    comparison with the OptiTrack rig.
%   'LogFile'      - Master log path. Default Results/Logs/GGA_RunsLog.mat.
%   'MarkerSize'   - Scatter marker area. Default 40.
%   'ViewAngle'    - 3D view angle [az el]. Default [45 30].
%   'SaveAs'       - Output filename (no extension). Default auto.

    %% Parse inputs
    defaultLog = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                          'Results', 'Logs', 'GGA_RunsLog.mat');

    p = inputParser;
    addParameter(p, 'TargetType',   1,        @isnumeric);
    addParameter(p, 'GridMode',     1,        @isnumeric);
    addParameter(p, 'Spacing',      1.0,      @isnumeric);
    addParameter(p, 'CostFunction', 3,        @isnumeric);
    addParameter(p, 'NumCameras',   7,        @isnumeric);
    addParameter(p, 'LogFile',      defaultLog, @ischar);
    addParameter(p, 'MarkerSize',   40,       @isnumeric);
    addParameter(p, 'ViewAngle',    [45 30],  @isnumeric);
    addParameter(p, 'SaveAs',       '',       @ischar);
    parse(p, varargin{:});
    opts = p.Results;

    %% Sanity check
    if opts.NumCameras ~= 7
        warning('plotHeatmap_GAvsOptiTrack:NonStandardCams', ...
            ['NumCameras = %d but the OptiTrack rig has 7 cameras. The ' ...
             'visual comparison is no longer apples-to-apples.'], opts.NumCameras);
    end

    sty = gaPlotStyle();
    ttNames = sty.TargetNames;            % {'UAV','UGV'}
    gmNames = sty.GridNames;              % {'Uniform','Normal'}
    ttStr   = ttNames{opts.TargetType};
    gmStr   = gmNames{opts.GridMode};

    %% Locate best GA run for the requested scenario
    if ~isfile(opts.LogFile)
        error('Log file not found: %s', opts.LogFile);
    end
    S = load(opts.LogFile, 'runLog');
    runLog = S.runLog;

    % Ensure expected fields exist (some legacy rows may be missing them)
    requiredFields = {'TargetType', 'GridMode', 'Spacing'};
    for f = 1:length(requiredFields)
        if ~isfield(runLog, requiredFields{f})
            [runLog.(requiredFields{f})] = deal(NaN);
        end
    end

    mask = ([runLog.NumCameras]        == opts.NumCameras)   & ...
           ([runLog.CostFunctionType]  == opts.CostFunction) & ...
           ([runLog.TargetType]        == opts.TargetType)   & ...
           ([runLog.GridMode]          == opts.GridMode)     & ...
           (abs([runLog.Spacing] - opts.Spacing) < 1e-6);

    candidates = runLog(mask);
    if isempty(candidates)
        error(['No GA runs found for TT=%d GM=%d Spacing=%.2f CF=%d ' ...
               'NumCameras=%d. Run batchRunGA first.'], ...
               opts.TargetType, opts.GridMode, opts.Spacing, ...
               opts.CostFunction, opts.NumCameras);
    end
    [gaCost, bestIdx] = min([candidates.BestCost]);
    bestRun = candidates(bestIdx);

    bestMat = resolveRunPath(bestRun.RunFilename, bestRun.NumCameras);
    if ~isfile(bestMat)
        error('GA result file not found: %s', bestMat);
    end
    L = load(bestMat, 'saveData');
    sd = L.saveData;

    specs   = sd.Specifications;
    gaChrom = sd.BestSolution.Chromosome;

    fprintf(['plotHeatmap_GAvsOptiTrack: %s scenario, GM=%s, sp=%.2f m, ' ...
             'CF=%d, %d cams.\n  GA-best run: %s (Cost=%.4f)\n'], ...
             ttStr, gmStr, opts.Spacing, opts.CostFunction, ...
             opts.NumCameras, bestMat, gaCost);

    %% Build OptiTrack chromosome (always 7 cams) and its CF3 cost on the
    %  same problem (GridMode + Spacing + TargetType) so the title is
    %  meaningful.
    optiChrom = buildOptiTrackChromosome();
    optiCost  = evaluateOptiTrackCost( ...
        'TargetType', opts.TargetType, ...
        'GridMode',   opts.GridMode, ...
        'Spacing',    opts.Spacing);
    cfField   = sprintf('CF%d', opts.CostFunction);
    if isfield(optiCost, cfField)
        optiCostVal = optiCost.(cfField);
    else
        optiCostVal = NaN;
    end

    %% Compute per-point camera-visibility for both configurations.
    %  Use the GA run's specs.Target so the comparison is on the same
    %  set of evaluation points.
    [covGA,   numCamsGA]   = perPointVisibility(gaChrom,   specs);
    [covOpti, numCamsOpti] = perPointVisibility(optiChrom, specs);

    if numCamsGA ~= numCamsOpti
        warning(['Camera counts differ: GA=%d, OptiTrack=%d. The colour ' ...
                 'scale will still be shared but the comparison should be ' ...
                 'discussed in text.'], numCamsGA, numCamsOpti);
    end

    %% Figure layout — 2 rows x 2 cols, shared scale within scenario
    figW = sty.FigWidthFull;            % side-by-side cols
    figH = sty.FigHeightTall * 1.1;     % 2 rows

    fig = figure('Name', sprintf('Coverage: GA vs OptiTrack — %s', ttStr), ...
        'Units', 'inches', ...
        'Position', [0.5, 0.5, figW, figH], ...
        'PaperPositionMode', 'auto', ...
        'Color', sty.BackgroundColor);

    globalMaxCams = max([numCamsGA, numCamsOpti, max(covGA), max(covOpti)]);
    cmap = parula(max(globalMaxCams + 1, 2));

    % Pre-pack panels so we can use one loop
    panels(1).chrom = gaChrom;    panels(1).cov = covGA;
    panels(1).name  = sprintf('GA-best (Cost: %.4f)', gaCost);
    panels(2).chrom = optiChrom;  panels(2).cov = covOpti;
    panels(2).name  = sprintf('OptiTrack ad-hoc (Cost: %.4f)', optiCostVal);

    TargetSpace = specs.Target;

    %% Plot top row (3D scatter) + bottom row (XY worst-case)
    for k = 1:2
        cov  = panels(k).cov;
        name = panels(k).name;

        zeroPct    = 100 * sum(cov == 0)  / specs.NumPoints;
        twoPlusPct = 100 * sum(cov >= 2)  / specs.NumPoints;
        avgCov     = mean(cov);
        statLine   = sprintf('Avg %.2f | 0-cam %.1f%% | 2+cam %.1f%%', ...
                              avgCov, zeroPct, twoPlusPct);

        % ---- Top row: 3D scatter ----
        ax3D = subplot(2, 2, k);
        scatter3(ax3D, TargetSpace(:,1), TargetSpace(:,2), TargetSpace(:,3), ...
            opts.MarkerSize, cov, 'filled', 'MarkerFaceAlpha', 0.85);
        colormap(ax3D, cmap);
        clim(ax3D, [0 globalMaxCams]);
        cb = colorbar(ax3D);
        cb.Label.String   = 'Visible cameras';
        cb.Label.FontSize = sty.FontSizeAxis;
        cb.Label.FontName = sty.FontName;
        cb.Ticks          = 0:globalMaxCams;

        axis(ax3D, 'equal');
        grid(ax3D, 'on');
        xlabel(ax3D, 'X (m)', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
        ylabel(ax3D, 'Y (m)', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
        zlabel(ax3D, 'Z (m)', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
        view(ax3D, opts.ViewAngle);
        title(ax3D, {name, statLine}, ...
            'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName, ...
            'FontWeight', 'normal');
        set(ax3D, 'FontSize', sty.FontSizeTick, 'FontName', sty.FontName);

        % ---- Bottom row: XY worst-case (min coverage over Z) ----
        axXY = subplot(2, 2, k + 2);
        plotMinCoverageByXY(axXY, TargetSpace, cov, opts.MarkerSize, cmap, globalMaxCams);
        cb2 = colorbar(axXY);
        cb2.Label.String   = 'Min visible cameras (over Z)';
        cb2.Label.FontSize = sty.FontSizeAxis;
        cb2.Label.FontName = sty.FontName;
        cb2.Ticks          = 0:globalMaxCams;

        axis(axXY, 'equal');
        grid(axXY, 'on');
        xlabel(axXY, 'X (m)', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
        ylabel(axXY, 'Y (m)', 'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName);
        title(axXY, sprintf('%s — XY worst-case', stripCostFromName(name)), ...
            'FontSize', sty.FontSizeAxis, 'FontName', sty.FontName, ...
            'FontWeight', 'normal');
        set(axXY, 'FontSize', sty.FontSizeTick, 'FontName', sty.FontName);
    end

    sgtitle(sprintf('%s — Coverage: GA-best vs OptiTrack ad-hoc (%s, sp = %.2f m, %d cams)', ...
        ttStr, gmStr, opts.Spacing, opts.NumCameras), ...
        'FontSize', sty.FontSizeTitle, 'FontWeight', 'bold', ...
        'FontName', sty.FontName);

    applyThesisStyle(fig);

    %% Export
    if isempty(opts.SaveAs)
        outName = sprintf('CoverageHeatmap_GAvsOptiTrack_%s_%dC_GM%d_sp%.0fcm', ...
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

function [cov, numCams] = perPointVisibility(chromosome, specs)
% Count cameras that can actually see each target point.
%
% A camera "sees" the point iff ALL THREE of the following hold:
%   (a) FOV check          — the perspective projection of the point
%                            lands inside the image plane bounds
%                            [1, W] x [1, H].
%   (b) Range check        — the point is within the camera's effective
%                            range (maxCameraRange for narrow-lens,
%                            maxCameraRangeWide for wide-lens). This
%                            matches what findVisibleCameras.m enforces
%                            and what dynamicOcclusion uses.
%   (c) In-front check     — the point's depth in the camera frame is
%                            positive (z_cam > 0). The Peter Corke
%                            CentralCamera.project() can return finite,
%                            in-bounds uv coordinates for points BEHIND
%                            the camera (perspective division by
%                            negative depth flips signs); without this
%                            check the count is over-stated.
%
% NOTE: the older plotCoverageHeatmap.m and computePointUncertainty.m
% only apply check (a). Their visibility counts are therefore
% optimistic. Use this helper (or the version inside findVisibleCameras)
% as the trustworthy source.

    numCams         = specs.Cams;
    resolution      = specs.Resolution;
    focalLength     = specs.Focal;
    focalLengthWide = specs.FocalWide;
    principalPoint  = specs.PrincipalPoint;
    T               = specs.Target;
    nPts            = size(T, 1);

    maxRange     = specs.PreComputed.maxCameraRange;
    maxRangeWide = specs.PreComputed.maxCameraRangeWide;
    focalWide    = specs.FocalWide;

    [cameras, camCenters] = setupCameras(chromosome, numCams, resolution, ...
        focalLength, focalLengthWide, principalPoint, specs.PixelSize);

    % Pre-compute each camera's world-frame optical axis (the +Z column
    % of its rotation matrix). For in-front check we test the sign of
    % (point - camCenter) . opticalAxis.
    opticAxes = zeros(3, numCams);
    for c = 1:numCams
        Rcw = cameras{c}.T.rotm;          % camera-to-world rotation
        opticAxes(:, c) = Rcw(:, 3);
    end

    cov = zeros(nPts, 1);
    for pt = 1:nPts
        point    = T(pt, :);
        visCount = 0;
        for c = 1:numCams
            % --- (c) In-front-of-camera ---
            viewVec = point(:) - camCenters(:, c);   % 3x1, world frame
            depth   = dot(viewVec, opticAxes(:, c)); % z in camera frame
            if depth <= 0
                continue;
            end

            % --- (a) FOV (projection inside image plane) ---
            uv = cameras{c}.project(point);
            if ~(uv(1) >= 1 && uv(1) <= resolution(1) && ...
                 uv(2) >= 1 && uv(2) <= resolution(2))
                continue;
            end

            % --- (b) Within effective range ---
            distance = norm(viewVec);
            if cameras{c}.f == focalWide
                effRange = maxRangeWide;
            else
                effRange = maxRange;
            end
            if distance > effRange || distance <= 0
                continue;
            end

            visCount = visCount + 1;
        end
        cov(pt) = visCount;
    end
end


function plotMinCoverageByXY(ax, TargetSpace, cov, markerSize, cmap, maxCams)
% For each unique (x,y), pick the MIN coverage seen across z. Surfaces the
% zero-coverage shells that 3D occlusion otherwise hides.
    xy = TargetSpace(:, 1:2);
    [~, ~, gIdx] = unique(round(xy * 1e3) / 1e3, 'rows');  % mm-resolution
    nUnique = max(gIdx);

    minCov = nan(nUnique, 1);
    xyU    = nan(nUnique, 2);
    for k = 1:nUnique
        members  = (gIdx == k);
        minCov(k) = min(cov(members));
        xyU(k, :) = mean(xy(members, :), 1);
    end

    scatter(ax, xyU(:,1), xyU(:,2), markerSize, minCov, 'filled');
    colormap(ax, cmap);
    clim(ax, [0 maxCams]);
end


function short = stripCostFromName(name)
% Title for XY panels: keep prefix ("GA-best" / "OptiTrack ad-hoc") and
% drop the parenthetical cost so the bottom titles read clean.
    pos = strfind(name, ' (');
    if isempty(pos)
        short = name;
    else
        short = name(1:pos(1)-1);
    end
end
